#!/usr/bin/env python3
"""
assign_clones.py

Assigns cells to clones using the MEDICC2 phylogenetic tree and scUnique events.

Strategies
----------
internal_node   Cut the tree at meaningful internal nodes (long branches,
                sufficient subtree size). This is the default and most
                biologically interpretable strategy.

distance        Hierarchical clustering on the tree pairwise distance matrix.
                Requires --distance_threshold.

event_profile   Cluster cells by Jaccard similarity of scUnique event profiles
                (which copy-number events each cell carries). Pure phenotypic
                clustering, tree-agnostic.

hybrid          Use tree topology to define initial clone groups, then
                refine by event similarity within each group (split if
                within-group Jaccard similarity drops below threshold).

Output
------
cell_clone_assignments.csv  : cell_id, clone_id, strategy, confidence
clone_summary.csv           : clone_id, n_cells, defining_events, mean_confidence
clone_events.csv            : clone_id, event (defining copy-number events)
clone_tree_annotated.new    : Newick tree with clone labels at internal nodes

Usage
-----
assign_clones.py --tree_pkl tree_data.pkl --node_table node_table.csv \\
    --scunique_events events.tsv --cell_metadata cell_metadata.csv \\
    --bam_manifest bam_manifest.csv \\
    --strategy internal_node --min_cells_per_clone 5 --min_branch_length 0.0 \\
    --event_similarity_thr 0.5 --small_clone_action merge \\
    --out_assignments cell_clone_assignments.csv \\
    --out_summary clone_summary.csv \\
    --out_events clone_events.csv \\
    --out_tree clone_tree_annotated.new
"""

from __future__ import annotations

import argparse
import json
import logging
import sys
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple

import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import fcluster, linkage
from scipy.spatial.distance import pdist, squareform

sys.path.insert(0, str(Path(__file__).parent))
from utils.tree_utils import (
    build_node_table,
    cut_tree_at_nodes,
    get_clade_leaves,
    get_leaf_names,
    get_node_by_label,
    load_pickle,
    pairwise_distance_matrix,
    save_newick,
    select_internal_nodes_by_branch_length,
)

logging.basicConfig(
    level=logging.INFO,
    format="[%(asctime)s] %(levelname)s: %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)
log = logging.getLogger(__name__)


# ─── CLI ─────────────────────────────────────────────────────────────────────

def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--tree_pkl",             required=True)
    p.add_argument("--node_table",           required=True)
    p.add_argument("--scunique_events",      required=True)
    p.add_argument("--cell_metadata",        required=True)
    p.add_argument("--bam_manifest",         required=True)
    p.add_argument("--strategy",             default="internal_node",
                   choices=["internal_node", "distance", "event_profile", "hybrid"])
    p.add_argument("--min_cells_per_clone",  type=int,   default=5)
    p.add_argument("--min_branch_length",    type=float, default=0.0)
    p.add_argument("--distance_threshold",   type=float, default=None)
    p.add_argument("--event_similarity_thr", type=float, default=0.5)
    p.add_argument("--max_clones",           type=int,   default=None)
    p.add_argument("--small_clone_action",   default="merge",
                   choices=["drop", "merge", "flag"])
    p.add_argument("--out_assignments",      required=True)
    p.add_argument("--out_summary",          required=True)
    p.add_argument("--out_events",           required=True)
    p.add_argument("--out_tree",             required=True)
    return p.parse_args()


# ─── Event profile helpers ───────────────────────────────────────────────────

def load_scunique_events(path: str) -> pd.DataFrame:
    """
    Load scUnique per-cell unique/recent events.
    Expected columns (adjust if your scUnique version differs):
        cell_id, chr, start, end, event_type  [gain | loss | neutral]

    TODO: Verify column names match your scUnique output format.
    """
    df = pd.read_csv(path, sep="\t")
    required = {"cell_id"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"scUnique events file missing columns: {missing}")
    log.info(f"Loaded {len(df)} scUnique events across {df['cell_id'].nunique()} cells")
    return df


def build_event_profile_matrix(
    events_df: pd.DataFrame,
    all_cells: List[str],
) -> pd.DataFrame:
    """
    Build a binary (cell × event) matrix.
    Each column is a unique (chr, start, end, event_type) combination.
    Cells with no events in a column get 0.

    Returns a DataFrame indexed by cell_id.
    """
    # Create a string key for each event
    key_cols = [c for c in ["chr", "chrom", "chromosome", "start", "end", "event_type"]
                if c in events_df.columns]
    if not key_cols:
        log.warning("Could not identify event coordinate columns; using first 3 columns")
        key_cols = events_df.columns.drop("cell_id").tolist()[:3]

    events_df = events_df.copy()
    events_df["event_key"] = events_df[key_cols].astype(str).agg("|".join, axis=1)

    # Pivot: rows = cell_id, cols = event_key, values = 1/0
    profile = (
        events_df.groupby(["cell_id", "event_key"])
        .size()
        .unstack(fill_value=0)
        .clip(upper=1)
    )

    # Add cells with zero events
    missing_cells = set(all_cells) - set(profile.index)
    if missing_cells:
        empty_rows = pd.DataFrame(
            0,
            index=list(missing_cells),
            columns=profile.columns,
        )
        profile = pd.concat([profile, empty_rows])

    profile = profile.loc[all_cells]  # reindex to canonical order
    log.info(
        f"Event profile matrix: {profile.shape[0]} cells × {profile.shape[1]} events"
    )
    return profile


def jaccard_distance_matrix(profile: pd.DataFrame) -> np.ndarray:
    """Pairwise Jaccard distances between cell event profiles."""
    return squareform(pdist(profile.values, metric="jaccard"))


# ─── Clone assignment strategies ────────────────────────────────────────────

def strategy_internal_node(
    tree,
    node_table: pd.DataFrame,
    min_cells: int,
    min_branch_length: float,
    max_clones: Optional[int],
) -> Dict[str, List[str]]:
    """
    Cut the tree at internal nodes with long incoming branches and
    sufficient subtree sizes.

    Biological rationale: Long branches in a MEDICC2 tree indicate many
    copy-number changes accumulated on that lineage — these define clonal
    populations. Shorter branches are more likely noise or rare transitions.
    """
    log.info(
        f"[internal_node] min_branch_length={min_branch_length}, "
        f"min_cells={min_cells}"
    )

    cut_nodes = select_internal_nodes_by_branch_length(
        tree, min_branch_length, min_cells
    )
    log.info(f"Selected {len(cut_nodes)} internal nodes as clone roots")

    if max_clones and len(cut_nodes) > max_clones:
        # Keep the nodes with the longest branches (most distinct clades)
        lengths = {}
        for label in cut_nodes:
            node = get_node_by_label(tree, label)
            lengths[label] = node.edge_length or 0.0
        cut_nodes = sorted(cut_nodes, key=lambda l: -lengths[l])[:max_clones]
        log.info(f"Capped to {max_clones} clones by branch length")

    if not cut_nodes:
        log.warning(
            "No internal nodes met criteria — assigning all cells to a single clone"
        )
        return {"clone_1": get_leaf_names(tree)}

    return cut_tree_at_nodes(tree, cut_nodes)


def strategy_distance(
    tree,
    distance_threshold: float,
    min_cells: int,
) -> Dict[str, List[str]]:
    """
    Hierarchical clustering on pairwise tree distances.
    Uses Ward linkage; cuts at *distance_threshold*.
    """
    if distance_threshold is None:
        raise ValueError("--distance_threshold is required for the distance strategy")
    log.info(f"[distance] threshold={distance_threshold}")

    dist_df = pairwise_distance_matrix(tree)
    cells = dist_df.index.tolist()
    condensed = squareform(dist_df.values)
    Z = linkage(condensed, method="ward")
    labels = fcluster(Z, t=distance_threshold, criterion="distance")

    clones: Dict[str, List[str]] = {}
    for cell, label in zip(cells, labels):
        clone_key = f"clone_{label}"
        clones.setdefault(clone_key, []).append(cell)
    return clones


def strategy_event_profile(
    all_cells: List[str],
    events_df: pd.DataFrame,
    event_similarity_thr: float,
    min_cells: int,
    max_clones: Optional[int],
) -> Dict[str, List[str]]:
    """
    Cluster cells by Jaccard similarity of their scUnique event profiles.
    Uses average linkage; cuts so within-cluster similarity ≥ event_similarity_thr.
    """
    log.info(f"[event_profile] similarity_threshold={event_similarity_thr}")
    profile = build_event_profile_matrix(events_df, all_cells)

    if profile.shape[1] == 0:
        log.warning("No events found — all cells → single clone")
        return {"clone_1": all_cells}

    dist_mat = jaccard_distance_matrix(profile)
    condensed = squareform(dist_mat)
    Z = linkage(condensed, method="average")
    t = 1.0 - event_similarity_thr  # distance = 1 - similarity
    labels = fcluster(Z, t=t, criterion="distance")

    clones: Dict[str, List[str]] = {}
    for cell, label in zip(all_cells, labels):
        clones.setdefault(f"clone_{label}", []).append(cell)

    if max_clones and len(clones) > max_clones:
        # Merge smallest clones until target reached
        clones = _merge_to_n(clones, max_clones)
    return clones


def strategy_hybrid(
    tree,
    all_cells: List[str],
    events_df: pd.DataFrame,
    min_cells: int,
    min_branch_length: float,
    event_similarity_thr: float,
    max_clones: Optional[int],
) -> Dict[str, List[str]]:
    """
    Two-stage hybrid:
      1. Internal-node grouping defines initial clones (tree-based).
      2. Each clone is checked for internal event homogeneity.
         If within-clone Jaccard similarity < threshold, the clone is split.
    """
    log.info("[hybrid] stage-1: internal_node grouping")
    initial_clones = strategy_internal_node(
        tree, None, min_cells=2, min_branch_length=min_branch_length,
        max_clones=max_clones,
    )

    profile = build_event_profile_matrix(events_df, all_cells)

    final_clones: Dict[str, List[str]] = {}
    split_counter = 0
    for clone_id, cells in initial_clones.items():
        if len(cells) < 2 or profile.shape[1] == 0:
            final_clones[clone_id] = cells
            continue

        sub = profile.loc[cells]
        dist_mat = jaccard_distance_matrix(sub)
        mean_sim = 1 - dist_mat.mean()

        if mean_sim >= event_similarity_thr or len(cells) < min_cells:
            final_clones[clone_id] = cells
        else:
            # Sub-cluster within this initial clone
            condensed = squareform(dist_mat)
            Z = linkage(condensed, method="average")
            t = 1.0 - event_similarity_thr
            sub_labels = fcluster(Z, t=t, criterion="distance")
            for cell, sub_label in zip(cells, sub_labels):
                sub_clone = f"{clone_id}_sub{sub_label}"
                final_clones.setdefault(sub_clone, []).append(cell)
            split_counter += 1
            log.info(f"  Split {clone_id} (mean_sim={mean_sim:.3f}) → {len(set(sub_labels))} sub-clones")

    log.info(f"[hybrid] {len(initial_clones)} initial → {len(final_clones)} final clones (split {split_counter})")
    return final_clones


# ─── Small-clone handling ────────────────────────────────────────────────────

def _merge_to_n(clones: Dict[str, List[str]], n: int) -> Dict[str, List[str]]:
    """Merge the two smallest clones repeatedly until ≤n clones remain."""
    while len(clones) > n:
        sorted_keys = sorted(clones, key=lambda k: len(clones[k]))
        smallest, second = sorted_keys[0], sorted_keys[1]
        clones[second] = clones[second] + clones.pop(smallest)
    return clones


def handle_small_clones(
    clones: Dict[str, List[str]],
    min_cells: int,
    action: str,
) -> Tuple[Dict[str, List[str]], List[str]]:
    """
    Handle clones smaller than *min_cells* according to *action*:
      drop  — remove those cells entirely
      merge — merge them into the largest clone
      flag  — keep but tag with a warning

    Returns (processed_clones, list_of_flagged_clone_ids).
    """
    small  = {k: v for k, v in clones.items() if len(v) < min_cells}
    normal = {k: v for k, v in clones.items() if len(v) >= min_cells}
    flags: List[str] = []

    if not small:
        return clones, flags

    log.warning(
        f"{len(small)} clone(s) below min_cells={min_cells}: "
        + ", ".join(f"{k}({len(v)})" for k, v in small.items())
    )

    if action == "drop":
        total_dropped = sum(len(v) for v in small.values())
        log.warning(f"Dropping {total_dropped} cells from small clones")
        return normal, flags

    elif action == "merge":
        if not normal:
            log.warning("All clones are small — keeping as single clone")
            merged_cells = [c for cells in clones.values() for c in cells]
            return {"clone_merged": merged_cells}, flags
        # Merge small cells into the largest clone
        largest = max(normal, key=lambda k: len(normal[k]))
        for cells in small.values():
            normal[largest].extend(cells)
        log.info(f"Merged {len(small)} small clone(s) into {largest}")
        return normal, flags

    elif action == "flag":
        flags = list(small.keys())
        return clones, flags

    return clones, flags


# ─── Output builders ─────────────────────────────────────────────────────────

def build_assignment_table(
    clones: Dict[str, List[str]],
    flagged: List[str],
    strategy: str,
) -> pd.DataFrame:
    rows = []
    for clone_id, cells in clones.items():
        for cell_id in cells:
            rows.append(
                {
                    "cell_id":    cell_id,
                    "clone_id":   clone_id,
                    "strategy":   strategy,
                    "flagged":    clone_id in flagged,
                    "confidence": 1.0,  # TODO: replace with soft assignment probability
                }
            )
    return pd.DataFrame(rows).sort_values(["clone_id", "cell_id"])


def build_clone_summary(
    clones: Dict[str, List[str]],
    events_df: pd.DataFrame,
) -> pd.DataFrame:
    rows = []
    for clone_id, cells in clones.items():
        # Find events shared by ≥50% of cells in this clone
        clone_events = events_df[events_df["cell_id"].isin(cells)]
        if "event_key" not in clone_events.columns and not clone_events.empty:
            key_cols = [c for c in ["chr", "chrom", "start", "end", "event_type"]
                        if c in clone_events.columns]
            clone_events = clone_events.copy()
            clone_events["event_key"] = clone_events[key_cols].astype(str).agg("|".join, axis=1)
        if not clone_events.empty and "event_key" in clone_events.columns:
            event_freq = clone_events.groupby("event_key")["cell_id"].nunique() / len(cells)
            defining = event_freq[event_freq >= 0.5].index.tolist()
        else:
            defining = []
        rows.append(
            {
                "clone_id":       clone_id,
                "n_cells":        len(cells),
                "n_defining_events": len(defining),
                "defining_events": ";".join(defining[:10]),  # cap for readability
            }
        )
    return pd.DataFrame(rows).sort_values("n_cells", ascending=False)


def build_clone_events_table(
    clones: Dict[str, List[str]],
    events_df: pd.DataFrame,
) -> pd.DataFrame:
    """Long-form table: one row per (clone_id, defining_event)."""
    rows = []
    key_cols = [c for c in ["chr", "chrom", "start", "end", "event_type"]
                if c in events_df.columns]
    for clone_id, cells in clones.items():
        clone_evts = events_df[events_df["cell_id"].isin(cells)].copy()
        if clone_evts.empty or not key_cols:
            continue
        clone_evts["event_key"] = clone_evts[key_cols].astype(str).agg("|".join, axis=1)
        freq = (
            clone_evts.groupby("event_key")["cell_id"].nunique() / len(cells)
        ).rename("fraction_cells")
        for ek, frac in freq.items():
            rows.append({"clone_id": clone_id, "event": ek, "fraction_cells": frac})
    return pd.DataFrame(rows)


def annotate_tree_with_clones(tree, clones: Dict[str, List[str]]) -> None:
    """Attach clone labels to internal tree nodes in-place."""
    cell_to_clone = {
        cell: clone_id
        for clone_id, cells in clones.items()
        for cell in cells
    }
    from utils.tree_utils import get_clade_leaves
    for node in tree.preorder_node_iter():
        if node.is_leaf():
            node.label = cell_to_clone.get(
                node.taxon.label if node.taxon else "", ""
            )
        else:
            clade_cells = get_clade_leaves(node)
            clone_ids = {cell_to_clone.get(c, "") for c in clade_cells}
            clone_ids.discard("")
            node.label = ",".join(sorted(clone_ids)) if clone_ids else ""


# ─── Main ────────────────────────────────────────────────────────────────────

def main() -> None:
    args = parse_args()

    # Load tree
    log.info(f"Loading tree from {args.tree_pkl}")
    tree = load_pickle(args.tree_pkl)
    all_cells = get_leaf_names(tree)
    log.info(f"Tree has {len(all_cells)} leaf cells")

    # Load supplementary data
    node_table   = pd.read_csv(args.node_table)
    events_df    = load_scunique_events(args.scunique_events)
    cell_meta    = pd.read_csv(args.cell_metadata)
    bam_manifest = pd.read_csv(args.bam_manifest)

    # Validate cell overlap
    bam_cells = set(bam_manifest["cell_id"].astype(str))
    tree_cells = set(all_cells)
    in_both = tree_cells & bam_cells
    log.info(
        f"Cells: {len(tree_cells)} in tree, {len(bam_cells)} in BAM manifest, "
        f"{len(in_both)} in both"
    )
    if len(in_both) < 2:
        log.error("Fewer than 2 cells overlap between tree and BAM manifest — aborting")
        sys.exit(1)
    ordered_cells = [cell for cell in all_cells if cell in in_both]

    # Run chosen clone assignment strategy
    if args.strategy == "internal_node":
        clones = strategy_internal_node(
            tree, node_table,
            min_cells=args.min_cells_per_clone,
            min_branch_length=args.min_branch_length,
            max_clones=args.max_clones,
        )
    elif args.strategy == "distance":
        clones = strategy_distance(
            tree,
            distance_threshold=args.distance_threshold,
            min_cells=args.min_cells_per_clone,
        )
    elif args.strategy == "event_profile":
        clones = strategy_event_profile(
            all_cells=ordered_cells,
            events_df=events_df,
            event_similarity_thr=args.event_similarity_thr,
            min_cells=args.min_cells_per_clone,
            max_clones=args.max_clones,
        )
    elif args.strategy == "hybrid":
        clones = strategy_hybrid(
            tree=tree,
            all_cells=ordered_cells,
            events_df=events_df,
            min_cells=args.min_cells_per_clone,
            min_branch_length=args.min_branch_length,
            event_similarity_thr=args.event_similarity_thr,
            max_clones=args.max_clones,
        )
    else:
        raise ValueError(f"Unknown strategy: {args.strategy}")

    # Handle small clones
    clones, flagged = handle_small_clones(
        clones,
        min_cells=args.min_cells_per_clone,
        action=args.small_clone_action,
    )

    # Rename clones sequentially for cleanliness
    clones = {f"clone_{i+1:03d}": cells
              for i, (_, cells) in enumerate(sorted(clones.items()))}

    log.info(f"Final: {len(clones)} clones")
    for cid, cells in clones.items():
        log.info(f"  {cid}: {len(cells)} cells")

    # Build output tables
    assignments = build_assignment_table(clones, flagged, args.strategy)
    summary     = build_clone_summary(clones, events_df)
    evt_table   = build_clone_events_table(clones, events_df)

    # Annotate tree
    annotate_tree_with_clones(tree, clones)
    save_newick(tree, args.out_tree)

    # Write outputs
    assignments.to_csv(args.out_assignments, index=False)
    summary.to_csv(args.out_summary, index=False)
    evt_table.to_csv(args.out_events, index=False)

    log.info(f"Assignments → {args.out_assignments}")
    log.info(f"Summary     → {args.out_summary}")
    log.info(f"Events      → {args.out_events}")
    log.info(f"Annotated tree → {args.out_tree}")


if __name__ == "__main__":
    main()
