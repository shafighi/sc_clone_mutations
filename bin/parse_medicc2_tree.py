#!/usr/bin/env python3
"""
parse_medicc2_tree.py

Parse a MEDICC2 Newick tree and produce:
  - Pickled dendropy.Tree object
  - Node metadata table (node_table.csv)
  - Edge table (edge_table.csv)
  - Pairwise cell distance matrix (pairwise_dist.csv)

If a MEDICC2 events-per-branch TSV is provided, branch events are annotated
onto the node table.

Usage:
    parse_medicc2_tree.py --tree_file tree.new --out_pkl tree_data.pkl \\
        --out_nodes node_table.csv --out_edges edge_table.csv \\
        --out_distances pairwise_dist.csv [--events_tsv events.tsv]
"""

import argparse
import logging
import sys
from pathlib import Path

import pandas as pd

# Add bin/ to path (Nextflow copies scripts into work dir)
sys.path.insert(0, str(Path(__file__).parent))
from utils.tree_utils import (
    build_edge_table,
    build_node_table,
    load_newick,
    pairwise_distance_matrix,
    save_pickle,
)

logging.basicConfig(
    level=logging.INFO,
    format="[%(asctime)s] %(levelname)s: %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)
log = logging.getLogger(__name__)


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--tree_file",    required=True, help="MEDICC2 Newick tree file")
    p.add_argument("--events_tsv",   default=None,  help="MEDICC2 events-per-branch TSV (optional)")
    p.add_argument("--out_pkl",      required=True, help="Output: pickled dendropy tree")
    p.add_argument("--out_nodes",    required=True, help="Output: node metadata CSV")
    p.add_argument("--out_edges",    required=True, help="Output: edge table CSV")
    p.add_argument("--out_distances",required=True, help="Output: pairwise distance CSV")
    return p.parse_args()


def load_medicc2_events(events_tsv: str) -> pd.DataFrame:
    """
    Load MEDICC2 branch events TSV.
    Expected columns: child_node, parent_node, chrom, start, end, cn_change, event_type
    TODO: adjust column names to match your MEDICC2 version output.
    """
    df = pd.read_csv(events_tsv, sep="\t")
    log.info(f"Loaded {len(df)} branch events from {events_tsv}")
    return df


def annotate_nodes_with_events(
    node_table: pd.DataFrame,
    events: pd.DataFrame,
) -> pd.DataFrame:
    """
    Add event summary columns to node_table by joining on node_id.
    """
    if events is None or events.empty:
        node_table["n_events"] = 0
        node_table["event_summary"] = ""
        return node_table

    # TODO: adjust the column name used for the child node in your MEDICC2 version
    child_col = "child_node" if "child_node" in events.columns else events.columns[0]
    event_counts = events.groupby(child_col).size().rename("n_events")
    node_table = node_table.merge(
        event_counts,
        left_on="node_id",
        right_index=True,
        how="left",
    )
    node_table["n_events"] = node_table["n_events"].fillna(0).astype(int)
    return node_table


def main() -> None:
    args = parse_args()

    # Load tree
    log.info(f"Loading MEDICC2 tree from {args.tree_file}")
    tree = load_newick(args.tree_file)

    n_leaves   = sum(1 for _ in tree.leaf_node_iter())
    n_internal = sum(1 for n in tree.preorder_node_iter() if not n.is_leaf())
    log.info(f"Tree: {n_leaves} leaves, {n_internal} internal nodes")

    if n_leaves < 2:
        log.error("Tree has fewer than 2 leaf cells — check input file format")
        sys.exit(1)

    # Build node and edge tables
    node_table = build_node_table(tree)
    edge_table = build_edge_table(tree)

    # Annotate with MEDICC2 branch events if provided
    events_df = None
    if args.events_tsv:
        log.info(f"Loading branch events from {args.events_tsv}")
        events_df = load_medicc2_events(args.events_tsv)
        node_table = annotate_nodes_with_events(node_table, events_df)

    # Compute pairwise distances
    log.info("Computing pairwise leaf distances …")
    dist_matrix = pairwise_distance_matrix(tree)

    # Save outputs
    save_pickle(tree, args.out_pkl)
    log.info(f"Saved pickled tree → {args.out_pkl}")

    node_table.to_csv(args.out_nodes, index=False)
    log.info(f"Saved node table ({len(node_table)} rows) → {args.out_nodes}")

    edge_table.to_csv(args.out_edges, index=False)
    log.info(f"Saved edge table ({len(edge_table)} rows) → {args.out_edges}")

    dist_matrix.to_csv(args.out_distances)
    log.info(f"Saved {n_leaves}×{n_leaves} distance matrix → {args.out_distances}")


if __name__ == "__main__":
    main()
