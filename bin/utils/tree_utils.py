"""
tree_utils.py — helper functions for MEDICC2 tree manipulation.

All tree objects use dendropy.Tree internally. Functions accept either
a dendropy.Tree or a path to a Newick file.
"""

from __future__ import annotations

import pickle
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple

import dendropy
import numpy as np
import pandas as pd


# ─── I/O ────────────────────────────────────────────────────────────────────

def load_newick(path: str | Path) -> dendropy.Tree:
    """Load a Newick tree file, tolerating MEDICC2's label conventions."""
    tree = dendropy.Tree.get(
        path=str(path),
        schema="newick",
        preserve_underscores=True,
        rooting="force-rooted",
    )
    _relabel_internal_nodes(tree)
    return tree


def save_pickle(tree: dendropy.Tree, path: str | Path) -> None:
    with open(path, "wb") as fh:
        pickle.dump(tree, fh, protocol=pickle.HIGHEST_PROTOCOL)


def load_pickle(path: str | Path) -> dendropy.Tree:
    with open(path, "rb") as fh:
        return pickle.load(fh)


def save_newick(tree: dendropy.Tree, path: str | Path) -> None:
    tree.write(path=str(path), schema="newick", suppress_rooting=False)


# ─── Tree annotation helpers ─────────────────────────────────────────────────

def _relabel_internal_nodes(tree: dendropy.Tree) -> None:
    """Give unnamed internal nodes auto-generated labels (n0, n1, ...)."""
    counter = 0
    for node in tree.preorder_node_iter():
        if node.taxon is None and not node.is_leaf():
            node.label = f"internal_{counter}"
            counter += 1
        elif node.taxon is not None:
            node.label = node.taxon.label


def get_leaf_names(tree: dendropy.Tree) -> List[str]:
    return [leaf.taxon.label for leaf in tree.leaf_node_iter()]


def get_node_by_label(tree: dendropy.Tree, label: str) -> Optional[dendropy.Node]:
    for node in tree.preorder_node_iter():
        if (node.taxon and node.taxon.label == label) or node.label == label:
            return node
    return None


# ─── Subtree operations ──────────────────────────────────────────────────────

def get_clade_leaves(node: dendropy.Node) -> List[str]:
    """Return all leaf labels in the subtree rooted at *node*."""
    return [leaf.taxon.label for leaf in node.leaf_iter()]


def count_clade_leaves(node: dendropy.Node) -> int:
    return sum(1 for _ in node.leaf_iter())


def subtree_root_for_cells(
    tree: dendropy.Tree, cell_ids: Set[str]
) -> Optional[dendropy.Node]:
    """Find the MRCA (most recent common ancestor) of a set of cells."""
    if not cell_ids:
        return None
    taxa = [tree.taxon_namespace.find_taxon(label=c) for c in cell_ids]
    taxa = [t for t in taxa if t is not None]
    if not taxa:
        return None
    return tree.mrca(taxa=taxa)


# ─── Distance computation ────────────────────────────────────────────────────

def pairwise_distance_matrix(tree: dendropy.Tree) -> pd.DataFrame:
    """
    Compute the phylogenetic distance between every pair of leaf nodes.
    Returns a symmetric DataFrame indexed by cell labels.
    """
    pdm = tree.phylogenetic_distance_matrix()
    leaves = list(tree.leaf_node_iter())
    labels = [leaf.taxon.label for leaf in leaves]
    n = len(labels)
    mat = np.zeros((n, n))
    for i, li in enumerate(leaves):
        for j, lj in enumerate(leaves):
            if i < j:
                d = pdm(li.taxon, lj.taxon)
                mat[i, j] = d
                mat[j, i] = d
    return pd.DataFrame(mat, index=labels, columns=labels)


# ─── Node / edge tables ──────────────────────────────────────────────────────

def build_node_table(tree: dendropy.Tree) -> pd.DataFrame:
    """Return a DataFrame with one row per tree node."""
    rows = []
    for node in tree.preorder_node_iter():
        label = node.taxon.label if node.taxon else node.label or ""
        parent = None
        if node.parent_node is not None:
            pn = node.parent_node
            parent = pn.taxon.label if pn.taxon else pn.label or ""
        rows.append(
            {
                "node_id":      label,
                "parent_id":    parent,
                "is_leaf":      node.is_leaf(),
                "edge_length":  node.edge_length if node.edge_length else 0.0,
                "n_leaves":     count_clade_leaves(node),
            }
        )
    return pd.DataFrame(rows)


def build_edge_table(tree: dendropy.Tree) -> pd.DataFrame:
    """Return a DataFrame with one row per tree edge (child → parent)."""
    rows = []
    for node in tree.preorder_node_iter():
        if node.parent_node is None:
            continue
        child_label  = node.taxon.label if node.taxon else node.label or ""
        parent_label = (node.parent_node.taxon.label
                        if node.parent_node.taxon
                        else node.parent_node.label or "")
        rows.append(
            {
                "child":       child_label,
                "parent":      parent_label,
                "edge_length": node.edge_length if node.edge_length else 0.0,
            }
        )
    return pd.DataFrame(rows)


# ─── Clone-cut helpers ───────────────────────────────────────────────────────

def cut_tree_into_clones(
    tree: dendropy.Tree,
    min_cells: int,
    min_branch_length: float = 0.0,
) -> Dict[str, List[str]]:
    """
    Partition a tree into clones by recursively cutting the longest qualifying
    branch, but only where BOTH resulting parts contain at least *min_cells*
    leaves. This guarantees every emitted clone is adequately sized, so small
    fragments never collapse back into a single clone.

    A branch qualifies if its length is >= *min_branch_length*. When
    *min_branch_length* <= 0 (or None), an automatic threshold of
    mean + 1 standard deviation of internal edge lengths is used, so only the
    longest branches (clear clone boundaries) are cut.

    Returns {clade_label: [cell_id, ...]}.
    """
    internal_edges = [
        node.edge_length or 0.0
        for node in tree.preorder_node_iter()
        if not node.is_leaf() and node.parent_node is not None
    ]
    threshold = min_branch_length
    if threshold is None or threshold <= 0.0:
        if internal_edges:
            arr = np.array(internal_edges, dtype=float)
            threshold = float(arr.mean() + arr.std())
        else:
            threshold = 0.0

    # Precompute the leaf set under every node (postorder, so children first).
    node_clade: Dict[dendropy.Node, frozenset] = {}
    for node in tree.postorder_node_iter():
        if node.is_leaf():
            node_clade[node] = frozenset([node.taxon.label])
        else:
            s: Set[str] = set()
            for ch in node.child_nodes():
                s |= node_clade[ch]
            node_clade[node] = frozenset(s)

    # Internal nodes whose incoming edge is long enough to be a cut candidate.
    candidates = [
        n for n in tree.preorder_node_iter()
        if not n.is_leaf() and n.parent_node is not None
        and (n.edge_length or 0.0) >= threshold
    ]

    clones: List[Set[str]] = []

    def partition(leafset: Set[str]) -> None:
        if len(leafset) < 2 * min_cells:
            clones.append(set(leafset))
            return
        best_node = None
        best_el = -1.0
        best_lower: Optional[Set[str]] = None
        for node in candidates:
            lower = node_clade[node] & leafset
            li = len(lower)
            if li >= min_cells and (len(leafset) - li) >= min_cells:
                el = node.edge_length or 0.0
                if el > best_el:
                    best_el, best_node, best_lower = el, node, set(lower)
        if best_node is None:
            clones.append(set(leafset))
            return
        partition(best_lower)
        partition(set(leafset) - best_lower)

    partition(set(get_leaf_names(tree)))
    return {f"clade_{i+1}": sorted(c) for i, c in enumerate(clones)}


def cut_tree_at_nodes(
    tree: dendropy.Tree,
    cut_nodes: List[str],
) -> Dict[str, List[str]]:
    """
    Given a set of internal node labels where the tree is cut,
    return {clone_label: [cell_id, ...]} mapping.

    Cells not under any cut node are assigned to their nearest ancestor cut node,
    or to 'root_clade' if no cut node is above them.
    """
    clones: Dict[str, List[str]] = {}
    assigned: Set[str] = set()

    # Resolve nodes and sort by clade size ascending so the deepest (most
    # specific) selected clades claim their cells before their ancestors.
    resolved = [(lbl, get_node_by_label(tree, lbl)) for lbl in cut_nodes]
    resolved = [(lbl, n) for lbl, n in resolved if n is not None]
    resolved.sort(key=lambda x: count_clade_leaves(x[1]))

    for node_label, node in resolved:
        leaves = get_clade_leaves(node)
        # Don't double-assign already-claimed cells
        new_leaves = [l for l in leaves if l not in assigned]
        if new_leaves:
            clones[node_label] = new_leaves
            assigned.update(new_leaves)

    # Assign remaining cells to a catch-all clone
    all_leaves = set(get_leaf_names(tree))
    unassigned = all_leaves - assigned
    if unassigned:
        clones["unassigned"] = list(unassigned)

    return clones


def select_internal_nodes_by_branch_length(
    tree: dendropy.Tree,
    min_branch_length: float,
    min_cells: int,
) -> List[str]:
    """
    Select internal nodes whose incoming edge marks a clone boundary: the edge
    is long (many accumulated copy-number changes) and the subtree holds at
    least *min_cells* leaves. The root is never selected (it has no incoming
    edge), which prevents the whole tree collapsing into a single clone.

    If *min_branch_length* <= 0 (or None), an automatic threshold is derived
    from the distribution of internal edge lengths (mean + 1 standard
    deviation), so only the longest branches — the clearest clone boundaries —
    are cut.
    """
    internal_edges = [
        node.edge_length or 0.0
        for node in tree.preorder_node_iter()
        if not node.is_leaf() and node.parent_node is not None
    ]

    threshold = min_branch_length
    if threshold is None or threshold <= 0.0:
        if internal_edges:
            arr = np.array(internal_edges, dtype=float)
            threshold = float(arr.mean() + arr.std())
        else:
            threshold = 0.0

    selected = []
    for node in tree.preorder_node_iter():
        if node.is_leaf() or node.parent_node is None:
            continue  # skip leaves and the root
        el = node.edge_length or 0.0
        if el > threshold and count_clade_leaves(node) >= min_cells:
            selected.append(node.label or (node.taxon.label if node.taxon else ""))
    return selected
