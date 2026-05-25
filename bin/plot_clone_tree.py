#!/usr/bin/env python3
"""
plot_clone_tree.py

Visualise the MEDICC2 tree with clone labels overlaid.
Produces a cladogram with leaf colours representing clone membership.

Usage:
    plot_clone_tree.py \\
        --tree_pkl tree_data.pkl \\
        --assignments cell_clone_assignments.csv \\
        --out_pdf clone_tree.pdf \\
        --out_png clone_tree.png
"""

from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
from utils.tree_utils import get_leaf_names, load_pickle

logging.basicConfig(
    level=logging.INFO,
    format="[%(asctime)s] %(levelname)s: %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)
log = logging.getLogger(__name__)

# Colorblind-friendly palette (up to 12 clones; cycles if more)
PALETTE = [
    "#4C72B0", "#DD8452", "#55A868", "#C44E52", "#8172B3",
    "#937860", "#DA8BC3", "#8C8C8C", "#CCB974", "#64B5CD",
    "#1f77b4", "#ff7f0e",
]


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--tree_pkl",    required=True)
    p.add_argument("--assignments", required=True)
    p.add_argument("--out_pdf",     required=True)
    p.add_argument("--out_png",     required=True)
    return p.parse_args()


def compute_node_positions(tree, leaf_names):
    """
    Compute (x, y) positions for all nodes in a left-to-right cladogram.
    x = cumulative branch length from root
    y = evenly spaced for leaves, average of children for internal nodes
    """
    positions = {}
    leaf_order = {name: i for i, name in enumerate(leaf_names)}

    def get_y(node):
        if node.is_leaf():
            return leaf_order.get(node.taxon.label if node.taxon else "", 0)
        child_ys = [get_y(c) for c in node.child_nodes()]
        y = np.mean(child_ys)
        return y

    def set_x(node, parent_x=0.0):
        el = node.edge_length or 0.0
        x = parent_x + el
        label = node.taxon.label if node.taxon else node.label or id(node)
        positions[label] = (x, get_y(node))
        for child in node.child_nodes():
            set_x(child, x)

    set_x(tree.seed_node)
    return positions


def plot_tree(tree, assignments_df, out_pdf, out_png):
    leaf_names = get_leaf_names(tree)
    cell_to_clone = dict(zip(assignments_df["cell_id"], assignments_df["clone_id"]))

    clone_ids = sorted(assignments_df["clone_id"].unique())
    clone_color = {cid: PALETTE[i % len(PALETTE)] for i, cid in enumerate(clone_ids)}

    positions = compute_node_positions(tree, leaf_names)

    fig, ax = plt.subplots(figsize=(10, max(6, len(leaf_names) * 0.25)))

    # Draw edges
    for node in tree.preorder_node_iter():
        node_label = node.taxon.label if node.taxon else node.label or str(id(node))
        if node_label not in positions:
            continue
        nx, ny = positions[node_label]
        for child in node.child_nodes():
            clabel = child.taxon.label if child.taxon else child.label or str(id(child))
            if clabel not in positions:
                continue
            cx, cy = positions[clabel]
            # Draw L-shaped edge
            ax.plot([nx, cx], [ny, ny], color="gray", lw=0.8, zorder=1)
            ax.plot([cx, cx], [ny, cy], color="gray", lw=0.8, zorder=1)

    # Draw leaf nodes coloured by clone
    for leaf in tree.leaf_node_iter():
        label = leaf.taxon.label if leaf.taxon else ""
        if label not in positions:
            continue
        x, y = positions[label]
        clone = cell_to_clone.get(label, "unassigned")
        color = clone_color.get(clone, "#aaaaaa")
        ax.scatter(x, y, c=color, s=20, zorder=3, linewidths=0)

    ax.set_xlabel("Branch length (MEDICC2 units)")
    ax.set_ylabel("Cells")
    ax.set_yticks([])
    ax.set_title("MEDICC2 tree — clone assignment")

    legend_patches = [
        mpatches.Patch(color=clone_color[cid], label=cid)
        for cid in clone_ids
    ]
    ax.legend(handles=legend_patches, loc="lower right",
              title="Clone", fontsize=7, title_fontsize=8)

    plt.tight_layout()
    fig.savefig(out_pdf, bbox_inches="tight")
    fig.savefig(out_png, dpi=150, bbox_inches="tight")
    plt.close(fig)
    log.info(f"Tree plots saved to {out_pdf} and {out_png}")


def main() -> None:
    args = parse_args()
    tree = load_pickle(args.tree_pkl)
    assignments = pd.read_csv(args.assignments)
    plot_tree(tree, assignments, args.out_pdf, args.out_png)


if __name__ == "__main__":
    main()
