#!/usr/bin/env python3
"""
compare_variants.py

Compare somatic variants across clones and callers.

Produces:
  - variant_matrix.csv      : variant × (clone, caller) presence/absence
  - caller_concordance.csv  : per-clone concordance statistics
  - private_shared_table.csv: variant classification (private/partial/shared)
  - per_clone_vcfs/         : per-clone merged VCFs for each caller
  - comparison_plots/       : UpSet plot, heatmap, Venn diagram (2-4 callers)

Input: all normalized VCF files (passed as a list on the command line).
       VCF filenames must follow the pattern: {clone_id}.{caller}.norm.vcf.gz

Usage:
    compare_variants.py \\
        --vcf_list clone_001.mutect2.norm.vcf.gz clone_001.strelka2.norm.vcf.gz ... \\
        --clone_summary clone_summary.csv \\
        --out_matrix variant_matrix.csv \\
        --out_concordance caller_concordance.csv \\
        --out_private private_shared_table.csv \\
        --out_per_clone per_clone_vcfs/ \\
        --out_plots comparison_plots/ \\
        --pass_only True
"""

from __future__ import annotations

import argparse
import logging
import re
import sys
from pathlib import Path
from typing import Dict, List, Tuple

import matplotlib
matplotlib.use("Agg")  # no display
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

sys.path.insert(0, str(Path(__file__).parent))
from utils.vcf_utils import (
    build_presence_absence_matrix,
    caller_concordance_per_clone,
    classify_variant_sharing,
    vcf_to_dataframe,
    vcf_to_variant_set,
)

logging.basicConfig(
    level=logging.INFO,
    format="[%(asctime)s] %(levelname)s: %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)
log = logging.getLogger(__name__)


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--vcf_list",        nargs="+", required=True)
    p.add_argument("--clone_summary",   required=True)
    p.add_argument("--out_matrix",      required=True)
    p.add_argument("--out_concordance", required=True)
    p.add_argument("--out_private",     required=True)
    p.add_argument("--out_per_clone",   required=True)
    p.add_argument("--out_plots",       required=True)
    p.add_argument("--pass_only",       type=lambda x: x.lower() == "true", default=True)
    return p.parse_args()


def parse_vcf_filename(path: str) -> Tuple[str, str]:
    """
    Extract clone_id and caller from filename convention:
      {clone_id}.{caller}.norm.vcf.gz
    """
    name = Path(path).name
    # Remove suffixes
    name = re.sub(r"\.norm\.vcf(\.gz)?$", "", name)
    name = re.sub(r"\.vcf(\.gz)?$", "", name)
    parts = name.split(".")
    if len(parts) < 2:
        raise ValueError(f"Cannot parse clone/caller from filename: {path}")
    caller   = parts[-1]
    clone_id = ".".join(parts[:-1])
    return clone_id, caller


def plot_concordance_heatmap(
    concordance_df: pd.DataFrame,
    out_dir: str,
) -> None:
    """Bar chart of per-clone concordance rates."""
    if concordance_df.empty:
        return
    fig, ax = plt.subplots(figsize=(max(6, len(concordance_df) * 0.5), 4))
    ax.bar(concordance_df["clone_id"], concordance_df["concordance_rate"])
    ax.set_xlabel("Clone")
    ax.set_ylabel("All-caller concordance rate")
    ax.set_title("Caller concordance per clone")
    ax.set_xticklabels(concordance_df["clone_id"], rotation=45, ha="right")
    ax.set_ylim(0, 1)
    plt.tight_layout()
    fig.savefig(f"{out_dir}/caller_concordance.pdf")
    fig.savefig(f"{out_dir}/caller_concordance.png", dpi=150)
    plt.close(fig)


def plot_variant_counts_per_clone(
    sharing_df: pd.DataFrame,
    out_dir: str,
) -> None:
    """Stacked bar: private / partial / shared variants per clone."""
    if "sharing_class" not in sharing_df.columns:
        return

    clone_cols = [c for c in sharing_df.columns
                  if c not in ("sharing_class", "n_clones_with_variant")]

    counts = pd.DataFrame()
    for clone in clone_cols:
        sub = sharing_df.groupby("sharing_class")[clone].sum()
        counts[clone] = sub

    counts = counts.T
    fig, ax = plt.subplots(figsize=(max(6, len(clone_cols) * 0.6), 4))
    bottom = np.zeros(len(counts))
    colors = {"private": "#e74c3c", "partial": "#f39c12", "shared": "#2ecc71"}
    for cls in ["private", "partial", "shared"]:
        if cls in counts.columns:
            ax.bar(
                counts.index,
                counts[cls],
                bottom=bottom,
                label=cls,
                color=colors.get(cls, "gray"),
            )
            bottom += counts[cls].values

    ax.set_xlabel("Clone")
    ax.set_ylabel("Number of variants")
    ax.set_title("Variant sharing across clones")
    ax.legend()
    ax.set_xticklabels(counts.index, rotation=45, ha="right")
    plt.tight_layout()
    fig.savefig(f"{out_dir}/variant_sharing.pdf")
    fig.savefig(f"{out_dir}/variant_sharing.png", dpi=150)
    plt.close(fig)


def plot_clone_similarity_heatmap(
    matrix_df: pd.DataFrame,
    out_dir: str,
) -> None:
    """Jaccard similarity heatmap of variant profiles across clones."""
    try:
        from utils.vcf_utils import build_presence_absence_matrix
        clone_presence = (
            matrix_df.groupby(level="clone_id", axis=1).max()
        )
        n = clone_presence.shape[1]
        sim = np.zeros((n, n))
        clones = clone_presence.columns.tolist()
        for i, c1 in enumerate(clones):
            for j, c2 in enumerate(clones):
                a = clone_presence[c1].values.astype(bool)
                b = clone_presence[c2].values.astype(bool)
                inter = (a & b).sum()
                union = (a | b).sum()
                sim[i, j] = inter / union if union > 0 else 0.0

        sim_df = pd.DataFrame(sim, index=clones, columns=clones)
        fig, ax = plt.subplots(figsize=(max(5, n), max(4, n)))
        sns.heatmap(sim_df, annot=True, fmt=".2f", cmap="Blues",
                    vmin=0, vmax=1, ax=ax, square=True)
        ax.set_title("Clone pairwise variant Jaccard similarity")
        plt.tight_layout()
        fig.savefig(f"{out_dir}/clone_similarity_heatmap.pdf")
        fig.savefig(f"{out_dir}/clone_similarity_heatmap.png", dpi=150)
        plt.close(fig)
    except Exception as e:
        log.warning(f"Could not generate similarity heatmap: {e}")


def main() -> None:
    args = parse_args()

    Path(args.out_per_clone).mkdir(parents=True, exist_ok=True)
    Path(args.out_plots).mkdir(parents=True, exist_ok=True)

    # Parse VCF filenames to extract clone_id and caller
    vcf_map: Dict[Tuple[str, str], str] = {}
    for vcf_path in args.vcf_list:
        if not Path(vcf_path).exists():
            log.warning(f"VCF not found, skipping: {vcf_path}")
            continue
        try:
            clone_id, caller = parse_vcf_filename(vcf_path)
            vcf_map[(clone_id, caller)] = vcf_path
        except ValueError as e:
            log.warning(str(e))
            continue

    if not vcf_map:
        log.error("No valid VCF files found — check --vcf_list arguments")
        sys.exit(1)

    log.info(
        f"Loaded {len(vcf_map)} VCFs: "
        f"{len({k[0] for k in vcf_map})} clones × "
        f"{len({k[1] for k in vcf_map})} callers"
    )

    # Build presence/absence matrix
    log.info("Building variant presence/absence matrix …")
    matrix_df = build_presence_absence_matrix(vcf_map, pass_only=args.pass_only)
    log.info(f"  Matrix: {len(matrix_df)} unique variants × {matrix_df.shape[1]} (clone,caller)")

    if matrix_df.empty:
        log.warning("No variants detected across any clone/caller — empty outputs")
        pd.DataFrame().to_csv(args.out_matrix)
        pd.DataFrame().to_csv(args.out_concordance)
        pd.DataFrame().to_csv(args.out_private)
        return

    matrix_df.reset_index().to_csv(args.out_matrix, index=False)

    # Per-clone concordance
    log.info("Computing caller concordance …")
    concordance_df = caller_concordance_per_clone(matrix_df)
    concordance_df.to_csv(args.out_concordance, index=False)

    # Sharing classification
    log.info("Classifying variant sharing (private / partial / shared) …")
    sharing_df = classify_variant_sharing(matrix_df)
    sharing_df.reset_index().to_csv(args.out_private, index=False)

    n_private = (sharing_df["sharing_class"] == "private").sum()
    n_shared  = (sharing_df["sharing_class"] == "shared").sum()
    n_partial = (sharing_df["sharing_class"] == "partial").sum()
    log.info(
        f"  Private: {n_private}  |  Partial: {n_partial}  |  Shared: {n_shared}"
    )

    # Per-clone VCF tables (TSV for easy inspection)
    clones = matrix_df.columns.get_level_values("clone_id").unique()
    for clone in clones:
        # Filter to variants present in this clone (any caller)
        clone_presence = matrix_df.xs(clone, axis=1, level="clone_id")
        detected_keys = set(matrix_df.index[clone_presence.max(axis=1) > 0])
        clone_rows = []
        for (clone_id, caller), vcf_path in vcf_map.items():
            if clone_id != clone:
                continue
            caller_df = vcf_to_dataframe(vcf_path, pass_only=args.pass_only)
            if caller_df.empty:
                continue
            keys = zip(caller_df["chrom"], caller_df["pos"], caller_df["ref"], caller_df["alt"])
            caller_df = caller_df[[key in detected_keys for key in keys]].copy()
            if caller_df.empty:
                continue
            caller_df.insert(0, "caller", caller)
            clone_rows.append(caller_df)

        clone_df = (
            pd.concat(clone_rows, ignore_index=True)
            if clone_rows
            else pd.DataFrame(columns=["caller", "chrom", "pos", "ref", "alt", "qual", "filter"])
        )
        clone_df.to_csv(f"{args.out_per_clone}/{clone}_variants.tsv", sep="\t", index=False)

    # Plots
    log.info("Generating comparison plots …")
    plot_concordance_heatmap(concordance_df, args.out_plots)
    plot_variant_counts_per_clone(sharing_df.reset_index(), args.out_plots)
    plot_clone_similarity_heatmap(matrix_df, args.out_plots)

    log.info("variant comparison complete.")


if __name__ == "__main__":
    main()
