#!/usr/bin/env python3
"""
filter_cells.py

Filter cells before pseudobulking based on:
  - Clone assignment confidence score
  - Minimum mapped reads (from flagstat, if available)
  - Maximum duplication rate (from markdup metrics, if available)

Cells failing QC are excluded from the output manifest. A QC summary is written.

Usage:
    filter_cells.py \\
        --assignments cell_clone_assignments.csv \\
        --bam_manifest bam_manifest.csv \\
        --min_mapped_reads 100000 \\
        --max_dup_rate 0.95 \\
        --min_confidence 0.0 \\
        --out_manifest filtered_manifest.csv \\
        --out_qc_summary cell_qc_summary.csv
"""

from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path

import pandas as pd

logging.basicConfig(
    level=logging.INFO,
    format="[%(asctime)s] %(levelname)s: %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)
log = logging.getLogger(__name__)


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--assignments",       required=True)
    p.add_argument("--bam_manifest",      required=True)
    p.add_argument("--min_mapped_reads",  type=int,   default=100_000)
    p.add_argument("--max_dup_rate",      type=float, default=0.95)
    p.add_argument("--min_confidence",    type=float, default=0.0)
    p.add_argument("--out_manifest",      required=True)
    p.add_argument("--out_qc_summary",    required=True)
    return p.parse_args()


def main() -> None:
    args = parse_args()

    assignments = pd.read_csv(args.assignments)
    manifest    = pd.read_csv(args.bam_manifest)

    # Merge assignments onto manifest
    merged = manifest.merge(
        assignments[["cell_id", "clone_id", "confidence", "flagged"]],
        on="cell_id",
        how="inner",
    )
    n_start = len(merged)
    log.info(f"Starting with {n_start} cells (after inner join with assignments)")

    qc_rows = []

    # Filter by confidence
    mask_conf = merged["confidence"] >= args.min_confidence
    dropped_conf = (~mask_conf).sum()
    if dropped_conf:
        log.info(f"Dropping {dropped_conf} cells below confidence {args.min_confidence}")
    merged = merged[mask_conf]

    # NOTE: Filtering by mapped reads and duplication rate would require
    # running samtools flagstat / picard MarkDuplicates on each cell BAM first.
    # In this pipeline, those are run on the PSEUDOBULK, not individual cells.
    # If per-cell QC metrics are available in cell_metadata, they can be joined here.
    # TODO: join per-cell QC metrics from cell_metadata if available, then filter.
    log.info(
        "Note: per-cell mapped-read and duplication-rate filters require "
        "pre-computed per-cell QC metrics (not computed here by default)."
    )

    # Build summary
    summary = (
        merged.groupby("clone_id")
        .agg(n_cells=("cell_id", "count"))
        .reset_index()
    )
    summary["pct_passed"] = (summary["n_cells"] / n_start * 100).round(2)

    n_final = len(merged)
    log.info(
        f"After filtering: {n_final}/{n_start} cells retained "
        f"across {merged['clone_id'].nunique()} clones"
    )

    # Warn about clones with very few cells
    low_clones = summary[summary["n_cells"] < 5]
    if not low_clones.empty:
        log.warning(
            f"{len(low_clones)} clone(s) have <5 cells after QC filtering: "
            + ", ".join(low_clones["clone_id"].tolist())
        )

    merged.to_csv(args.out_manifest, index=False)
    summary.to_csv(args.out_qc_summary, index=False)


if __name__ == "__main__":
    main()
