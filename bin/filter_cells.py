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
    p.add_argument("--min_cells_for_calling", type=int, default=2,
                   help="clones with fewer cells are not sent to mutation calling")
    p.add_argument("--min_cells_reliable",    type=int, default=20,
                   help="clones with fewer cells are still called but flagged LOW confidence")
    p.add_argument("--out_manifest",       required=True)
    p.add_argument("--out_qc_summary",     required=True)
    p.add_argument("--out_reliability",    required=True)
    p.add_argument("--out_reliability_md", required=True)
    return p.parse_args()


def write_reliability_md(rel: pd.DataFrame, path: str,
                         min_call: int, min_reliable: int) -> None:
    """Human-readable per-clone reliability report (the 'why' for the user)."""
    lines = [
        "# Clone reliability for mutation calling", "",
        f"A clone is **called** if it has at least **{min_call}** cells, and treated "
        f"as **reliable** if it has at least **{min_reliable}** cells. Clones below "
        f"the reliable threshold are still called — their VCFs are produced — but "
        f"their variant lists are exploratory, not confident.", "",
        "| Clone | Cells | ~Depth vs largest | Called | Reliable | Note |",
        "|-------|------:|------------------:|:------:|:--------:|------|",
    ]
    for _, r in rel.sort_values("n_cells", ascending=False).iterrows():
        reliable_cell = ("**yes**" if r["reliable"]
                         else ("low" if r["called"] else "—"))
        lines.append(
            f"| {r['clone_id']} | {int(r['n_cells'])} | {r['rel_coverage']:.0%} | "
            f"{'yes' if r['called'] else 'no'} | {reliable_cell} | {r['reason']} |"
        )
    lines += [
        "", "## Why small clones give unreliable calls", "",
        "- **Pseudobulk depth scales with cell count.** Coverage is roughly "
        "proportional to the number of merged cells, so a clone with few cells has "
        "too few reads at most positions to call somatic SNVs confidently.",
        "- **Tumor-only calling.** With no matched normal, germline variants and "
        "sequencing artefacts are not subtracted; they inflate the call set, and the "
        "effect is worst at low depth.",
        "- **scDNA-seq sparsity.** Per-cell coverage is sparse and uneven, so merging "
        "a handful of cells leaves coverage gaps across the genome.",
    ]
    with open(path, "w") as fh:
        fh.write("\n".join(lines) + "\n")


def main() -> None:
    args = parse_args()

    assignments = pd.read_csv(args.assignments)
    manifest    = pd.read_csv(args.bam_manifest)

    # ── Reconcile the two cell sets BEFORE the join, so dropped cells are loud
    #    instead of silently disappearing in an inner join. ──────────────────
    bam_cells    = set(manifest["cell_id"])
    assign_cells = set(assignments["cell_id"])
    missing_bam    = assign_cells - bam_cells   # have a clone, but no BAM to merge
    missing_assign = bam_cells - assign_cells   # have a BAM, but no clone label
    if missing_bam:
        ex = ", ".join(sorted(missing_bam)[:5])
        log.warning(
            f"{len(missing_bam)} cell(s) have a clone assignment but NO BAM in the "
            f"manifest -> excluded from pseudobulk (e.g. {ex})"
        )
    if missing_assign:
        ex = ", ".join(sorted(missing_assign)[:5])
        log.warning(
            f"{len(missing_assign)} cell(s) have a BAM but NO clone assignment -> "
            f"excluded (e.g. {ex})"
        )

    # Merge assignments onto manifest (inner: a cell needs both a clone and a BAM)
    merged = manifest.merge(
        assignments[["cell_id", "clone_id", "confidence", "flagged"]],
        on="cell_id",
        how="inner",
    )
    n_start = len(merged)
    log.info(
        f"Starting with {n_start} cells across {merged['clone_id'].nunique()} "
        f"clone(s) (cells present in BOTH manifest and assignments)"
    )

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

    # ── Calling gate + reliability ───────────────────────────────────────────
    # Two different questions: CAN a clone be called (does it have a usable
    # pseudobulk) and is that call RELIABLE (enough depth to trust). We call
    # broadly and flag reliability honestly rather than silently dropping clones.
    clone_counts = merged.groupby("clone_id")["cell_id"].count().sort_index()
    max_cells = int(clone_counts.max()) if len(clone_counts) else 0

    rel_rows, called_clones = [], []
    for clone_id, n in clone_counts.items():
        n = int(n)
        called   = n >= args.min_cells_for_calling
        reliable = n >= args.min_cells_reliable
        rel_cov  = (n / max_cells) if max_cells else 0.0
        if not called:
            reason = (f"only {n} cells (< {args.min_cells_for_calling}) — too few for "
                      f"a usable pseudobulk; not called")
            log.warning(f"  SKIP      {clone_id}: {reason}")
        elif reliable:
            reason = f"{n} cells — sufficient pseudobulk depth"
            called_clones.append(clone_id)
            log.info(f"  RELIABLE  {clone_id}: {n} cells")
        else:
            reason = (f"only {n} cells (< {args.min_cells_reliable}); ~{rel_cov:.0%} of "
                      f"the largest clone's depth — LOW-confidence calls (sparse "
                      f"coverage; tumor-only artefacts/germline not subtracted)")
            called_clones.append(clone_id)
            log.warning(f"  LOW-CONF  {clone_id}: {reason}")
        rel_rows.append({"clone_id": clone_id, "n_cells": n,
                         "rel_coverage": round(rel_cov, 4),
                         "called": called, "reliable": reliable, "reason": reason})

    reliability = pd.DataFrame(rel_rows)
    reliability.to_csv(args.out_reliability, index=False)
    write_reliability_md(reliability, args.out_reliability_md,
                         args.min_cells_for_calling, args.min_cells_reliable)

    if not called_clones:
        log.error(f"No clone has >= {args.min_cells_for_calling} cells — nothing to call.")

    # filtered_manifest feeds MERGE_BAMS: every CALLED clone goes forward.
    out = merged[merged["clone_id"].isin(called_clones)]
    n_reliable = int(reliability["reliable"].sum())
    log.info(f"Pseudobulk manifest: {len(out)} cells across {len(called_clones)} "
             f"called clone(s) ({n_reliable} reliable) of "
             f"{merged['clone_id'].nunique()} total")

    # Summary reports ALL clones with called/reliable flags for transparency.
    summary = (merged.groupby("clone_id")
               .agg(n_cells=("cell_id", "count")).reset_index())
    summary["pct_passed"] = (summary["n_cells"] / n_start * 100).round(2)
    summary = summary.merge(reliability[["clone_id", "called", "reliable"]],
                            on="clone_id", how="left")

    out.to_csv(args.out_manifest, index=False)
    summary.to_csv(args.out_qc_summary, index=False)


if __name__ == "__main__":
    main()
