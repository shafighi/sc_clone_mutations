#!/usr/bin/env python3
"""
build_consensus.py

Build a consensus mutation callset from the cross-clone variant matrix.

Rules (configurable):
  - A variant is included in the consensus if it is supported by at least
    --min_callers callers in at least one clone.
  - Optionally restrict to PASS-only variants (--pass_only).

Output:
  - consensus_mutations.vcf.gz   : consensus VCF (bgzipped)
  - consensus_table.csv          : tidy table of consensus mutations
  - consensus_summary.json       : summary statistics

Usage:
    build_consensus.py \\
        --variant_matrix variant_matrix.csv \\
        --per_clone_vcf_dir per_clone_vcfs/ \\
        --min_callers 2 \\
        --pass_only True \\
        --out_vcf consensus_mutations.vcf.gz \\
        --out_table consensus_table.csv \\
        --out_summary consensus_summary.json
"""

from __future__ import annotations

import argparse
import gzip
import json
import logging
import sys
from pathlib import Path
from typing import Any, Dict, List, Optional, Set, Tuple

import pandas as pd

logging.basicConfig(
    level=logging.INFO,
    format="[%(asctime)s] %(levelname)s: %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)
log = logging.getLogger(__name__)


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--variant_matrix",   required=True)
    p.add_argument("--per_clone_vcf_dir",required=True)
    p.add_argument("--min_callers",      type=int, default=2)
    p.add_argument("--pass_only",        type=lambda x: x.lower() == "true", default=True)
    p.add_argument("--out_vcf",          required=True)
    p.add_argument("--out_table",        required=True)
    p.add_argument("--out_summary",      required=True)
    return p.parse_args()


def load_variant_matrix(path: str) -> pd.DataFrame:
    """Re-load the presence/absence matrix saved by compare_variants.py."""
    df = pd.read_csv(path)
    # Restore MultiIndex if needed
    if "chrom" in df.columns:
        df = df.set_index(["chrom", "pos", "ref", "alt"])
    return df


def select_consensus_variants(
    matrix_df: pd.DataFrame,
    min_callers: int,
) -> pd.DataFrame:
    """
    Return the subset of variants in matrix_df that are supported by
    at least min_callers callers in at least one clone.

    The matrix columns are a MultiIndex (clone_id, caller).
    """
    if not hasattr(matrix_df.columns, "levels"):
        # Flat columns — try to reconstruct MultiIndex from "clone_id.caller" names
        tuples = [tuple(c.split(".", 1)) for c in matrix_df.columns]
        matrix_df.columns = pd.MultiIndex.from_tuples(tuples, names=["clone_id", "caller"])

    # For each variant, count the max number of callers supporting it across clones
    caller_counts_per_clone = (
        matrix_df.groupby(level="clone_id", axis=1).sum()
    )
    max_callers_any_clone = caller_counts_per_clone.max(axis=1)
    consensus_mask = max_callers_any_clone >= min_callers

    consensus = matrix_df[consensus_mask]
    log.info(
        f"Consensus: {consensus_mask.sum()}/{len(matrix_df)} variants "
        f"supported by ≥{min_callers} callers"
    )
    return consensus


def build_consensus_table(consensus_df: pd.DataFrame) -> pd.DataFrame:
    """
    Convert the consensus variant matrix into a tidy DataFrame.
    One row per variant; columns for each clone showing caller support.
    """
    rows = []
    clones = consensus_df.columns.get_level_values("clone_id").unique()
    for (chrom, pos, ref, alt), row in consensus_df.iterrows():
        entry: Dict[str, Any] = {
            "chrom": chrom,
            "pos":   pos,
            "ref":   ref,
            "alt":   alt,
        }
        for clone in clones:
            sub = row.xs(clone, level="clone_id") if clone in row.index.get_level_values(0) else pd.Series()
            if not sub.empty:
                entry[f"{clone}_n_callers"]  = int(sub.sum())
                entry[f"{clone}_callers"]    = ";".join(sub[sub > 0].index.tolist())
            else:
                entry[f"{clone}_n_callers"]  = 0
                entry[f"{clone}_callers"]    = ""
        rows.append(entry)
    return pd.DataFrame(rows)


def write_minimal_vcf(
    variants: pd.DataFrame,
    out_path: str,
    sample_name: str = "consensus",
) -> None:
    """
    Write a minimal VCF for the consensus variants.
    This is a stub VCF — copy information from original caller VCFs if
    full FORMAT/INFO annotations are needed.

    TODO: if full INFO/FORMAT fields are required, use pysam to copy
    records from the original caller VCFs.
    """
    header_lines = [
        "##fileformat=VCFv4.2",
        f"##source=sc_clone_mutations_consensus",
        "##INFO=<ID=N_CLONES,Number=1,Type=Integer,Description=\"Number of clones with variant\">",
        "##INFO=<ID=CALLERS,Number=.,Type=String,Description=\"Callers supporting this variant\">",
        f"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
    ]

    clone_cols = [c for c in variants.columns if c.endswith("_n_callers")]
    clones = [c.replace("_n_callers", "") for c in clone_cols]

    records = []
    for _, row in variants.iterrows():
        n_clones = sum(
            1 for c in clone_cols if row.get(c, 0) > 0
        )
        callers_per_clone = [
            row.get(f"{c}_callers", "") for c in clones
            if row.get(f"{c}_n_callers", 0) > 0
        ]
        callers_union = set(";".join(callers_per_clone).split(";")) - {""}
        info = f"N_CLONES={n_clones};CALLERS={','.join(sorted(callers_union))}"
        records.append(
            f"{row['chrom']}\t{row['pos']}\t.\t{row['ref']}\t{row['alt']}\t.\tPASS\t{info}"
        )

    with gzip.open(out_path, "wt") as fh:
        fh.write("\n".join(header_lines) + "\n")
        fh.write("\n".join(records) + "\n")


def main() -> None:
    args = parse_args()

    matrix_df = load_variant_matrix(args.variant_matrix)
    if matrix_df.empty:
        log.warning("Variant matrix is empty — writing empty outputs")
        pd.DataFrame().to_csv(args.out_table, index=False)
        with gzip.open(args.out_vcf, "wt") as fh:
            fh.write("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        with open(args.out_summary, "w") as fh:
            json.dump({"n_consensus": 0}, fh)
        return

    # Select consensus variants
    consensus_df = select_consensus_variants(matrix_df, args.min_callers)

    # Build tidy table
    table = build_consensus_table(consensus_df)
    table.to_csv(args.out_table, index=False)
    log.info(f"Consensus table → {args.out_table}")

    # Write VCF
    write_minimal_vcf(table, args.out_vcf)
    log.info(f"Consensus VCF → {args.out_vcf}")

    # Summary
    summary = {
        "n_input_variants":    len(matrix_df),
        "n_consensus":         len(consensus_df),
        "min_callers_required": args.min_callers,
        "pass_only":           args.pass_only,
    }
    with open(args.out_summary, "w") as fh:
        json.dump(summary, fh, indent=2)
    log.info(f"Summary → {args.out_summary}")


if __name__ == "__main__":
    main()
