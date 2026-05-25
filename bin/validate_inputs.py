#!/usr/bin/env python3
"""
validate_inputs.py

Validates all input manifests before the pipeline runs, failing fast with
helpful error messages rather than cryptic process failures mid-run.

Checks:
  - Required columns present in each manifest
  - BAM / BAI files exist and are non-empty
  - Cell IDs consistent between BAM manifest and cell metadata
  - MEDICC2 tree is parseable Newick
  - scUnique events file has required columns
  - Normal manifest (if provided) is valid

Usage:
    validate_inputs.py \\
        --bam_manifest bam_manifest.csv \\
        --cell_metadata cell_metadata.csv \\
        --medicc2_tree tree.new \\
        --scunique_events events.tsv \\
        [--normal_manifest normal.csv] \\
        --out_bam_manifest bam_manifest_validated.csv \\
        --out_cell_metadata cell_metadata_validated.csv \\
        --out_report validation_report.json \\
        --min_mapped_reads 100000
"""

from __future__ import annotations

import argparse
import json
import logging
import sys
from pathlib import Path
from typing import Any, Dict, List, Optional

import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
from utils.tree_utils import load_newick

logging.basicConfig(
    level=logging.INFO,
    format="[%(asctime)s] %(levelname)s: %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)
log = logging.getLogger(__name__)


BAM_MANIFEST_REQUIRED_COLS  = {"cell_id", "bam_path", "sample_id", "patient_id"}
CELL_METADATA_REQUIRED_COLS  = {"cell_id"}
NORMAL_MANIFEST_REQUIRED_COLS = {"sample_id", "patient_id", "bam_path"}
SCUNIQUE_REQUIRED_COLS        = {"cell_id"}


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--bam_manifest",      required=True)
    p.add_argument("--cell_metadata",     required=True)
    p.add_argument("--medicc2_tree",      required=True)
    p.add_argument("--scunique_events",   required=True)
    p.add_argument("--normal_manifest",   default=None)
    p.add_argument("--out_bam_manifest",  required=True)
    p.add_argument("--out_cell_metadata", required=True)
    p.add_argument("--out_report",        required=True)
    p.add_argument("--min_mapped_reads",  type=int, default=100_000)
    p.add_argument("--skip_bam_check",    action="store_true",
                   help="Skip BAM file existence checks (for smoke tests with dummy paths)")
    return p.parse_args()


# ─── Validation helpers ───────────────────────────────────────────────────────

class ValidationError(Exception):
    pass


def check_required_columns(df: pd.DataFrame, required: set, name: str) -> None:
    missing = required - set(df.columns)
    if missing:
        raise ValidationError(
            f"{name} is missing required columns: {sorted(missing)}\n"
            f"  Found columns: {sorted(df.columns.tolist())}"
        )


def validate_bam_manifest(
    path: str,
    min_mapped_reads: int,
    skip_bam_check: bool = False,
) -> pd.DataFrame:
    log.info(f"Validating BAM manifest: {path}")
    df = pd.read_csv(path)
    check_required_columns(df, BAM_MANIFEST_REQUIRED_COLS, "BAM manifest")

    errors: List[str] = []
    missing_bams: List[str] = []
    missing_bais: List[str] = []

    if not skip_bam_check:
        for _, row in df.iterrows():
            bam = Path(str(row["bam_path"]))
            if not bam.exists():
                missing_bams.append(str(bam))
                continue
            if bam.stat().st_size == 0:
                errors.append(f"BAM is empty: {bam}")

            # Infer BAI path
            bai = Path(str(row.get("bai_path", f"{bam}.bai")))
            if not bai.exists():
                bai_alt = bam.with_suffix(".bai")
                if not bai_alt.exists():
                    missing_bais.append(str(bam))
                else:
                    df.loc[df["bam_path"] == row["bam_path"], "bai_path"] = str(bai_alt)
            else:
                df.loc[df["bam_path"] == row["bam_path"], "bai_path"] = str(bai)

        if missing_bams:
            errors.append(
                f"{len(missing_bams)} BAM file(s) not found:\n  "
                + "\n  ".join(missing_bams[:10])
                + ("  ..." if len(missing_bams) > 10 else "")
            )
        if missing_bais:
            log.warning(
                f"{len(missing_bais)} BAM index (.bai) file(s) missing — "
                "Nextflow will need to index them or they must be present."
            )

        if errors:
            raise ValidationError("\n".join(errors))
    else:
        log.warning("BAM file existence check skipped (--skip_bam_check)")

    # Deduplicate cell IDs
    dupes = df[df.duplicated("cell_id", keep=False)]["cell_id"].unique()
    if len(dupes) > 0:
        raise ValidationError(
            f"Duplicate cell_id(s) in BAM manifest: {list(dupes[:10])}"
        )

    log.info(f"  {len(df)} cells validated ({df['sample_id'].nunique()} samples)")
    return df


def validate_cell_metadata(path: str) -> pd.DataFrame:
    log.info(f"Validating cell metadata: {path}")
    df = pd.read_csv(path)
    check_required_columns(df, CELL_METADATA_REQUIRED_COLS, "cell metadata")
    dupes = df[df.duplicated("cell_id", keep=False)]["cell_id"].unique()
    if len(dupes) > 0:
        raise ValidationError(f"Duplicate cell_id(s) in cell metadata: {list(dupes[:10])}")
    log.info(f"  {len(df)} cells in metadata")
    return df


def validate_medicc2_tree(path: str) -> Dict[str, Any]:
    log.info(f"Validating MEDICC2 tree: {path}")
    try:
        tree = load_newick(path)
        n_leaves = sum(1 for _ in tree.leaf_node_iter())
        n_nodes  = sum(1 for _ in tree.preorder_node_iter())
    except Exception as e:
        raise ValidationError(f"Cannot parse MEDICC2 tree: {e}")
    if n_leaves < 2:
        raise ValidationError(f"Tree has only {n_leaves} leaf — need ≥2 cells")
    log.info(f"  Tree OK: {n_leaves} leaves, {n_nodes} total nodes")
    return {"n_leaves": n_leaves, "n_nodes": n_nodes}


def validate_scunique_events(path: str) -> Dict[str, Any]:
    log.info(f"Validating scUnique events: {path}")
    sep = "\t" if path.endswith(".tsv") else ","
    df  = pd.read_csv(path, sep=sep)
    check_required_columns(df, SCUNIQUE_REQUIRED_COLS, "scUnique events")
    n_cells = df["cell_id"].nunique()
    n_events = len(df)
    log.info(f"  {n_events} events across {n_cells} cells")
    return {"n_cells": n_cells, "n_events": n_events}


def validate_normal_manifest(path: str) -> pd.DataFrame:
    if not path or path in ("NO_FILE", "null", ""):
        return pd.DataFrame()
    log.info(f"Validating normal manifest: {path}")
    df = pd.read_csv(path)
    check_required_columns(df, NORMAL_MANIFEST_REQUIRED_COLS, "normal manifest")
    errors: List[str] = []
    for _, row in df.iterrows():
        bam = Path(str(row["bam_path"]))
        if not bam.exists():
            errors.append(str(bam))
    if errors:
        raise ValidationError(
            f"{len(errors)} normal BAM(s) not found:\n  " + "\n  ".join(errors[:5])
        )
    log.info(f"  {len(df)} normal sample(s) validated")
    return df


def cross_validate(
    bam_df: pd.DataFrame,
    cell_meta_df: pd.DataFrame,
    tree_info: Dict[str, Any],
) -> Dict[str, Any]:
    """Check consistency across inputs."""
    bam_cells  = set(bam_df["cell_id"].astype(str))
    meta_cells = set(cell_meta_df["cell_id"].astype(str))

    in_bam_not_meta   = bam_cells - meta_cells
    in_meta_not_bam   = meta_cells - bam_cells
    in_both           = bam_cells & meta_cells

    if in_bam_not_meta:
        log.warning(
            f"{len(in_bam_not_meta)} cells in BAM manifest but not in cell metadata"
        )
    if in_meta_not_bam:
        log.warning(
            f"{len(in_meta_not_bam)} cells in cell metadata but not in BAM manifest"
        )

    overlap_pct = 100 * len(in_both) / max(len(bam_cells), 1)
    if overlap_pct < 50:
        raise ValidationError(
            f"Only {overlap_pct:.1f}% of BAM cells match cell metadata — "
            "check that cell_id naming is consistent."
        )

    return {
        "n_bam_cells":           len(bam_cells),
        "n_meta_cells":          len(meta_cells),
        "n_overlap":             len(in_both),
        "pct_overlap":           round(overlap_pct, 2),
        "n_in_bam_not_meta":     len(in_bam_not_meta),
        "n_in_meta_not_bam":     len(in_meta_not_bam),
        "n_tree_leaves":         tree_info["n_leaves"],
    }


# ─── Main ────────────────────────────────────────────────────────────────────

def main() -> None:
    args = parse_args()
    errors: List[str] = []
    report: Dict[str, Any] = {"status": "ok", "warnings": [], "errors": []}

    # Validate each input
    try:
        bam_df       = validate_bam_manifest(
            args.bam_manifest, args.min_mapped_reads,
            skip_bam_check=args.skip_bam_check,
        )
    except ValidationError as e:
        errors.append(f"BAM manifest: {e}")
        bam_df = pd.DataFrame()

    try:
        cell_meta_df = validate_cell_metadata(args.cell_metadata)
    except ValidationError as e:
        errors.append(f"Cell metadata: {e}")
        cell_meta_df = pd.DataFrame()

    try:
        tree_info    = validate_medicc2_tree(args.medicc2_tree)
    except ValidationError as e:
        errors.append(f"MEDICC2 tree: {e}")
        tree_info = {}

    try:
        event_info   = validate_scunique_events(args.scunique_events)
        report["scunique"] = event_info
    except ValidationError as e:
        errors.append(f"scUnique events: {e}")

    try:
        validate_normal_manifest(args.normal_manifest)
    except ValidationError as e:
        errors.append(f"Normal manifest: {e}")

    # Cross-validate
    if bam_df is not None and not bam_df.empty and not cell_meta_df.empty and tree_info:
        try:
            cross_info = cross_validate(bam_df, cell_meta_df, tree_info)
            report["cross_validation"] = cross_info
        except ValidationError as e:
            errors.append(f"Cross-validation: {e}")

    # Fail if any hard errors
    if errors:
        report["status"] = "failed"
        report["errors"] = errors
        with open(args.out_report, "w") as fh:
            json.dump(report, fh, indent=2)
        for err in errors:
            log.error(err)
        sys.exit(1)

    # Write validated outputs (write even if warnings)
    bam_df.to_csv(args.out_bam_manifest, index=False)
    cell_meta_df.to_csv(args.out_cell_metadata, index=False)

    report["status"]     = "ok"
    report["tree"]       = tree_info
    report["n_bam_cells"] = len(bam_df)

    with open(args.out_report, "w") as fh:
        json.dump(report, fh, indent=2)

    log.info("All input validation checks passed.")


if __name__ == "__main__":
    main()
