#!/usr/bin/env python3
"""
generate_report.py

Generate a custom HTML/Markdown summary report for sc_clone_mutations runs.

Sections:
  1. Run overview (inputs, parameters)
  2. Clone assignment summary
  3. Pseudobulk QC (loaded from mosdepth/flagstat outputs if present)
  4. Variant calling counts per caller and per clone
  5. Consensus mutation summary
  6. Cross-clone variant sharing
  7. Warnings and caveats

Usage:
    generate_report.py \\
        --validation_report validation_report.json \\
        --clone_summary clone_summary.csv \\
        --consensus_table consensus_table.csv \\
        --cross_clone_matrix variant_matrix.csv \\
        --pipeline_version 1.0.0 \\
        --out_html sc_clone_mutations_report.html \\
        --out_md   sc_clone_mutations_report.md
"""

from __future__ import annotations

import argparse
import datetime
import json
import logging
import sys
from pathlib import Path
from typing import Any, Dict, List, Optional

import pandas as pd

logging.basicConfig(
    level=logging.INFO,
    format="[%(asctime)s] %(levelname)s: %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)
log = logging.getLogger(__name__)


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--validation_report",  required=True)
    p.add_argument("--clone_summary",      required=True)
    p.add_argument("--consensus_table",    required=True)
    p.add_argument("--cross_clone_matrix", required=True)
    p.add_argument("--pipeline_version",   default="unknown")
    p.add_argument("--out_html",           required=True)
    p.add_argument("--out_md",             required=True)
    return p.parse_args()


def load_json(path: str) -> Dict[str, Any]:
    try:
        with open(path) as fh:
            return json.load(fh)
    except Exception as e:
        log.warning(f"Could not load {path}: {e}")
        return {}


def load_csv_safe(path: str) -> pd.DataFrame:
    try:
        return pd.read_csv(path)
    except Exception as e:
        log.warning(f"Could not load {path}: {e}")
        return pd.DataFrame()


def df_to_html(df: pd.DataFrame, max_rows: int = 50) -> str:
    if df.empty:
        return "<p><em>No data available.</em></p>"
    return df.head(max_rows).to_html(
        index=False, border=0, classes="table table-striped table-sm",
        justify="left",
    )


def df_to_md(df: pd.DataFrame, max_rows: int = 50) -> str:
    if df.empty:
        return "_No data available._\n"
    return df.head(max_rows).to_markdown(index=False) + "\n"


def build_markdown(
    val_report: Dict,
    clone_summary: pd.DataFrame,
    consensus_table: pd.DataFrame,
    variant_matrix: pd.DataFrame,
    version: str,
) -> str:
    ts = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    lines: List[str] = []

    lines.append(f"# sc_clone_mutations Report")
    lines.append(f"\n**Pipeline version:** {version}  ")
    lines.append(f"**Generated:** {ts}\n")

    # ── 1. Run overview ──────────────────────────────────────────────────────
    lines.append("## 1. Input Validation")
    status = val_report.get("status", "unknown")
    lines.append(f"- Status: **{status}**")
    if "n_bam_cells" in val_report:
        lines.append(f"- Cells in BAM manifest: {val_report['n_bam_cells']}")
    if "tree" in val_report:
        lines.append(f"- Tree leaves: {val_report['tree'].get('n_leaves', 'N/A')}")
    if "scunique" in val_report:
        sc = val_report["scunique"]
        lines.append(
            f"- scUnique events: {sc.get('n_events', 'N/A')} events "
            f"in {sc.get('n_cells', 'N/A')} cells"
        )
    if val_report.get("errors"):
        lines.append("\n### Validation Errors")
        for err in val_report["errors"]:
            lines.append(f"- ⚠ {err}")
    lines.append("")

    # ── 2. Clone summary ─────────────────────────────────────────────────────
    lines.append("## 2. Clone Assignment Summary")
    if not clone_summary.empty:
        lines.append(f"- Number of clones: **{len(clone_summary)}**")
        total_cells = clone_summary["n_cells"].sum() if "n_cells" in clone_summary.columns else "N/A"
        lines.append(f"- Total cells assigned: {total_cells}")
        lines.append("")
        lines.append(df_to_md(clone_summary))
    else:
        lines.append("_Clone summary not available._\n")

    # ── 3. Consensus mutations ────────────────────────────────────────────────
    lines.append("## 3. Consensus Mutation Summary")
    if not consensus_table.empty:
        lines.append(f"- Consensus variants: **{len(consensus_table)}**")
        lines.append("")
        lines.append(df_to_md(consensus_table.head(20)))
        if len(consensus_table) > 20:
            lines.append(f"_... and {len(consensus_table) - 20} more variants_\n")
    else:
        lines.append("_No consensus mutations available._\n")

    # ── 4. Cross-clone sharing ────────────────────────────────────────────────
    lines.append("## 4. Cross-Clone Variant Sharing")
    if not variant_matrix.empty:
        lines.append(f"- Total unique variant positions: {len(variant_matrix)}")
        if "sharing_class" in variant_matrix.columns:
            sharing_counts = variant_matrix["sharing_class"].value_counts()
            for cls, count in sharing_counts.items():
                lines.append(f"  - {cls}: {count}")
    else:
        lines.append("_Variant matrix not available._")
    lines.append("")

    # ── 5. Warnings ───────────────────────────────────────────────────────────
    warnings: List[str] = []
    if not clone_summary.empty and "n_cells" in clone_summary.columns:
        low = clone_summary[clone_summary["n_cells"] < 5]
        for _, row in low.iterrows():
            warnings.append(f"Clone {row['clone_id']} has only {row['n_cells']} cells — low statistical power")

    if warnings:
        lines.append("## 5. Warnings")
        for w in warnings:
            lines.append(f"- ⚠ {w}")
        lines.append("")

    return "\n".join(lines)


HTML_TEMPLATE = """<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>sc_clone_mutations Report</title>
<link rel="stylesheet"
  href="https://cdn.jsdelivr.net/npm/bootstrap@5.3.0/dist/css/bootstrap.min.css">
<style>
  body {{ font-family: 'Segoe UI', Arial, sans-serif; padding: 2rem; }}
  h1   {{ color: #2c3e50; }}
  h2   {{ color: #34495e; border-bottom: 2px solid #ecf0f1; padding-bottom: 0.3rem; }}
  .badge-ok      {{ background-color: #27ae60; }}
  .badge-failed  {{ background-color: #e74c3c; }}
  .table-sm td, .table-sm th {{ font-size: 0.85rem; }}
  .warning-box {{ background: #fff3cd; border-left: 4px solid #ffc107; padding: 0.8rem 1rem; }}
</style>
</head>
<body>
<div class="container-fluid">

<h1>sc_clone_mutations Report</h1>
<p class="text-muted">Pipeline version: <strong>{version}</strong> &nbsp;|&nbsp;
Generated: <strong>{timestamp}</strong></p>
<hr>

<h2>1. Input Validation</h2>
<p>Status: <span class="badge {badge_class}">{status}</span></p>
{validation_details}
{validation_errors}

<h2>2. Clone Assignment</h2>
{clone_table}

<h2>3. Consensus Mutations</h2>
{consensus_table}

<h2>4. Cross-Clone Variant Sharing</h2>
{sharing_section}

{warnings_section}

</div>
</body>
</html>"""


def build_html(
    val_report: Dict,
    clone_summary: pd.DataFrame,
    consensus_table: pd.DataFrame,
    variant_matrix: pd.DataFrame,
    version: str,
) -> str:
    ts = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    status = val_report.get("status", "unknown")
    badge_class = "badge-ok" if status == "ok" else "badge-failed"

    val_details_lines = []
    if "n_bam_cells" in val_report:
        val_details_lines.append(f"<li>Cells in BAM manifest: {val_report['n_bam_cells']}</li>")
    if "tree" in val_report:
        val_details_lines.append(f"<li>Tree leaves: {val_report['tree'].get('n_leaves','N/A')}</li>")
    validation_details = f"<ul>{''.join(val_details_lines)}</ul>" if val_details_lines else ""

    errors = val_report.get("errors", [])
    if errors:
        err_items = "".join(f"<li>{e}</li>" for e in errors)
        validation_errors = f"<div class='warning-box'><strong>Errors:</strong><ul>{err_items}</ul></div>"
    else:
        validation_errors = ""

    clone_table_html = df_to_html(clone_summary) if not clone_summary.empty else "<p>Not available.</p>"

    consensus_html = (
        df_to_html(consensus_table.head(30), max_rows=30) if not consensus_table.empty
        else "<p>No consensus mutations.</p>"
    )

    if not variant_matrix.empty and "sharing_class" in variant_matrix.columns:
        sc = variant_matrix["sharing_class"].value_counts().reset_index()
        sc.columns = ["Category", "Count"]
        sharing_section = df_to_html(sc)
    else:
        sharing_section = "<p>Not available.</p>"

    warnings_items = []
    if not clone_summary.empty and "n_cells" in clone_summary.columns:
        low = clone_summary[clone_summary["n_cells"] < 5]
        for _, row in low.iterrows():
            warnings_items.append(
                f"<li>Clone <strong>{row['clone_id']}</strong> has only "
                f"{row['n_cells']} cells — low statistical power.</li>"
            )
    if warnings_items:
        warnings_section = (
            "<h2>5. Warnings</h2>"
            f"<div class='warning-box'><ul>{''.join(warnings_items)}</ul></div>"
        )
    else:
        warnings_section = ""

    return HTML_TEMPLATE.format(
        version=version,
        timestamp=ts,
        status=status,
        badge_class=badge_class,
        validation_details=validation_details,
        validation_errors=validation_errors,
        clone_table=clone_table_html,
        consensus_table=consensus_html,
        sharing_section=sharing_section,
        warnings_section=warnings_section,
    )


def main() -> None:
    args = parse_args()

    val_report     = load_json(args.validation_report)
    clone_summary  = load_csv_safe(args.clone_summary)
    consensus      = load_csv_safe(args.consensus_table)
    variant_matrix = load_csv_safe(args.cross_clone_matrix)

    md_content   = build_markdown(val_report, clone_summary, consensus, variant_matrix, args.pipeline_version)
    html_content = build_html(val_report, clone_summary, consensus, variant_matrix, args.pipeline_version)

    with open(args.out_md, "w") as fh:
        fh.write(md_content)
    log.info(f"Markdown report → {args.out_md}")

    with open(args.out_html, "w") as fh:
        fh.write(html_content)
    log.info(f"HTML report → {args.out_html}")


if __name__ == "__main__":
    main()
