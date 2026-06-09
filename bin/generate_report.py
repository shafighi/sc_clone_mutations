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
import os
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
    p.add_argument("--clone_reliability",  default=None)
    p.add_argument("--clone_tree_png",     default=None)
    p.add_argument("--signatures",         nargs="*", default=[],
                   help="per-clone *.exposures.tsv files")
    p.add_argument("--signature_plots",    nargs="*", default=[],
                   help="per-clone *.signatures.png files")
    p.add_argument("--pipeline_version",   default="unknown")
    p.add_argument("--out_html",           required=True)
    p.add_argument("--out_md",             required=True)
    return p.parse_args()


def img_data_uri(path: Optional[str]) -> Optional[str]:
    """Read a PNG and return a base64 data URI so the HTML is self-contained."""
    if not path or not os.path.exists(path):
        return None
    try:
        import base64
        with open(path, "rb") as fh:
            b64 = base64.b64encode(fh.read()).decode("ascii")
        return f"data:image/png;base64,{b64}"
    except Exception as e:
        log.warning(f"Could not embed image {path}: {e}")
        return None


def aggregate_signatures(paths: List[str]) -> pd.DataFrame:
    """
    Collapse per-clone *.exposures.tsv into one tidy row per clone:
    clone_id, verdict (clone_label), and the top bootstrap-stable signatures.
    """
    rows = []
    for p in paths or []:
        try:
            df = pd.read_csv(p, sep="\t")   # exposures files are tab-separated
        except Exception as e:
            log.warning(f"Could not load {p}: {e}")
            continue
        if df.empty or "clone_id" not in df.columns:
            continue
        clone = str(df["clone_id"].iloc[0])
        verdict = str(df["clone_label"].iloc[0]) if "clone_label" in df.columns else "n/a"
        stable = df[df["stable"]] if "stable" in df.columns else df
        stable = stable.sort_values("fraction", ascending=False) if "fraction" in stable.columns else stable
        top = "; ".join(
            f"{r['signature']} ({r['fraction']:.0%})"
            for _, r in stable.head(5).iterrows()
        ) if not stable.empty else "none stable"
        rows.append({"clone_id": clone, "verdict": verdict, "top_signatures": top})
    return pd.DataFrame(rows).sort_values("clone_id") if rows else pd.DataFrame()


def clone_banner(clone_summary: pd.DataFrame,
                 reliability: Optional[pd.DataFrame]) -> tuple:
    """
    Loud, top-of-report status on how many clones were resolved/called.
    Returns (level, icon, message). level in {danger, warn, ok}.
    The danger case (one clone) is the 'where did my clones go?' guard.
    """
    n_clones = len(clone_summary) if (clone_summary is not None and not clone_summary.empty) else 0
    has_rel  = reliability is not None and not reliability.empty
    n_called   = (int(reliability["called"].sum())
                  if has_rel and "called" in reliability.columns else n_clones)
    n_reliable = (int(reliability["reliable"].sum())
                  if has_rel and "reliable" in reliability.columns else None)

    if n_clones <= 1:
        return ("danger", "⛔",
                "Only ONE clone was resolved — no subclonal structure, so there is "
                "nothing to compare across clones. If you expected subclones, then no "
                "tree branch carried enough copy-number events to split them: lower "
                "min_event_count (the clone-boundary threshold) or inspect the tree. "
                "See clone_summary.csv.")
    if n_reliable == 0:
        return ("warn", "⚠",
                f"{n_clones} clones resolved and {n_called} called — but NONE meet the "
                f"reliable cell-count threshold, so every call is exploratory "
                f"(low pseudobulk depth). See the reliability table below.")
    tail = f" · {n_reliable} reliable" if n_reliable is not None else ""
    return ("ok", "✅",
            f"{n_clones} clones resolved · {n_called} called{tail}. "
            f"Small clones are called and flagged below.")


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
    # Build the table by hand rather than via DataFrame.to_markdown(), which
    # requires the optional 'tabulate' dependency that may be absent.
    sub  = df.head(max_rows)
    cols = [str(c) for c in sub.columns]
    header = "| " + " | ".join(cols) + " |"
    sep    = "| " + " | ".join("---" for _ in cols) + " |"
    rows   = [
        "| " + " | ".join("" if pd.isna(v) else str(v) for v in row) + " |"
        for row in sub.itertuples(index=False, name=None)
    ]
    return "\n".join([header, sep, *rows]) + "\n"


def build_markdown(
    val_report: Dict,
    clone_summary: pd.DataFrame,
    consensus_table: pd.DataFrame,
    variant_matrix: pd.DataFrame,
    version: str,
    reliability: Optional[pd.DataFrame] = None,
    signatures: Optional[pd.DataFrame] = None,
) -> str:
    ts = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    lines: List[str] = []

    lines.append(f"# sc_clone_mutations Report")
    lines.append(f"\n**Pipeline version:** {version}  ")
    lines.append(f"**Generated:** {ts}\n")

    # ── Clone banner (loud 'where did my clones go?' status) ─────────────────
    _level, _icon, _msg = clone_banner(clone_summary, reliability)
    lines.append(f"> {_icon} **Clones:** {_msg}\n")

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

    # ── Clone reliability for mutation calling ───────────────────────────────
    lines.append("## Clone Reliability for Mutation Calling")
    if reliability is not None and not reliability.empty:
        lines.append(df_to_md(reliability))
        lines.append("_Clones below the reliable threshold are still called but their "
                     "variant lists are exploratory (low pseudobulk depth; tumor-only "
                     "artefacts/germline not subtracted)._\n")
    else:
        lines.append("_Reliability table not available._\n")

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

    # ── Mutational signatures ────────────────────────────────────────────────
    lines.append("## Mutational Signatures (per clone)")
    if signatures is not None and not signatures.empty:
        lines.append(df_to_md(signatures))
        lines.append("_'verdict' is the audit label; only bootstrap-stable signatures "
                     "are listed. See the HTML report for per-clone SBS96 spectra and "
                     "exposure plots._\n")
    else:
        lines.append("_Signatures not run (enable with --run_signatures and "
                     "--cosmic_signatures)._\n")

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

{clone_banner}

<h2>1. Input Validation</h2>
<p>Status: <span class="badge {badge_class}">{status}</span></p>
{validation_details}
{validation_errors}

<h2>2. Clone Assignment</h2>
{clone_table}

<h2>Clone Tree</h2>
{tree_section}

<h2>Clone Reliability for Mutation Calling</h2>
{reliability_section}

<h2>3. Consensus Mutations</h2>
{consensus_table}

<h2>4. Cross-Clone Variant Sharing</h2>
{sharing_section}

<h2>Mutational Signatures</h2>
{signatures_section}

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
    reliability: Optional[pd.DataFrame] = None,
    tree_png_uri: Optional[str] = None,
    signatures: Optional[pd.DataFrame] = None,
    signature_plot_uris: Optional[List[str]] = None,
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

    _banner_colors = {"danger": ("#f8d7da", "#dc3545"),
                      "warn":   ("#fff3cd", "#ffc107"),
                      "ok":     ("#d1e7dd", "#198754")}
    level, icon, msg = clone_banner(clone_summary, reliability)
    bg, border = _banner_colors[level]
    clone_banner_html = (
        f'<div style="background:{bg};border-left:6px solid {border};'
        f'padding:0.9rem 1.1rem;margin:1.2rem 0;border-radius:4px;font-size:1.05rem;">'
        f'{icon} <strong>Clones:</strong> {msg}</div>'
    )

    tree_section = (f'<img src="{tree_png_uri}" style="max-width:100%;height:auto;">'
                    if tree_png_uri else "<p>Tree image not available.</p>")

    if reliability is not None and not reliability.empty:
        reliability_section = (
            df_to_html(reliability) +
            "<p class='text-muted'>Clones below the reliable threshold are still "
            "called, but their variant lists are exploratory (low pseudobulk depth; "
            "tumor-only artefacts/germline not subtracted).</p>"
        )
    else:
        reliability_section = "<p>Not available.</p>"

    if signatures is not None and not signatures.empty:
        imgs = "".join(
            f'<div style="margin:1rem 0;"><img src="{u}" '
            f'style="max-width:100%;height:auto;"></div>'
            for u in (signature_plot_uris or [])
        )
        signatures_section = df_to_html(signatures) + imgs
    else:
        signatures_section = ("<p>Signatures not run (enable with "
                              "<code>--run_signatures</code> and "
                              "<code>--cosmic_signatures</code>).</p>")

    return HTML_TEMPLATE.format(
        version=version,
        timestamp=ts,
        status=status,
        badge_class=badge_class,
        validation_details=validation_details,
        validation_errors=validation_errors,
        clone_banner=clone_banner_html,
        clone_table=clone_table_html,
        tree_section=tree_section,
        reliability_section=reliability_section,
        consensus_table=consensus_html,
        sharing_section=sharing_section,
        signatures_section=signatures_section,
        warnings_section=warnings_section,
    )


def main() -> None:
    args = parse_args()

    val_report     = load_json(args.validation_report)
    clone_summary  = load_csv_safe(args.clone_summary)
    consensus      = load_csv_safe(args.consensus_table)
    variant_matrix = load_csv_safe(args.cross_clone_matrix)
    reliability    = load_csv_safe(args.clone_reliability) if args.clone_reliability else pd.DataFrame()
    signatures     = aggregate_signatures(args.signatures)
    tree_png_uri   = img_data_uri(args.clone_tree_png)
    sig_plot_uris  = [u for u in (img_data_uri(p) for p in args.signature_plots) if u]

    md_content   = build_markdown(
        val_report, clone_summary, consensus, variant_matrix, args.pipeline_version,
        reliability=reliability, signatures=signatures)
    html_content = build_html(
        val_report, clone_summary, consensus, variant_matrix, args.pipeline_version,
        reliability=reliability, tree_png_uri=tree_png_uri,
        signatures=signatures, signature_plot_uris=sig_plot_uris)

    with open(args.out_md, "w") as fh:
        fh.write(md_content)
    log.info(f"Markdown report → {args.out_md}")

    with open(args.out_html, "w") as fh:
        fh.write(html_content)
    log.info(f"HTML report → {args.out_html}")


if __name__ == "__main__":
    main()
