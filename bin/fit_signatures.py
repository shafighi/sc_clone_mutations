#!/usr/bin/env python3
"""
fit_signatures.py

Per-clone mutational-signature refit with a built-in audit (precision-modeling).

Pipeline:  VCF + reference FASTA -> SBS96 spectrum -> non-negative refit to a fixed
COSMIC SBS reference -> bootstrap stability -> AUDIT -> labelled outputs.

The audit decides the *permitted claim level* for the clone and can return
'insufficient_evidence' / 'poor_reconstruction' / 'confounded' instead of forcing a
signature label. See docs/mutational_signatures.md for the semantic contract.

Usage:
    fit_signatures.py \\
        --clone_id clone_001 \\
        --vcf clone_001.mutect2.norm.vcf.gz \\
        --fasta hsa.GRCh37_g1kp2.fa \\
        --cosmic COSMIC_v3.4_SBS_GRCh37.txt \\
        --min_snv 50 --min_cosine 0.85 --n_boot 200 --seed 0 --tumor_only true \\
        --out_prefix clone_001
"""
from __future__ import annotations

import argparse
import json
import logging
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
from utils import signature_utils as su

logging.basicConfig(
    level=logging.INFO,
    format="[%(asctime)s] %(levelname)s: %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)
log = logging.getLogger(__name__)


def _bool(x: str) -> bool:
    return str(x).lower() in ("1", "true", "yes", "t", "y")


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--clone_id",   required=True)
    p.add_argument("--vcf",        required=True)
    p.add_argument("--fasta",      required=True)
    p.add_argument("--cosmic",     required=True, help="COSMIC SBS reference (TSV/CSV)")
    p.add_argument("--pass_only",  type=_bool, default=True)
    p.add_argument("--tumor_only", type=_bool, default=True)
    p.add_argument("--min_snv",    type=int,   default=50)
    p.add_argument("--min_cosine", type=float, default=0.85)
    p.add_argument("--n_boot",     type=int,   default=200)
    p.add_argument("--seed",       type=int,   default=0)
    p.add_argument("--out_prefix", required=True)
    return p.parse_args()


def plot_signatures(clone_id, counts, exposures, sig_names, audit, out_png, top_n=8):
    """SBS96 spectrum + top signature exposures. Non-fatal; needs matplotlib."""
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except Exception as e:
        log.warning(f"[{clone_id}] matplotlib unavailable, skipping plot: {e}")
        return

    sub_colors = ["#03bcee", "#000000", "#e42926", "#cbcacb", "#a1cf63", "#eccfc6"]
    subs = ["C>A", "C>G", "C>T", "T>A", "T>C", "T>G"]
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 7),
                                   gridspec_kw={"height_ratios": [2, 1]})
    colors = [sub_colors[i // 16] for i in range(96)]
    ax1.bar(range(96), counts, color=colors, width=0.85)
    ax1.set_xlim(-0.5, 95.5)
    ax1.set_ylabel("Mutations")
    ax1.set_xticks([])
    ax1.set_title(f"{clone_id} — SBS96 spectrum (n={audit['n_snv']}, "
                  f"cosine={audit['reconstruction_cosine']:.3f}, "
                  f"verdict={audit['label']})")
    ymax = ax1.get_ylim()[1]
    for i, s in enumerate(subs):
        ax1.text(i * 16 + 7.5, ymax * 0.95, s, ha="center", fontweight="bold")

    order = sorted(range(len(sig_names)), key=lambda k: exposures[k], reverse=True)[:top_n]
    total = float(exposures.sum()) or 1.0
    names = [sig_names[k] for k in order]
    fracs = [float(exposures[k]) / total for k in order]
    stable = {s["signature"]: s["stable"] for s in audit["per_signature"]}
    bar_colors = ["#2c7fb8" if stable.get(n) else "#bdbdbd" for n in names]
    ax2.bar(range(len(names)), fracs, color=bar_colors)
    ax2.set_xticks(range(len(names)))
    ax2.set_xticklabels(names, rotation=45, ha="right")
    ax2.set_ylabel("Exposure fraction")
    ax2.set_title("Top signatures (blue = bootstrap-stable, grey = unstable)")
    plt.tight_layout()
    fig.savefig(out_png, dpi=130)
    plt.close(fig)


def main() -> None:
    args = parse_args()
    channels = su.sbs96_channels()

    log.info(f"[{args.clone_id}] building SBS96 spectrum from {args.vcf}")
    counts, stats = su.build_sbs96_from_vcf(args.vcf, args.fasta, pass_only=args.pass_only)
    n_snv = int(counts.sum())
    log.info(f"[{args.clone_id}] {n_snv} usable SNVs "
             f"(records={stats['records']}, non_pass={stats['non_pass']}, "
             f"non_snv={stats['non_snv']}, multiallelic={stats['multiallelic']}, "
             f"context_skipped={stats['context_skipped']})")

    pd.DataFrame({"channel": channels, "count": counts.astype(int)}).to_csv(
        f"{args.out_prefix}.sbs96_counts.tsv", sep="\t", index=False)

    sig_names, W = su.load_signature_reference(args.cosmic)
    log.info(f"[{args.clone_id}] reference: {len(sig_names)} COSMIC SBS signatures")

    exposures = su.refit_exposures(counts, W, seed=args.seed)
    boot = su.bootstrap_exposures(counts, W, n_boot=args.n_boot, seed=args.seed)
    audit = su.audit_fit(
        counts, W, sig_names, exposures, boot,
        min_snv=args.min_snv, min_cosine=args.min_cosine,
        tumor_only=args.tumor_only,
    )
    audit["clone_id"] = args.clone_id
    audit["vcf_stats"] = stats
    audit["config"] = {
        "vcf": args.vcf, "fasta": args.fasta, "cosmic": args.cosmic,
        "pass_only": args.pass_only, "tumor_only": args.tumor_only,
        "min_snv": args.min_snv, "min_cosine": args.min_cosine,
        "n_boot": args.n_boot, "seed": args.seed,
        "reference_signatures": len(sig_names),
    }

    # Tidy per-signature table (sorted by exposure), + the clone-level verdict.
    rows = sorted(audit["per_signature"], key=lambda r: r["exposure"], reverse=True)
    df = pd.DataFrame(rows)
    df.insert(0, "clone_id", args.clone_id)
    df["clone_label"] = audit["label"]
    df.to_csv(f"{args.out_prefix}.exposures.tsv", sep="\t", index=False)

    boot_df = pd.DataFrame({
        "clone_id": args.clone_id,
        "signature": sig_names,
        "boot_mean_frac": boot["mean"],
        "boot_lo": boot["lo"],
        "boot_hi": boot["hi"],
    })
    boot_df.to_csv(f"{args.out_prefix}.exposures_bootstrap.tsv", sep="\t", index=False)

    with open(f"{args.out_prefix}.signatures_audit.json", "w") as fh:
        json.dump(audit, fh, indent=2)

    plot_signatures(args.clone_id, counts, exposures, sig_names, audit,
                    f"{args.out_prefix}.signatures.png")

    # Loud, honest verdict — never bury a failure label.
    verdict = audit["label"]
    cos = audit["reconstruction_cosine"]
    if verdict == "ok":
        top = [r for r in rows if r["stable"]][:5]
        names = ", ".join(f"{r['signature']}({r['fraction']:.0%})" for r in top) or "none stable"
        log.info(f"[{args.clone_id}] VERDICT=ok  n_snv={n_snv}  cosine={cos:.3f}  "
                 f"stable: {names}")
    else:
        log.warning(f"[{args.clone_id}] VERDICT={verdict}  n_snv={n_snv}  cosine={cos:.3f} "
                    f"-> exposures reported but NOT actionable (see docs/mutational_signatures.md)")


if __name__ == "__main__":
    main()
