#!/usr/bin/env python3
"""
aggregate_signatures.py

Combine the per-clone signature exposures into a single cross-clone view:
  - signatures_combined.tsv : clone x signature exposure-fraction matrix
  - signatures_combined.png : stacked-bar of exposure fractions per clone,
                              with low-confidence clones marked (*).

Honest by construction (precision-modeling): a clone whose audit verdict is not
'ok' is still shown, but flagged, so the figure never implies confidence it does
not have.
"""
from __future__ import annotations

import argparse
import logging
import sys

import numpy as np
import pandas as pd

logging.basicConfig(
    level=logging.INFO,
    format="[%(asctime)s] %(levelname)s: %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)
log = logging.getLogger(__name__)


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--exposures", nargs="+", required=True,
                   help="per-clone *.exposures.tsv files")
    p.add_argument("--min_fraction", type=float, default=0.05,
                   help="signatures below this fraction in every clone are pooled as 'Other'")
    p.add_argument("--out_tsv", required=True)
    p.add_argument("--out_png", required=True)
    return p.parse_args()


def load_long(paths) -> pd.DataFrame:
    frames = []
    for p in paths:
        try:
            df = pd.read_csv(p, sep="\t")
        except Exception as e:
            log.warning(f"skip {p}: {e}")
            continue
        if "clone_id" in df.columns and "signature" in df.columns:
            frames.append(df)
    return pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()


def plot_combined(mat: pd.DataFrame, verdict: dict, out_png: str) -> None:
    """Stacked bar of signature fractions per clone. Non-fatal; needs matplotlib."""
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except Exception as e:
        log.warning(f"matplotlib unavailable, skipping plot: {e}")
        return

    clones = list(mat.index)
    sigs = list(mat.columns)
    fig, ax = plt.subplots(figsize=(max(6, len(clones) * 1.3), 6))
    bottom = np.zeros(len(clones))
    cmap = plt.get_cmap("tab20")
    for i, s in enumerate(sigs):
        vals = mat[s].to_numpy(dtype=float)
        ax.bar(range(len(clones)), vals, bottom=bottom, label=s,
               color="#cccccc" if s == "Other" else cmap(i % 20))
        bottom += vals
    ax.set_xticks(range(len(clones)))
    ax.set_xticklabels([f"{c}{'' if verdict.get(c) == 'ok' else ' *'}" for c in clones],
                       rotation=45, ha="right")
    ax.set_ylabel("Signature exposure fraction")
    ax.set_ylim(0, 1)
    ax.set_title("Cross-clone signature exposures  (* = low-confidence verdict)")
    ax.legend(bbox_to_anchor=(1.02, 1), loc="upper left", fontsize=8)
    fig.savefig(out_png, dpi=130, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    args = parse_args()
    long = load_long(args.exposures)
    if long.empty:
        log.warning("No exposures to aggregate — writing empty table")
        pd.DataFrame().to_csv(args.out_tsv, sep="\t", index=False)
        return

    verdict = (long.groupby("clone_id")["clone_label"].first().to_dict()
               if "clone_label" in long.columns else {})

    full = (long.pivot_table(index="clone_id", columns="signature",
                             values="fraction", aggfunc="first")
            .fillna(0.0).sort_index())

    # Keep signatures that are non-trivial in some clone or bootstrap-stable anywhere;
    # pool the rest as 'Other' so each bar still sums to ~1.
    keep = set(full.columns[full.max(axis=0) >= args.min_fraction])
    if "stable" in long.columns:
        keep |= set(long.loc[long["stable"].astype(bool), "signature"].unique())
    keep = [c for c in full.columns if c in keep]
    dropped = [c for c in full.columns if c not in keep]

    mat = full[keep].copy() if keep else full.copy()
    if dropped:
        mat["Other"] = full[dropped].sum(axis=1)

    mat.to_csv(args.out_tsv, sep="\t")
    log.info(f"Combined matrix: {mat.shape[0]} clones x {mat.shape[1]} signatures "
             f"-> {args.out_tsv}")
    plot_combined(mat, verdict, args.out_png)


if __name__ == "__main__":
    main()
