"""
vcf_utils.py — lightweight helpers for VCF comparison and manipulation.

Uses pysam for VCF I/O (via htslib). Avoids loading entire VCFs into memory.
"""

from __future__ import annotations

import gzip
import re
from collections import defaultdict
from pathlib import Path
from typing import Dict, Iterable, Iterator, List, Optional, Set, Tuple

import pandas as pd
import pysam


# ─── Types ───────────────────────────────────────────────────────────────────

# A variant key: (chrom, pos, ref, alt)
VariantKey = Tuple[str, int, str, str]


# ─── VCF reading ─────────────────────────────────────────────────────────────

def iter_pass_variants(
    vcf_path: str | Path,
    pass_only: bool = True,
) -> Iterator[pysam.VariantRecord]:
    """Yield variant records, optionally restricting to PASS filter."""
    with pysam.VariantFile(str(vcf_path)) as vcf:
        for rec in vcf.fetch():
            if pass_only and rec.filter.keys() and "PASS" not in rec.filter:
                continue
            yield rec


def variant_key(rec: pysam.VariantRecord) -> VariantKey:
    """Canonical key for matching variants across callers."""
    return (rec.chrom, rec.pos, rec.ref, ",".join(rec.alts or []))


def vcf_to_variant_set(
    vcf_path: str | Path,
    pass_only: bool = True,
) -> Set[VariantKey]:
    """Return the set of (chrom, pos, ref, alt) tuples in a VCF."""
    return {variant_key(rec) for rec in iter_pass_variants(vcf_path, pass_only)}


def vcf_to_dataframe(
    vcf_path: str | Path,
    extra_fields: Optional[List[str]] = None,
    pass_only: bool = True,
) -> pd.DataFrame:
    """
    Convert a VCF to a tidy DataFrame with columns:
      chrom, pos, ref, alt, qual, filter, [extra INFO/FORMAT fields]
    """
    rows = []
    with pysam.VariantFile(str(vcf_path)) as vcf:
        for rec in vcf.fetch():
            if pass_only and rec.filter.keys() and "PASS" not in rec.filter:
                continue
            row: Dict = {
                "chrom":  rec.chrom,
                "pos":    rec.pos,
                "ref":    rec.ref,
                "alt":    ",".join(rec.alts or []),
                "qual":   rec.qual,
                "filter": ";".join(rec.filter.keys()),
            }
            if extra_fields:
                for field in extra_fields:
                    try:
                        row[field] = rec.info.get(field)
                    except Exception:
                        row[field] = None
            rows.append(row)
    return pd.DataFrame(rows)


# ─── Chromosome name normalisation ───────────────────────────────────────────

_CHR_MAP = {str(i): f"chr{i}" for i in range(1, 23)}
_CHR_MAP.update({"X": "chrX", "Y": "chrY", "M": "chrM", "MT": "chrM"})
_UCSC_TO_ENSEMBL = {v: k for k, v in _CHR_MAP.items()}


def normalize_chrom(chrom: str, target: str = "ucsc") -> str:
    """
    Harmonize chromosome names to UCSC (chr1) or Ensembl (1) style.
    Handles chr1/1, chrM/MT, etc.
    """
    if target == "ucsc":
        return _CHR_MAP.get(chrom, chrom if chrom.startswith("chr") else f"chr{chrom}")
    elif target == "ensembl":
        return _UCSC_TO_ENSEMBL.get(chrom, chrom.lstrip("chr"))
    return chrom


# ─── Multi-VCF comparison ────────────────────────────────────────────────────

def build_presence_absence_matrix(
    vcf_map: Dict[Tuple[str, str], str | Path],
    pass_only: bool = True,
) -> pd.DataFrame:
    """
    Build a variant × (clone, caller) presence/absence matrix.

    Parameters
    ----------
    vcf_map : dict mapping (clone_id, caller) → vcf_path
    pass_only : restrict to PASS variants

    Returns
    -------
    DataFrame indexed by (chrom, pos, ref, alt), columns = (clone_id, caller)
    """
    variant_sets: Dict[Tuple[str, str], Set[VariantKey]] = {}
    for (clone_id, caller), path in vcf_map.items():
        variant_sets[(clone_id, caller)] = vcf_to_variant_set(path, pass_only)

    all_variants = sorted(set().union(*variant_sets.values()))
    cols = sorted(vcf_map.keys())
    data = {col: [int(v in variant_sets[col]) for v in all_variants] for col in cols}

    idx = pd.MultiIndex.from_tuples(all_variants, names=["chrom", "pos", "ref", "alt"])
    df = pd.DataFrame(data, index=idx)
    df.columns = pd.MultiIndex.from_tuples(df.columns, names=["clone_id", "caller"])
    return df


def caller_concordance_per_clone(
    presence_matrix: pd.DataFrame,
) -> pd.DataFrame:
    """
    For each clone, compute how many variants are shared across callers.
    Returns per-clone concordance statistics.
    """
    rows = []
    clones = presence_matrix.columns.get_level_values("clone_id").unique()
    for clone in clones:
        sub = presence_matrix.xs(clone, axis=1, level="clone_id")
        n_callers = sub.shape[1]
        n_any       = (sub.sum(axis=1) > 0).sum()
        n_all       = (sub.sum(axis=1) == n_callers).sum()
        n_two_plus  = (sub.sum(axis=1) >= 2).sum()
        rows.append(
            {
                "clone_id":         clone,
                "n_callers":        n_callers,
                "n_variants_any":   n_any,
                "n_variants_all":   n_all,
                "n_variants_2plus": n_two_plus,
                "concordance_rate": n_all / n_any if n_any > 0 else 0.0,
            }
        )
    return pd.DataFrame(rows)


def classify_variant_sharing(
    presence_matrix: pd.DataFrame,
) -> pd.DataFrame:
    """
    Classify each variant as:
      - shared: detected in all clones (by at least one caller each)
      - partial: detected in >1 but <all clones
      - private: detected in exactly one clone
    """
    clone_presence = (
        presence_matrix.groupby(level="clone_id", axis=1)
        .max()  # any caller in the clone detected it
    )
    n_clones = clone_presence.shape[1]
    clone_counts = clone_presence.sum(axis=1)

    labels = pd.cut(
        clone_counts,
        bins=[0, 1, n_clones - 1, n_clones],
        labels=["private", "partial", "shared"],
        include_lowest=True,
    )
    result = clone_presence.copy()
    result["sharing_class"] = labels
    result["n_clones_with_variant"] = clone_counts
    return result
