"""
signature_utils.py — SBS96 construction and audited COSMIC signature refitting.

Design follows the precision-modeling skill: the math is transparent (no black-box
refit), every biological claim has an audit, and the audit can return
"insufficient evidence" rather than forcing a signature label.

The numerical core (refit, cosine, bootstrap, audits) depends only on numpy so it
can be unit-tested on synthetic data without a reference genome, pysam, or COSMIC.
VCF/FASTA parsing imports pysam lazily, so importing this module never requires it.
"""

from __future__ import annotations

import logging
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np

log = logging.getLogger(__name__)

_BASES = ("A", "C", "G", "T")
_SUBS = ("C>A", "C>G", "C>T", "T>A", "T>C", "T>G")
_COMP = {"A": "T", "C": "G", "G": "C", "T": "A", "N": "N"}

# COSMIC v3 "possible sequencing artefact" signatures — exposure landing here is a
# confounder flag in tumor-only data, not biology. Conservative default set.
ARTIFACT_SBS = frozenset(
    ["SBS27", "SBS43", "SBS45", "SBS46", "SBS47", "SBS48", "SBS49", "SBS50",
     "SBS51", "SBS52", "SBS53", "SBS54", "SBS55", "SBS56", "SBS57", "SBS58",
     "SBS59", "SBS60", "SBS95"]
)
CLOCK_SBS = frozenset(["SBS1", "SBS5"])


def sbs96_channels() -> List[str]:
    """The 96 SBS channels in canonical COSMIC order, e.g. 'A[C>A]A'."""
    return [f"{five}[{sub}]{three}"
            for sub in _SUBS for five in _BASES for three in _BASES]


def revcomp(seq: str) -> str:
    return "".join(_COMP.get(b, "N") for b in reversed(seq))


def classify_snv(ref: str, alt: str, triplet: str) -> Optional[str]:
    """
    Map an SNV + its 5'-ref-3' trinucleotide to an SBS96 channel label.

    Returns None for non-SNVs, ambiguous bases, or context whose centre base does
    not equal *ref* (a guard against mis-aligned FASTA lookups).
    """
    if len(ref) != 1 or len(alt) != 1 or ref == alt:
        return None
    ref, alt, triplet = ref.upper(), alt.upper(), triplet.upper()
    if len(triplet) != 3 or any(b not in _BASES for b in triplet):
        return None
    if triplet[1] != ref:
        return None
    # Fold to pyrimidine reference (C/T) so each mutation maps to one channel.
    if ref in ("A", "G"):
        ref, alt, triplet = _COMP[ref], _COMP[alt], revcomp(triplet)
    sub = f"{ref}>{alt}"
    if sub not in _SUBS:
        return None
    return f"{triplet[0]}[{sub}]{triplet[2]}"


def build_sbs96_from_vcf(
    vcf_path: str,
    fasta_path: str,
    pass_only: bool = True,
) -> Tuple[np.ndarray, Dict[str, int]]:
    """
    Build a 96-channel SBS count vector from a VCF + reference FASTA.

    Imports pysam lazily (only needed for the real pipeline, not unit tests).
    Returns (counts[96], stats) where stats records how records were handled.
    """
    import pysam  # lazy

    channels = sbs96_channels()
    index = {c: i for i, c in enumerate(channels)}
    counts = np.zeros(96, dtype=float)
    stats = {"records": 0, "snv": 0, "non_snv": 0, "non_pass": 0,
             "multiallelic": 0, "context_skipped": 0, "counted": 0}

    fa = pysam.FastaFile(fasta_path)
    fa_contigs = set(fa.references)
    vcf = pysam.VariantFile(vcf_path)
    for rec in vcf:
        stats["records"] += 1
        if pass_only and rec.filter.keys() not in ([], ["PASS"]):
            stats["non_pass"] += 1
            continue
        alts = rec.alts or ()
        if len(alts) != 1:
            stats["multiallelic"] += 1
            continue
        ref, alt = rec.ref or "", alts[0] or ""
        if len(ref) != 1 or len(alt) != 1:
            stats["non_snv"] += 1
            continue
        stats["snv"] += 1
        contig = rec.contig
        if contig not in fa_contigs:
            # tolerate chr-prefix mismatch
            alt_contig = contig[3:] if contig.startswith("chr") else "chr" + contig
            contig = alt_contig if alt_contig in fa_contigs else None
        if contig is None:
            stats["context_skipped"] += 1
            continue
        pos0 = rec.pos - 1  # pysam.pos is 1-based; fetch is 0-based half-open
        try:
            triplet = fa.fetch(contig, pos0 - 1, pos0 + 2)
        except (ValueError, KeyError):
            stats["context_skipped"] += 1
            continue
        channel = classify_snv(ref, alt, triplet)
        if channel is None:
            stats["context_skipped"] += 1
            continue
        counts[index[channel]] += 1.0
        stats["counted"] += 1
    return counts, stats


def load_signature_reference(path: str) -> Tuple[List[str], np.ndarray]:
    """
    Load a COSMIC SBS reference matrix (TSV/CSV) into (signature_names, W[96, K]).

    The channel column may be named Type/MutationType/MutationsType or be the first
    column. Rows are reordered to canonical SBS96 order and each signature column is
    normalised to sum to 1 (a probability distribution over the 96 channels).
    """
    import csv

    with open(path, newline="") as fh:
        sniff = fh.readline()
        delim = "\t" if "\t" in sniff else ","
        fh.seek(0)
        reader = csv.reader(fh, delimiter=delim)
        header = next(reader)
        rows = [r for r in reader if r and any(c.strip() for c in r)]

    chan_col = 0
    for i, name in enumerate(header):
        if name.strip().lower() in ("type", "mutationtype", "mutationstype", "channel"):
            chan_col = i
            break
    sig_names = [h.strip() for j, h in enumerate(header) if j != chan_col]

    table: Dict[str, List[float]] = {}
    for r in rows:
        ch = r[chan_col].strip()
        table[ch] = [float(v) for j, v in enumerate(r) if j != chan_col]

    channels = sbs96_channels()
    missing = [c for c in channels if c not in table]
    if missing:
        raise ValueError(
            f"Reference is missing {len(missing)} SBS96 channels (e.g. {missing[:3]}); "
            f"check that channel labels use the 'A[C>A]A' convention."
        )
    W = np.array([table[c] for c in channels], dtype=float)  # 96 x K
    colsum = W.sum(axis=0)
    colsum[colsum == 0] = 1.0
    W = W / colsum
    return sig_names, W


def refit_exposures(
    counts: np.ndarray,
    W: np.ndarray,
    max_iter: int = 10000,
    tol: float = 1e-7,
    seed: int = 0,
) -> np.ndarray:
    """
    Refit non-negative exposures h (length K) so that W·h approximates *counts*,
    minimising generalised KL divergence (the right objective for Poisson counts).

    Lee–Seung multiplicative updates with W fixed. Deterministic given *seed*.
    Returns exposures scaled so they sum to the total mutation count.
    """
    counts = np.asarray(counts, dtype=float)
    total = counts.sum()
    if total <= 0:
        return np.zeros(W.shape[1], dtype=float)

    rng = np.random.default_rng(seed)
    h = rng.uniform(0.5, 1.5, size=W.shape[1]) * (total / W.shape[1])
    wcol = W.sum(axis=0)               # per-signature column mass (≈1 after norm)
    wcol[wcol == 0] = 1e-12
    prev = h.copy()
    for it in range(max_iter):
        recon = W @ h + 1e-12
        h = h * ((W.T @ (counts / recon)) / wcol)
        if it % 25 == 0:
            if np.max(np.abs(h - prev)) / (np.max(prev) + 1e-12) < tol:
                break
            prev = h.copy()
    # Rescale to conserve total mutation mass (interpretable as #mutations/signature)
    s = h.sum()
    if s > 0:
        h = h * (total / s)
    return h


def reconstruct(h: np.ndarray, W: np.ndarray) -> np.ndarray:
    return W @ np.asarray(h, dtype=float)


def cosine_similarity(a: np.ndarray, b: np.ndarray) -> float:
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    na, nb = np.linalg.norm(a), np.linalg.norm(b)
    if na == 0 or nb == 0:
        return 0.0
    return float(np.dot(a, b) / (na * nb))


def bootstrap_exposures(
    counts: np.ndarray,
    W: np.ndarray,
    n_boot: int = 200,
    seed: int = 0,
    ci: float = 0.95,
) -> Dict[str, np.ndarray]:
    """
    Multinomial bootstrap over mutations: resample the spectrum n_boot times, refit,
    and summarise per-signature exposure fraction. Returns mean and CI arrays (len K).
    A signature whose CI lower bound is ~0 is not robustly present.
    """
    counts = np.asarray(counts, dtype=float)
    total = int(round(counts.sum()))
    K = W.shape[1]
    if total <= 0:
        z = np.zeros(K)
        return {"mean": z, "lo": z.copy(), "hi": z.copy()}
    p = counts / counts.sum()
    rng = np.random.default_rng(seed)
    fracs = np.zeros((n_boot, K), dtype=float)
    for b in range(n_boot):
        resampled = rng.multinomial(total, p).astype(float)
        h = refit_exposures(resampled, W, seed=seed + b + 1)
        s = h.sum()
        fracs[b] = h / s if s > 0 else h
    lo_q, hi_q = (1 - ci) / 2, 1 - (1 - ci) / 2
    return {
        "mean": fracs.mean(axis=0),
        "lo": np.quantile(fracs, lo_q, axis=0),
        "hi": np.quantile(fracs, hi_q, axis=0),
    }


def audit_fit(
    counts: np.ndarray,
    W: np.ndarray,
    sig_names: Sequence[str],
    exposures: np.ndarray,
    boot: Dict[str, np.ndarray],
    *,
    min_snv: int = 50,
    min_cosine: float = 0.85,
    min_stable_frac: float = 0.01,
    tumor_only: bool = True,
    artifact_sbs: frozenset = ARTIFACT_SBS,
) -> Dict[str, object]:
    """
    Turn a fit into an auditable verdict. Returns a JSON-serialisable dict including
    the overall label, per-signature stability calls, and confounder flags.

    The overall label is the model's permitted claim level:
      insufficient_evidence < poor_reconstruction < confounded < ok
    """
    counts = np.asarray(counts, dtype=float)
    n_snv = int(round(counts.sum()))
    recon = reconstruct(exposures, W)
    cos = cosine_similarity(counts, recon)

    total_exp = exposures.sum()
    frac = exposures / total_exp if total_exp > 0 else exposures
    lo = boot["lo"]

    per_sig = []
    artifact_share = 0.0
    clock_share = 0.0
    for k, name in enumerate(sig_names):
        stable = bool(lo[k] > min_stable_frac)
        per_sig.append({
            "signature": name,
            "exposure": float(exposures[k]),
            "fraction": float(frac[k]),
            "boot_mean_frac": float(boot["mean"][k]),
            "boot_lo": float(boot["lo"][k]),
            "boot_hi": float(boot["hi"][k]),
            "stable": stable,
        })
        if name in artifact_sbs:
            artifact_share += float(frac[k])
        if name in CLOCK_SBS:
            clock_share += float(frac[k])

    if n_snv < min_snv:
        label = "insufficient_evidence"
    elif cos < min_cosine:
        label = "poor_reconstruction"
    elif tumor_only and artifact_share >= 0.5:
        label = "confounded"
    else:
        label = "ok"

    return {
        "label": label,
        "n_snv": n_snv,
        "reconstruction_cosine": cos,
        "artifact_exposure_fraction": artifact_share,
        "clocklike_exposure_fraction": clock_share,
        "tumor_only": bool(tumor_only),
        "gates": {"min_snv": min_snv, "min_cosine": min_cosine,
                  "min_stable_frac": min_stable_frac},
        "n_stable_signatures": int(sum(1 for s in per_sig if s["stable"])),
        "per_signature": per_sig,
    }
