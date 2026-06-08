# Mutational signatures (per-clone SBS96 refit)

Follows the `precision-modeling` skill: the biological claim is part of the model
contract, every claim is audited, and the model is allowed to answer
**"insufficient evidence"** instead of being forced to emit a signature label.

## Semantic contract

**What is learned.** For one clone, a vector of non-negative **exposures** over a
*fixed* reference set of COSMIC SBS signatures, such that
`exposures · reference ≈ observed SBS96 spectrum`.

**What the model MAY claim.** "This clone's somatic SNV spectrum is consistent
with a mixture of COSMIC signatures {Sᵢ} in these proportions, at reconstruction
cosine = r over n_snv mutations." Nothing more.

**What the model MUST NOT claim.** That a signature is biologically *active* or
*causal*. A non-zero exposure is necessary but not sufficient for that claim; it
also requires (a) sufficient mutation burden, (b) good reconstruction, and
(c) bootstrap stability. Tumor-only calling means the spectrum can contain
germline/artifact mutations, so etiology cannot be asserted without that caveat.

**Null interpretation (failure label).** If any gate fails, the clone is labelled
`insufficient_evidence` / `poor_reconstruction` / `confounded` and the per-signature
exposures are reported as *non-actionable* — never surfaced as findings.

## Model–meaning table

| Claim | Model object | Enforced by | Audited by | Failure label |
|-------|--------------|-------------|------------|---------------|
| A signature mixture explains the spectrum | exposures (≥0) | non-negativity in refit | reconstruction cosine vs observed | `poor_reconstruction` |
| There is enough signal to fit | mutation count | min-burden gate | `n_snv` | `insufficient_evidence` |
| A signature is robustly present | per-signature exposure | bootstrap resampling | bootstrap CI lower bound > 0 | `unstable` |
| Exposures reflect somatic biology | SBS96 spectrum | tumor-only / PoN caveat + artifact-signature audit | `tumor_only` flag, artifact-SBS share | `confounded` |

## Nonsense checks (pre-coding)

- **Low burden:** below `min_snv` (default 50) mutations, refit is noise → `insufficient_evidence`.
- **Reconstruction:** cosine(observed, reconstructed) below `min_cosine` (default 0.85) → `poor_reconstruction`.
- **Stability:** exposures whose bootstrap lower CI ≈ 0 are not asserted (`unstable`).
- **Confounder:** tumor-only input → flag; report the fraction of exposure landing on
  COSMIC "possible-sequencing-artifact" signatures (SBS27/43/45–60 by default).
- **Collapse to clock:** a fit dominated by SBS1/SBS5 (clock-like) is the *null*; it is
  reported but explicitly not over-interpreted.

## Outputs (auditable, under `results/signatures/<clone>/`)

- `sbs96_counts.tsv` — the 96-channel spectrum.
- `exposures.tsv` — per-signature exposure (counts and fraction).
- `exposures_bootstrap.tsv` — per-signature bootstrap mean + 95% CI.
- `audit.json` — n_snv, cosine, gates, labels, artifact share, config echo, seed.
- `signatures_report.tsv` — one tidy row per (clone, signature) with the stability call.

All thresholds, seeds, the reference file and the bootstrap count are explicit in
`nextflow.config` and echoed into `audit.json`. Generated files live under `results/`
and are never hand-edited.

## Reference data

COSMIC SBS reference is **not** bundled. Provide it via `--cosmic_signatures`
(GRCh37/v3.4 recommended). See `assets/cosmic/README.md` for how to obtain it.
The loader validates and reorders rows to the canonical SBS96 channel order, so any
COSMIC release with a `Type`/`MutationType` channel column works.
