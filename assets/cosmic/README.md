# COSMIC SBS reference signatures

The mutational-signatures step (`--run_signatures`) refits each clone's SBS96
spectrum to a **fixed** COSMIC SBS reference. That reference is **not** bundled
(licence + size); provide it with `--cosmic_signatures <path>`.

## Get the file

Download the single-base-substitution (SBS) reference from COSMIC, matching your
reference build. This pipeline calls variants on **GRCh37**, so use the GRCh37 SBS
file (v3.4 recommended):

- https://cancer.sanger.ac.uk/signatures/downloads/  → "SBS" → GRCh37

Example (v3.4):
```
COSMIC_v3.4_SBS_GRCh37.txt
```

## Expected format

A TSV (or CSV) with one row per SBS96 channel and one column per signature:

```
Type        SBS1      SBS2      SBS3    ...
A[C>A]A     0.000886  0.000234  0.0201  ...
A[C>A]C     0.000228  0.000045  0.0166  ...
...
T[T>G]T     0.000392  0.000012  0.0095  ...
```

- The channel column may be named `Type`, `MutationType`, `MutationsType`, or
  `channel` (or simply be the first column).
- Channel labels must use the `5'[REF>ALT]3'` convention (e.g. `A[C>A]A`).
- All 96 channels must be present. Columns are renormalised to sum to 1 on load,
  and rows are reordered to the canonical SBS96 order automatically.

## Run

```
nextflow run main.nf -profile slurm \
    --scunique_dir ... \
    --run_signatures \
    --cosmic_signatures /path/to/COSMIC_v3.4_SBS_GRCh37.txt
```

Outputs land in `results/signatures/<clone>/`. Read the verdict in
`<clone>.signatures_audit.json` (`label` field) — and see
`docs/mutational_signatures.md` for what each label permits you to claim.
