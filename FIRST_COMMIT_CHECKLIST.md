# First-Commit Checklist

Generated during a scaffold hardening pass on 2026-05-25.

## What Works

- `bash tests/smoke_test.sh` passes on the bundled `examples/data` files with placeholder BAM paths because it uses `--skip_bam_check`.
- Python syntax check passes with `python3 -m compileall -q bin`.
- Python CLI help loads for `validate_inputs.py`, `parse_medicc2_tree.py`, `assign_clones.py`, `filter_cells.py`, `plot_clone_tree.py`, `compare_variants.py`, `build_consensus.py`, and `generate_report.py`.
- Example BAM manifest, cell metadata, scUnique events, and MEDICC2 tree are internally consistent for the Python validation path: 12 BAM-manifest cells, 12 metadata cells, 22 scUnique events, and 12 tree leaves.
- `schemas/bam_manifest.json` matches the example BAM manifest header: `cell_id,bam_path,bai_path,sample_id,patient_id`.

## What Fails Or Is Not Yet Verified

- `nextflow` is not installed in the current shell, so Nextflow config/script parsing and workflow execution were not verified.
- Strict input validation without `--skip_bam_check` fails on all 12 example BAM paths because they are placeholders under `/path/to/cells`.
- `conf/test.config` points to missing `tests/fixtures` files; only `tests/smoke_test.sh` and `examples/data` exist.
- The documented `nextflow run main.nf -profile test,docker` path is not credible until tiny FASTA/BAM fixtures or a validation-only test profile exist.
- Python environment mismatch: the `python` used by the smoke test has the core smoke dependencies, but it is missing some packages listed in `containers/requirements_python.txt` (`biopython`, `plotly`, `pyyaml`, `rich`) and has versions outside some pins. `python3` is missing most pipeline dependencies.
- Container image names still contain placeholders such as `ghcr.io/TODO/scclone-python:1.0.0`, and `containers/build.sh` defaults to `ghcr.io/TODO_ORG`.
- Matched-normal joining is explicitly not implemented in `subworkflows/mutation_calling.nf`; paired mode currently maps clones through tumor-only-shaped tuples.
- `modules/local/merge_caller_vcfs` is included but unused.
- `modules/local/annotate_vep` is a stub requiring a pre-downloaded cache and unverified flags.

## Fixed In This Pass

- Aliased the three `NORMALIZE_VCF` invocations in `subworkflows/mutation_calling.nf` so the same DSL2 process is not invoked multiple times under one component name.
- Removed unsupported-looking `type: 'dir'` output options from local MultiQC and comparison modules.
- Made optional normal-manifest handling in `VALIDATE_MANIFESTS` truthiness-based instead of assuming a sentinel filename.
- Removed the likely invalid Debian package name `bgzip` from the Python Dockerfile install list; `tabix` provides the command.
- Fixed `compare_variants.py` per-clone TSV writing so it writes variants detected for the current clone/caller instead of reusing the first VCF for every clone.
- Fixed nondeterministic cell ordering in `assign_clones.py` for `event_profile` and `hybrid` strategies by preserving MEDICC2 tree leaf order after BAM-manifest intersection.

## Before Initial Commit

1. Install or provide Nextflow locally, then run at least `nextflow config -profile test` and a no-container parse/lint path if available.
2. Replace `conf/test.config` with real tiny fixtures or rename it to clarify that it is aspirational.
3. Add tiny synthetic BAM/BAI/reference fixtures if the project wants a real Nextflow smoke test; otherwise document that `tests/smoke_test.sh` is Python-only.
4. Replace placeholder container registries and manifest metadata before advertising Docker/Singularity commands.
5. Decide whether initial commit supports tumor-only only; if yes, fail fast when `normal_manifest` is supplied with `tumor_only=false`.
6. Add a minimal test for `compare_variants.py` and `build_consensus.py` using tiny synthetic VCFs.
7. Reconcile README examples with actual runnable files and clearly label commands that require user-provided references/BAMs.
