---
name: sc-clone-mutations-scaffold
description: Production Nextflow pipeline scaffold for clone-aware somatic mutation analysis built in this project
metadata:
  type: project
---

Full pipeline scaffold built 2026-05-15 for clone-aware somatic mutation analysis from scDNA-seq (starting from scUnique + MEDICC2 tree + per-cell BAMs).

**Why:** Research-grade pipeline to call somatic mutations per clone, compare across clones, and generate publication-ready QC.

**How to apply:** All code is in `/Users/darvis01/Documents/sc_clone_mutations`. Smoke test passes with `bash tests/smoke_test.sh`.

Key design decisions:
- Default clone strategy is `internal_node` (tree branch length-based)
- Three callers: Mutect2 + Strelka2 + FreeBayes (consensus = ≥2 callers)
- Python container: `ghcr.io/TODO/scclone-python:1.0.0` (Dockerfile.python), needs `TODO_ORG` replaced
- All third-party tools use Biocontainers or official images (GATK, samtools, picard, mosdepth, bcftools)
- Strelka2 tumor-only mode is a known caveat (uses germline workflow)
- scUnique column name auto-detection: TODO verify against user's actual output format
- MEDICC2 events TSV column names: TODO verify `child_node` column name matches user's version

Still needs before real run:
1. Replace `ghcr.io/TODO/scclone-python:1.0.0` with real registry URL in all module main.nf files
2. Run `bash containers/build.sh --push` to build and publish Python container
3. Fill real paths in `examples/params/full_run.yml`
4. Verify scUnique column names match pipeline expectations

[[sc-clone-mutations-user]]
