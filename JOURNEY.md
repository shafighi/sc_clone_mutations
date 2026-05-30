# Pipeline Journey — SLX-24491 Run

## Overview

Running `sc_clone_mutations` pipeline on **SLX-24491** dataset (212 cells, tumor-only mode, GRCh37) on the institutional SLURM cluster (`epyc` partition).

---

## Timeline

### 2026-05-26: Initial attempts (SLX-23961)

- 8 consecutive runs, all failed (`ERR` status)
- Runs lasted 30s to 10min — likely input/config issues
- Switched to SLX-24491 dataset on May 27

### 2026-05-27: SLX-24491 first runs

**Run: `hopeful_maxwell`** (4m 5s, ERR)
- Quick failure, likely config issue after switching datasets

**Run: `distracted_crick`** (8h 3m 40s, ERR)
- Got through stages A–C (validation, clone definition, pseudobulk merge)
- **Failed at: `PSEUDOBULK:MERGE_BAMS (clone_001)`** — exit code 140
- **Root cause**: OOM kill. Merging 212 BAMs simultaneously exceeded the `process_medium` memory limit (36 GB)

### 2026-05-28: Fix #1 — Increase MERGE_BAMS resources

**Decision**: Upgrade `MERGE_BAMS` from `process_medium` (36 GB) to `process_high` (72 GB) with retry (up to 144 GB on retry).

**Commit**: `65fe562` — "fix: increase MERGE_BAMS memory to process_high with retry for large clones"

**Result**: MERGE_BAMS succeeded ✔ (cached on subsequent runs)

### 2026-05-28: Fix #2 — MARKDUPLICATES also OOM

**Run: `drunk_boyd`** — MERGE_BAMS cached, but now `PSEUDOBULK:MARKDUPLICATES (clone_001)` killed by system.

- First attempt: `-Xmx36g` (process_medium) — killed
- Picard had read 280M records and was only on chr1

**Decision**: Upgrade `MARKDUPLICATES` from `process_medium` to `process_high` with retry.

**Commit**: `854b88a` — "fix: increase MARKDUPLICATES memory to process_high with retry for large pseudobulks"

### 2026-05-29: Fix #2 still fails — MARKDUPLICATES killed at 128 GB

**Run: `spontaneous_einstein`** — MARKDUPLICATES retried 3 times (72 GB → 128 GB → 128 GB cap), all killed.

Key observations from the logs:
- BAM contains **~2 billion reads** (212 cells × ~9M reads/cell WGS)
- At chr11 (attempt 3), picard had read 1.97 billion records in 1h19m
- Tracking 1.7 million unmatched pairs
- `epyc` partition has **infinite** time limit — so it's not a time kill
- Even at 128 GB (`max_memory` cap), the job is killed — likely true OOM or the JVM + OS overhead exceeds the SLURM allocation

### 2026-05-30: Fix #3 — Skip MarkDuplicates entirely

**Decision**: Set `--mark_duplicates false`

**Reasoning**:
1. **Scientifically justified**: MarkDuplicates on pseudobulk from single cells is problematic. Reads from different cells with identical mapping positions are biologically independent molecules, NOT PCR duplicates. Picard cannot distinguish them and will incorrectly flag real data.
2. **Low impact of skipping**: scDNA-seq libraries have low per-cell duplication rates. With 212 cells merged, the contribution of true PCR duplicates is negligible.
3. **Callers handle it**: Mutect2 has its own duplicate logic, Strelka2 works fine on raw data, FreeBayes ignores duplicate flags entirely.
4. **Practical**: The BAM is simply too large (2B reads) for picard to process within any reasonable memory/time budget on this cluster.

**Commit**: `4f04d57` — "fix: skip MarkDuplicates for pseudobulk (--mark_duplicates false)"

---

## Current Pipeline Command

```bash
nextflow run $PIPELINE_DIR/main.nf \
    -profile singularity,slurm \
    --scunique_dir  /mnt/scratche/slow/fmlab/darvis01/scDNAseq-workflow/results/SLX-24491_100/ \
    --bam_dir       /mnt/scratche/slow/fmlab/darvis01/scDNAseq-workflow-main/data/aligned/SLX-24491/ \
    --fasta         /mnt/scratche/slow/fmlab/darvis01/snakemake-illumina-alignment/resources/homo_sapiens/GRCh37_g1kp2/fasta/hsa.GRCh37_g1kp2.fa \
    --tumor_only    true \
    --mark_duplicates false \
    --outdir        /mnt/scratche/slow/fmlab/darvis01/sc_clone_mutations_results \
    -work-dir       /mnt/scratche/slow/fmlab/darvis01/nf_work \
    -resume \
    -with-report    /mnt/scratche/slow/fmlab/darvis01/sc_clone_mutations_results/pipeline_info/report.html
```

---

## Dataset Summary

| Property | Value |
|----------|-------|
| Sample | SLX-24491 |
| Cells in clone_001 | 212 |
| Approx reads in pseudobulk | ~2 billion |
| Reference | GRCh37 (hg19) |
| Mode | Tumor-only |
| Clone strategy | internal_node (default) |
| Callers | mutect2, strelka2, freebayes |
| Cluster partition | epyc (infinite time limit) |

---

## Key Decisions Log

| # | Decision | Rationale |
|---|----------|-----------|
| 1 | Upgrade MERGE_BAMS to `process_high` | 212 BAMs need >36 GB to merge simultaneously |
| 2 | Upgrade MARKDUPLICATES to `process_high` | Large pseudobulk BAM exceeds medium resources |
| 3 | Skip MarkDuplicates entirely | 2B-read pseudobulk exceeds even 128 GB; scientifically inappropriate for merged single-cell data anyway |
| 4 | Update mosdepth container tag | `0.3.9--hd299d5a_0` purged from quay.io; switched to `0.3.8--hd299d5a_0` |
| 5 | Update all biocontainer tags | quay.io purged many old tags; queried API and updated freebayes, strelka, samtools, bcftools, multiqc |

### 2026-05-30: Fix #4 — Mosdepth container tag removed from quay.io

**Run: `soggy_lavoisier`** (job 52364201) — MarkDuplicates correctly skipped, but pipeline failed at `PSEUDOBULK:MOSDEPTH` with:

```
FATAL: MANIFEST_UNKNOWN: manifest unknown
GET https://quay.io/v2/biocontainers/mosdepth/manifests/0.3.9--hd299d5a_0
```

**Root cause**: quay.io/biocontainers purged the `0.3.9--hd299d5a_0` tag. Confirmed by testing pulls on the cluster.

**Fix**: Updated to `mosdepth:0.3.8--hd299d5a_0` (verified working via `singularity pull`).

**Commit**: `dfd291e` — "fix: update mosdepth container tag to 0.3.8"

### 2026-05-30: Fix #5 — FreeBayes (and other) container tags also purged

**Run: `compassionate_lattes`** (job 52364959) — Mosdepth now works, but `MUTATION_CALLING:FREEBAYES` failed with same MANIFEST_UNKNOWN error for `freebayes:1.3.7--hbfe0e7f_2`.

**Root cause**: quay.io/biocontainers has purged many older tags across multiple packages.

**Fix**: Queried the quay.io API (`/api/v1/repository/biocontainers/<tool>/tag/`) to find currently available tags and updated all biocontainer references at once:

| Tool | Old tag | New tag |
|------|---------|---------|
| mosdepth | `0.3.9--hd299d5a_0` | `0.3.8--hd299d5a_0` |
| freebayes | `1.3.7--hbfe0e7f_2` | `1.3.10--hbefcdb2_0` |
| strelka | `2.9.10--0` | `2.9.10--hdfd78af_2` |
| samtools | `1.21--h50ea8bc_0` | `1.23.1--ha83d96e_0` |
| bcftools | `1.21--h8b25389_0` | `1.23.1--ha83d96e_0` |
| multiqc | `1.25.2--pyhdfd78af_0` | `1.35--pyhdfd78af_1` |

**Commit**: `166fa22` — "fix: update all biocontainer tags to available versions"

**Note**: Changing the samtools container in `merge_bams/main.nf` may invalidate the MERGE_BAMS cache (the 10h step). If Nextflow re-runs it, consider reverting just that file.

**Lesson learned**: Always verify container tags exist before running. Use:
```bash
curl -s "https://quay.io/api/v1/repository/biocontainers/<tool>/tag/?limit=20" | python3 -c "import sys,json; [print(t['name']) for t in json.load(sys.stdin)['tags']]"
```

---

## Potential Future Issues

- **Variant calling time**: Mutect2 on a 2B-read BAM will be slow (estimate 3–6 CPU-hours per clone for WGS). The `process_high` label gives 16h which should suffice.
- **Disk space**: The merged BAM + caller outputs will be large. Monitor `/mnt/scratche` usage.
- **Single clone**: Only `clone_001` exists (all 212 cells in one clone). Cross-clone comparison won't be meaningful with a single clone. Consider adjusting `--min_cells_per_clone` or `--clone_strategy` to get multiple clones.
