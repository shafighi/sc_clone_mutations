---
name: nfcore-pipeline
description: Build, debug, and operate nf-core-style Nextflow (DSL2) pipelines on HPC/SLURM with containers. Use when creating modules/subworkflows, wiring resource labels, diagnosing failed runs, or fixing -resume cache and clone/grouping logic.
---

# nf-core-style Nextflow pipelines

Practical guidance distilled from building and debugging a clone-aware
genomics pipeline (Nextflow DSL2 + Singularity + SLURM). Applies to any
nf-core-style pipeline: `main.nf` → subworkflows → `modules/local/*` with
process resource labels in `conf/base.config`.

## Structure conventions
- `main.nf` orchestrates `subworkflows/*.nf`, which `include` `modules/local/<name>/main.nf`.
- Each module declares a `label` (`process_low|medium|high|...`) defined in `conf/base.config`.
- Cluster caps live in `conf/resources.config` (`max_cpus`, `max_memory`, `max_time`)
  applied via a `check_max(...)` helper.
- Keep per-process containers pinned (`quay.io/biocontainers/<tool>:<version>`).

## Resource labels and auto-scaling retries
- `base.config` typically multiplies by `task.attempt`
  (e.g. `memory = { check_max( 36.GB * task.attempt, 'memory') }`).
- For steps that can OOM on large inputs, add retry so attempts scale up:
  ```groovy
  process BIG_STEP {
      label 'process_high'
      errorStrategy 'retry'
      maxRetries 2
  }
  ```
- Exit status `140` (or "terminated by the external system" with empty
  exit) almost always means SLURM killed the job for **memory or time** —
  not a code bug. Bump the label / add retry, or reduce input size.
- Raising memory has a ceiling: it's capped by `max_memory`. If a step is
  killed even at the cap, the real fix is smaller inputs (see below), not
  more RAM.

## The -resume cache: critical gotchas
- `-resume` keys each task on inputs + the **process script block** +
  **container reference** — NOT on the contents of helper scripts in `bin/`.
- Therefore editing a `bin/*.py` (or any externally-called script) does
  **not** invalidate the cache. `-resume` will reuse the stale result.
  To force a changed `bin/` step to re-run: run **without** `-resume`
  (or change the process script/inputs/container so the hash changes).
- The resume cache lives in `.nextflow/cache/` in the **launch directory**,
  keyed by session. Changing `-work-dir` alone does NOT bypass it — old
  work files are still referenced by hash.
- Changing a container tag (e.g. `samtools:1.21` → `1.23.1`) **does**
  invalidate that task and everything downstream, forcing expensive
  re-runs. Pin containers and avoid bumping tags mid-project.
- A published output file in `--outdir` is not proof a step re-ran; check
  the file timestamp and the run log, not just its contents.

## Debugging a failed run (fast path)
1. `nextflow log` — list runs and STATUS (`OK`/`ERR`).
2. Inspect the error from the relevant run:
   `grep -A10 "ERROR\|Caused by" .nextflow.log | tail -40`
3. Identify the failing process and its work dir from the
   "Work dir:" line, then read the real stderr:
   `cat <workdir>/.command.err` and `.command.log`.
4. Reproduce in place: `cd <workdir> && bash .command.run`.
5. For task status fields: `nextflow log <run> -f name,status,exit,workdir`
   (note: `error` is not a valid field name).

## Input size is often the real root cause
- When a pipeline fans inputs into per-group jobs (e.g. merge cells →
  pseudobulk per clone), an upstream **grouping bug** can dump everything
  into one giant group, making every downstream step (merge, dedup,
  caller) intractable. Symptoms: a single huge work dir, multi-hour
  merges, OOM that more RAM can't fix.
- Diagnose by inspecting the grouping output (e.g. a `*_summary.csv`)
  BEFORE blaming resources. One group with all items = grouping bug.
- Fix the grouping logic; don't keep raising memory/time.

## Tree / clustering grouping logic (domain)
- When cutting a tree into groups by branch length, **never select the
  root** (it has no incoming edge) — selecting it collapses all leaves
  into one group.
- Use a strict threshold (`>`), and when assigning nested selected clades,
  process **deepest (smallest) clades first** so ancestors don't claim
  cells before their descendants.
- A sensible automatic threshold: `mean + 1*std` of internal edge lengths,
  so only the longest branches (clear boundaries) are cut.
- Verify grouping logic with a tiny synthetic tree (long inter-group
  branches, short intra-group branches) and assert `n_groups > 1`.

## SLURM submission pattern
- The submit job runs only the Nextflow head process; keep it small
  (`--mem=8G`, 1–2 cpus). Heavy work is submitted as separate SLURM jobs
  by the `slurm` executor with the labels' resources.
- Give the head job ample `--time` (long pipelines run many hours/days).
- Set a writable Singularity cache and work dir on scratch:
  `export NXF_SINGULARITY_CACHEDIR=/scratch/$USER/singularity_cache`,
  `-work-dir /scratch/$USER/nf_work`.
- `tput: libtinfow.so.6` errors in logs are harmless terminal-capability
  warnings, not pipeline failures.

## Domain caveat: dedup on merged single-cell pseudobulk
- Marking duplicates on a BAM merged from many single cells is usually
  wrong: identical start/end reads from different cells are independent
  molecules, not PCR duplicates. Prefer skipping it (`--mark_duplicates
  false`) for pseudobulk; somatic callers handle duplicates themselves.

## Keep a journey/plans log
- Record each failure, its root cause, the decision, and the commit hash.
  It prevents re-patching the same symptom and makes handoffs clean.
