# sc_clone_mutations

**Clone-aware somatic mutation analysis from single-cell DNA sequencing data**

Given scUnique copy-number outputs, a MEDICC2 phylogenetic tree, and per-cell BAM files, this pipeline:

1. Assigns cells to biologically meaningful clones using the tree topology and/or copy-number event profiles
2. Constructs clone-level pseudobulk BAMs by merging reads from cells in each clone
3. Calls somatic mutations in each pseudobulk with Mutect2, Strelka2, and FreeBayes
4. Compares mutations across clones (private vs shared) and across callers (concordance)
5. Builds a consensus callset and generates publication-ready QC reports

---

## Table of Contents

1. [Requirements](#1-requirements)
2. [Repository structure](#2-repository-structure)
3. [Inputs](#3-inputs)
4. [Quick start](#4-quick-start)
5. [Workflow stages in detail](#5-workflow-stages-in-detail)
6. [Clone assignment strategies](#6-clone-assignment-strategies)
7. [Container usage and HPC setup](#7-container-usage-and-hpc-setup)
8. [Configuration reference](#8-configuration-reference)
9. [Outputs](#9-outputs)
10. [Example commands](#10-example-commands)
11. [Extending the pipeline](#11-extending-the-pipeline)
12. [Assumptions and caveats](#12-assumptions-and-caveats)

---

## 1. Requirements

| Tool | Version | Used for |
|------|---------|---------|
| Nextflow | ≥23.04 | Workflow orchestration |
| Docker **or** Singularity/Apptainer | any recent | Container runtime |
| Java | ≥11 | Required by Nextflow |

All bioinformatics tools (GATK, Strelka2, samtools, etc.) run inside containers — **no manual cluster installation required**.

Install Nextflow:
```bash
curl -s https://get.nextflow.io | bash
# Or with conda: conda install -c bioconda nextflow
```

---

## 2. Repository structure

```
sc_clone_mutations/
├── README.md                   ← this file
├── nextflow.config             ← all parameters + profiles
├── main.nf                     ← top-level workflow
│
├── conf/
│   ├── base.config             ← resource labels (low/medium/high)
│   ├── resources.config        ← cluster-wide resource caps
│   ├── slurm.config            ← SLURM-specific settings
│   └── test.config             ← smoke-test parameters
│
├── subworkflows/
│   ├── validate_inputs.nf      ← Stage A: input validation
│   ├── clone_definition.nf     ← Stage B: clone assignment
│   ├── pseudobulk.nf           ← Stage C: BAM construction
│   ├── mutation_calling.nf     ← Stage D: variant calling
│   ├── variant_analysis.nf     ← Stage E: cross-clone comparison
│   └── reporting.nf            ← Stage F: MultiQC + HTML report
│
├── modules/local/
│   ├── validate_manifests/     ← validate_inputs.py wrapper
│   ├── parse_medicc2_tree/     ← parse_medicc2_tree.py wrapper
│   ├── assign_clones/          ← assign_clones.py wrapper
│   ├── plot_clone_tree/        ← plot_clone_tree.py wrapper
│   ├── filter_cells/           ← filter_cells.py wrapper
│   ├── merge_bams/             ← samtools merge + sort + index
│   ├── markduplicates/         ← picard MarkDuplicates
│   ├── mosdepth/               ← mosdepth coverage
│   ├── samtools_flagstat/      ← samtools flagstat
│   ├── mutect2/                ← GATK Mutect2
│   ├── filter_mutect2/         ← GATK FilterMutectCalls
│   ├── strelka2/               ← Strelka2 somatic/germline
│   ├── freebayes/              ← FreeBayes
│   ├── normalize_vcf/          ← bcftools norm
│   ├── compare_variants/       ← compare_variants.py wrapper
│   ├── build_consensus/        ← build_consensus.py wrapper
│   ├── multiqc/                ← MultiQC
│   └── custom_report/          ← generate_report.py wrapper
│
├── bin/
│   ├── parse_medicc2_tree.py   ← MEDICC2 tree parsing
│   ├── assign_clones.py        ← clone assignment (4 strategies)
│   ├── validate_inputs.py      ← input validation
│   ├── filter_cells.py         ← per-cell QC filter
│   ├── plot_clone_tree.py      ← tree visualisation
│   ├── compare_variants.py     ← cross-clone VCF comparison
│   ├── build_consensus.py      ← consensus callset
│   ├── generate_report.py      ← HTML/Markdown report
│   └── utils/
│       ├── tree_utils.py       ← dendropy helpers
│       └── vcf_utils.py        ← pysam/pandas VCF helpers
│
├── containers/
│   ├── Dockerfile.python       ← Python analysis container
│   ├── requirements_python.txt ← pinned Python dependencies
│   └── build.sh                ← container build script
│
├── assets/
│   └── multiqc_config.yml      ← MultiQC appearance settings
│
├── examples/
│   ├── data/                   ← synthetic 12-cell example dataset
│   └── params/                 ← example parameter YAML files
│
├── schemas/
│   └── bam_manifest.json       ← JSON schema for BAM manifest
│
└── tests/
    └── smoke_test.sh           ← Python-level smoke test (no Nextflow)
```

---

## 3. Inputs

### 3.1 BAM manifest (required)

CSV with one row per cell:

```
cell_id,bam_path,bai_path,sample_id,patient_id
CELL_001,/data/cells/CELL_001.bam,/data/cells/CELL_001.bam.bai,SAMPLE01,PATIENT01
```

- `cell_id` must match the leaf labels in the MEDICC2 Newick tree
- `bai_path` is optional (inferred as `{bam_path}.bai` if absent)
- `sample_id` groups cells from the same tumour sample
- `patient_id` is used to join cells with their matched normal

### 3.2 Cell metadata (required)

CSV with at minimum `cell_id`. Extra columns (ploidy, rpc, pass_qc, etc.) are carried through:

```
cell_id,sample_id,patient_id,ploidy,pass_qc
CELL_001,SAMPLE01,PATIENT01,2.1,TRUE
```

### 3.3 MEDICC2 tree (required)

Standard Newick format output from MEDICC2 (`*_final_tree.new`).
Internal node labels are optional but helpful.

### 3.4 scUnique unique events (required)

TSV with at minimum `cell_id` plus coordinate columns. Expected:

```
cell_id	chr	start	end	event_type	copy_number
CELL_001	chr5	20000000	80000000	gain	4
```

Column names can vary — see `assign_clones.py` for the auto-detection logic.
**TODO**: verify that your scUnique output column names match what the script expects.

### 3.5 Matched normal manifest (optional)

CSV with one row per normal sample:

```
sample_id,patient_id,bam_path,bai_path
NORMAL01,PATIENT01,/data/normals/NORMAL01.bam,
```

If absent, the pipeline runs in tumor-only mode (`--tumor_only true`).

### 3.6 Reference genome resources

| Parameter | File | Required |
|-----------|------|---------|
| `--fasta` | Reference FASTA | Yes |
| `--fai` | FASTA index | Auto-inferred |
| `--dict` | Sequence dictionary | Auto-inferred |
| `--germline_resource` | gnomAD or similar VCF.gz + TBI | Strongly recommended for Mutect2 |
| `--panel_of_normals` | PoN VCF.gz + TBI | Recommended for Mutect2 |
| `--intervals` | Target BED / interval list | Optional (whole genome if absent) |

GATK resource bundles for hg38 are available from the
[GATK Resource Bundle](https://gatk.broadinstitute.org/hc/en-us/articles/360035890811).

---

## 4. Quick start

### Smoke test (Python only, no containers)

```bash
# Install Python dependencies
pip install -r containers/requirements_python.txt

# Run smoke test on example data
bash tests/smoke_test.sh
```

### Local Docker run

```bash
nextflow run main.nf \
    -profile docker \
    --bam_manifest   examples/data/bam_manifest.csv \
    --cell_metadata  examples/data/cell_metadata.csv \
    --medicc2_tree   examples/data/medicc2_tree.new \
    --scunique_events examples/data/scunique_events.tsv \
    --fasta          /ref/hg38/hg38.fa \
    --tumor_only     true \
    --outdir         results/
```

### HPC with Singularity + SLURM

```bash
nextflow run main.nf \
    -profile singularity,slurm \
    -params-file examples/params/full_run.yml \
    -resume
```

---

## 5. Workflow stages in detail

```
INPUT MANIFESTS
      │
      ▼
[A] VALIDATE_INPUTS
    validate_inputs.py
    - Check required columns, file existence, cell_id consistency
    - Output: validated CSVs, validation_report.json
      │
      ▼
[B] CLONE_DEFINITION
    parse_medicc2_tree.py  → pickled tree + node/edge tables
    assign_clones.py       → cell_clone_assignments.csv
    plot_clone_tree.py     → clone_tree.pdf / .png
      │
      ▼
[C] PSEUDOBULK
    filter_cells.py       → filter by QC thresholds
    samtools merge        → one BAM per clone
    picard MarkDuplicates → mark (not remove) duplicates
    mosdepth              → coverage QC
    samtools flagstat     → alignment QC
      │
      ▼
[D] MUTATION_CALLING                 (parallelised per clone)
    Mutect2 + FilterMutectCalls  →  *.mutect2.filtered.vcf.gz
    Strelka2                     →  *.strelka2.vcf.gz
    FreeBayes                    →  *.freebayes.vcf.gz
    bcftools norm (all)          →  *.norm.vcf.gz
      │
      ▼
[E] VARIANT_ANALYSIS
    compare_variants.py  → variant_matrix.csv, concordance, sharing
    build_consensus.py   → consensus_mutations.vcf.gz + table
    [optional] VEP       → annotated VCF
      │
      ▼
[F] REPORTING
    MultiQC              → multiqc_report.html
    generate_report.py   → sc_clone_mutations_report.html + .md
```

---

## 6. Clone assignment strategies

Set with `--clone_strategy`. Default: `internal_node`.

### `internal_node` (default)

Cuts the MEDICC2 tree at internal nodes whose incoming branch length
exceeds `--min_branch_length` and whose subtree contains ≥ `--min_cells_per_clone` cells.

**Biological rationale**: Long branches in a MEDICC2 tree reflect many
accumulated copy-number changes on a lineage — this defines a clonal
population that diverged substantially from its ancestor.

Best for: well-resolved trees with clear branching structure.

### `distance`

Hierarchical clustering (Ward linkage) on MEDICC2 pairwise tree distances.
Requires `--distance_threshold`.

Best for: trees without strong internal node support or when you want
distance-based grouping independent of branch topology.

### `event_profile`

Clusters cells by Jaccard similarity of their scUnique copy-number event
profiles. Tree-agnostic — purely phenotypic.

Best for: cases where the phylogenetic tree is uncertain or when you want
to validate tree-based clones against copy-number event similarity.

### `hybrid`

Two-stage approach:
1. Internal-node grouping (tree-based)
2. Each resulting group is checked for internal event homogeneity;
   groups with low within-clone Jaccard similarity are split further.

Best for: datasets where the tree is reliable for large-scale structure
but event profiles refine within-clade heterogeneity.

### Small clone handling (`--small_clone_action`)

| Value | Behaviour |
|-------|-----------|
| `merge` (default) | Merge small clones into the nearest/largest clone |
| `drop` | Remove cells from small clones entirely |
| `flag` | Keep small clones but tag them with a warning |

---

## 7. Container usage and HPC setup

### Docker (local)

```bash
# Build the Python analysis container
bash containers/build.sh

# Or pull existing public image (after pushing to GHCR):
docker pull ghcr.io/YOUR_ORG/scclone-python:1.0.0
```

Callers use public Biocontainers images — Nextflow pulls them automatically.

### Singularity / Apptainer (HPC, no root)

```bash
# Set cache directory (must be writable without root)
export NXF_SINGULARITY_CACHEDIR=/scratch/${USER}/singularity_cache
mkdir -p $NXF_SINGULARITY_CACHEDIR

# Pull the Python container as SIF
apptainer pull ${NXF_SINGULARITY_CACHEDIR}/scclone-python.sif \
    docker://ghcr.io/YOUR_ORG/scclone-python:1.0.0

# Run with singularity profile
nextflow run main.nf -profile singularity,slurm -params-file params.yml
```

Nextflow auto-converts Docker images to SIF format when running under
the `singularity` profile. Common HPC bind paths:

```groovy
// Add to nextflow.config if needed:
singularity {
    runOptions = '--bind /data:/data --bind /ref:/ref --bind /scratch:/scratch'
}
```

### SLURM submission script

```bash
#!/bin/bash
#SBATCH --job-name=sc_clone_mutations
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=4:00:00
#SBATCH --output=nf_submit_%j.log

module load apptainer/1.3.0
export NXF_SINGULARITY_CACHEDIR=/scratch/${USER}/singularity_cache

nextflow run /path/to/sc_clone_mutations/main.nf \
    -profile singularity,slurm \
    -params-file params.yml \
    -resume \
    -with-report results/pipeline_info/report.html
```

### Troubleshooting common HPC issues

| Issue | Fix |
|-------|-----|
| `Failed to pull singularity image` | Check `NXF_SINGULARITY_CACHEDIR` is writable; try `apptainer pull` manually |
| `Work directory not writable` | Set `-work-dir /scratch/${USER}/nf_work` |
| `SLURM out of memory` | Increase `--max_memory` or edit `conf/resources.config` |
| `Container: operation not permitted` | Add `--no-home` to `singularity.runOptions` |
| Temp files filling root | Set `NXF_TEMP=/scratch/${USER}/nf_tmp` |

---

## 8. Configuration reference

All parameters can be set via command line (`--param value`) or YAML
(`-params-file params.yml`). Key parameters:

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--clone_strategy` | `internal_node` | Clone assignment strategy |
| `--min_cells_per_clone` | 5 | Minimum cells to form a clone |
| `--min_branch_length` | 0.0 | Minimum tree branch length for clone cut |
| `--distance_threshold` | null | Required for `distance` strategy |
| `--event_similarity_thr` | 0.5 | Jaccard threshold for event-based strategies |
| `--small_clone_action` | `merge` | How to handle small clones |
| `--mark_duplicates` | `true` | Run picard MarkDuplicates on pseudobulks |
| `--tumor_only` | `false` | Tumor-only calling (no normal) |
| `--callers` | all three | Which callers to run |
| `--consensus_min_callers` | 2 | Callers required for consensus |
| `--annotate` | `false` | Run VEP annotation |

---

## 9. Outputs

```
results/
├── pipeline_info/              ← Nextflow run reports and traces
├── validation/
│   ├── bam_manifest_validated.csv
│   ├── cell_metadata_validated.csv
│   └── validation_report.json
│
├── clone_definition/
│   ├── cell_clone_assignments.csv   ← cell_id → clone_id mapping
│   ├── clone_summary.csv            ← per-clone statistics
│   ├── clone_events.csv             ← defining copy-number events per clone
│   ├── clone_tree_annotated.new     ← Newick with clone labels
│   └── figures/
│       ├── clone_tree.pdf
│       └── clone_tree.png
│
├── pseudobulk/
│   ├── bams/                        ← clone_NNN.markdup.bam + .bai
│   ├── cell_qc/                     ← filter_cells QC summary
│   └── qc/
│       ├── mosdepth/                ← coverage summaries
│       └── flagstat/                ← alignment stats
│
├── variant_calling/
│   ├── mutect2/
│   │   ├── raw/                     ← raw Mutect2 VCFs
│   │   └── filtered/                ← FilterMutectCalls output
│   ├── strelka2/                    ← Strelka2 VCFs
│   ├── freebayes/                   ← FreeBayes VCFs
│   └── */normalized/                ← bcftools norm output
│
├── variant_analysis/
│   ├── variant_matrix.csv           ← presence/absence × (clone, caller)
│   ├── caller_concordance.csv       ← per-clone concordance stats
│   ├── private_shared_table.csv     ← variant sharing classification
│   ├── per_clone_vcfs/              ← per-clone variant TSVs
│   ├── comparison_plots/            ← concordance, sharing, similarity plots
│   └── consensus/
│       ├── consensus_mutations.vcf.gz  ← consensus callset
│       ├── consensus_table.csv
│       └── consensus_summary.json
│
└── reports/
    ├── multiqc_report.html          ← MultiQC QC report
    └── sc_clone_mutations_report.html  ← custom summary report
```

---

## 10. Example commands

### Local Docker run

```bash
nextflow run main.nf \
    -profile docker \
    --bam_manifest      bam_manifest.csv \
    --cell_metadata     cell_metadata.csv \
    --medicc2_tree      sample_final_tree.new \
    --medicc2_events    sample_events_per_branch.tsv \
    --scunique_events   unique_events.tsv \
    --fasta             /ref/hg38/hg38.fa \
    --germline_resource /ref/hg38/af-only-gnomad.hg38.vcf.gz \
    --panel_of_normals  /ref/hg38/1000g_pon.hg38.vcf.gz \
    --normal_manifest   normal_manifest.csv \
    --outdir            results/
```

### Singularity local run

```bash
export NXF_SINGULARITY_CACHEDIR=/tmp/singularity_cache

nextflow run main.nf \
    -profile singularity \
    --bam_manifest      bam_manifest.csv \
    --cell_metadata     cell_metadata.csv \
    --medicc2_tree      sample_final_tree.new \
    --scunique_events   unique_events.tsv \
    --fasta             /ref/hg38/hg38.fa \
    --tumor_only        true \
    --outdir            results/
```

### SLURM cluster execution

```bash
export NXF_SINGULARITY_CACHEDIR=/scratch/${USER}/singularity_cache

nextflow run /path/to/sc_clone_mutations/main.nf \
    -profile singularity,slurm \
    -params-file examples/params/full_run.yml \
    -work-dir    /scratch/${USER}/nf_work \
    -resume
```

### Smoke test (no Nextflow, no containers)

```bash
pip install -r containers/requirements_python.txt
bash tests/smoke_test.sh
```

### Changing clone strategy at runtime

```bash
# Use event-profile clustering instead of tree topology
nextflow run main.nf -profile docker \
    --clone_strategy      event_profile \
    --event_similarity_thr 0.6 \
    --min_cells_per_clone 3 \
    ...
```

---

## 11. Extending the pipeline

### Adding a new caller

1. Create `modules/local/my_caller/main.nf` (follow the FreeBayes module as template)
2. Add the caller name to `params.callers` list in `nextflow.config`
3. Add an `if ('my_caller' in callers)` block in `subworkflows/mutation_calling.nf`
4. Add the caller container image

### Adding signature analysis

The consensus variant table and per-clone VCFs output by this pipeline are
directly usable as input to:
- **SigProfilerExtractor** (de novo signature discovery)
- **SigProfilerAssignment** (signature attribution)
- **mSigAct** (Bayesian signature assignment)

Clone-level VAF profiles from the variant matrix can also inform phylogeny-aware
SNV modelling (e.g. mapping mutations onto the MEDICC2 tree).

### Adding VEP annotation

```bash
nextflow run main.nf -profile docker \
    --annotate    true \
    --vep_cache   /ref/vep_cache \
    --vep_species homo_sapiens \
    --vep_genome_build GRCh38 \
    ...
```

The VEP annotation module (`modules/local/annotate_vep/main.nf`) requires a
pre-downloaded VEP cache directory.

---

## 12. Assumptions and caveats

### Cell identity

**Assumption**: MEDICC2 leaf labels exactly match `cell_id` in the BAM manifest.
If your MEDICC2 run uses a different cell naming convention (e.g. prefixed or
suffixed), normalise one side before running this pipeline.

### Pseudobulk and somatic calling

**Assumption**: Merging reads from cells assigned to the same clone creates a
pseudobulk BAM that behaves like a low-coverage bulk tumour BAM for variant
calling purposes.

**Caveat**: This is an approximation. Pseudobulk BAMs from few cells (< 10–20)
will have very low coverage and produce unreliable somatic calls regardless of
the caller. Clones with < 5 cells are flagged and should be interpreted with
great caution. The `--min_cells_per_clone` parameter directly controls this.

### Duplication marking

**Caveat**: scDNA-seq reads from different cells with the same start/end
coordinates are biologically distinct (they come from different cells, not PCR
duplicates of the same molecule). Picard MarkDuplicates does not know this and
will incorrectly flag cross-cell duplicates. Consider:

- Setting `--mark_duplicates false` (not marking at all)
- Or adding clone-aware read groups before merging (TODO in `merge_bams/main.nf`)

Current default: mark but do not remove, so callers can apply their own logic.

### Strelka2 tumor-only mode

**Caveat**: Strelka2 does not natively support true tumor-only somatic calling.
In `--tumor_only` mode, the pipeline runs Strelka2's germline workflow, which
calls variants relative to a germline prior — not true somatic calls.
For tumor-only analysis, **Mutect2 is the primary recommended caller**.

### FreeBayes somatic calling

**Caveat**: FreeBayes is not a dedicated somatic caller. It is included as an
additional supporting caller. Use Mutect2 and Strelka2 (paired mode) as primary
evidence. FreeBayes variants should be treated as supporting evidence only.

### Reference genome

Default is **hg38**. The genome build is not explicitly enforced — it is your
responsibility to ensure all inputs (BAMs, reference, VCF resources) use the
same build and chromosome naming convention (UCSC chr1 vs Ensembl 1).
`normalize_vcf.nf` runs `bcftools norm` which will catch some inconsistencies.

### Computational cost

Pseudobulk BAM construction and variant calling scale with the number of clones
and the depth of each pseudobulk. A run with 10 clones × 3 callers on a 30 Gb
genome requires roughly:

- 3–6 CPU-hours per clone for Mutect2
- 1–2 CPU-hours per clone for Strelka2
- 1–4 CPU-hours per clone for FreeBayes
- 10–40 GB RAM per Mutect2 job

Plan SLURM resource allocations accordingly via `conf/resources.config`.

### TODO markers in code

Search for `# TODO:` in the codebase for domain-specific tuning points, including:
- scUnique column name adaptation (varies by scUnique version)
- MEDICC2 branch event TSV column name adaptation
- Read-group tagging in pseudobulk BAMs
- Matched normal joining logic in `mutation_calling.nf`
- VEP annotation module (stub only)
