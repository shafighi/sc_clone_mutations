#!/bin/bash
#SBATCH --job-name=sc_clone_mutations
#SBATCH --partition=epyc
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=48:00:00
#SBATCH --output=sc_clone_mutations_%j.log

# ─── Step 1: Load required modules ───────────────────────────────────────────
# Adjust module names to match your cluster (check with: module avail)
module load apptainer 2>/dev/null || module load singularity 2>/dev/null || true
module load java 2>/dev/null || true

# ─── Step 2: Set up Nextflow (install if not present) ─────────────────────────
export NXF_HOME=/mnt/scratche/slow/fmlab/darvis01/.nextflow
mkdir -p $NXF_HOME

if ! command -v nextflow &> /dev/null; then
    if [ ! -f "$NXF_HOME/nextflow" ]; then
        echo "Installing Nextflow..."
        cd $NXF_HOME
        curl -s https://get.nextflow.io | bash
        cd -
    fi
    export PATH=$NXF_HOME:$PATH
fi

# ─── Step 3: Set up Singularity cache ────────────────────────────────────────
export NXF_SINGULARITY_CACHEDIR=/mnt/scratche/slow/fmlab/darvis01/singularity_cache
mkdir -p $NXF_SINGULARITY_CACHEDIR

# ─── Step 4: Set work directory on scratch ────────────────────────────────────
WORK_DIR=/mnt/scratche/slow/fmlab/darvis01/nf_work
mkdir -p $WORK_DIR

# ─── Step 5: Clone/update the pipeline ───────────────────────────────────────
PIPELINE_DIR=/mnt/scratche/slow/fmlab/darvis01/sc_clone_mutations
if [ ! -d "$PIPELINE_DIR" ]; then
    git clone https://github.com/shafighi/sc_clone_mutations.git $PIPELINE_DIR
else
    cd $PIPELINE_DIR && git pull && cd -
fi

# ─── Step 6: Run the pipeline ────────────────────────────────────────────────
nextflow run $PIPELINE_DIR/main.nf \
    -profile singularity,slurm \
    --scunique_dir  /mnt/scratche/slow/fmlab/darvis01/scDNAseq-workflow/results/SLX-24491_100/ \
    --bam_dir       /mnt/scratche/slow/fmlab/darvis01/scDNAseq-workflow-main/data/aligned/SLX-24491/ \
    --fasta         /mnt/scratche/slow/fmlab/darvis01/snakemake-illumina-alignment/resources/homo_sapiens/GRCh37_g1kp2/fasta/hsa.GRCh37_g1kp2.fa \
    --tumor_only    true \
    --outdir        /mnt/scratche/slow/fmlab/darvis01/sc_clone_mutations_results \
    -work-dir       $WORK_DIR \
    -resume \
    -with-report    /mnt/scratche/slow/fmlab/darvis01/sc_clone_mutations_results/pipeline_info/report.html
