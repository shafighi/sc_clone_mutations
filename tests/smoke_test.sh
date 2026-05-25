#!/usr/bin/env bash
# tests/smoke_test.sh
#
# Smoke test: runs Python scripts on example data without Nextflow/containers.
# Useful for rapid development and CI checks.
#
# Requirements:
#   pip install -r containers/requirements_python.txt
#   dendropy, pandas, scipy, scikit-learn, pysam, matplotlib, seaborn
#
# Usage:
#   bash tests/smoke_test.sh

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
BIN="${ROOT}/bin"
DATA="${ROOT}/examples/data"
TMPDIR="$(mktemp -d /tmp/scclone_smoke_XXXXXX)"

echo "=== sc_clone_mutations smoke test ==="
echo "Root:    ${ROOT}"
echo "Scratch: ${TMPDIR}"
echo ""

# Add bin to PATH for inter-script imports
export PYTHONPATH="${BIN}:${PYTHONPATH:-}"

# ── Step 1: Validate inputs ───────────────────────────────────────────────────
echo "[1/6] validate_inputs.py"
python "${BIN}/validate_inputs.py" \
    --bam_manifest      "${DATA}/bam_manifest.csv" \
    --cell_metadata     "${DATA}/cell_metadata.csv" \
    --medicc2_tree      "${DATA}/medicc2_tree.new" \
    --scunique_events   "${DATA}/scunique_events.tsv" \
    --out_bam_manifest  "${TMPDIR}/bam_manifest_validated.csv" \
    --out_cell_metadata "${TMPDIR}/cell_metadata_validated.csv" \
    --out_report        "${TMPDIR}/validation_report.json" \
    --min_mapped_reads  0 \
    --skip_bam_check    # example data uses placeholder paths
echo "    OK"

# ── Step 2: Parse MEDICC2 tree ────────────────────────────────────────────────
echo "[2/6] parse_medicc2_tree.py"
python "${BIN}/parse_medicc2_tree.py" \
    --tree_file      "${DATA}/medicc2_tree.new" \
    --out_pkl        "${TMPDIR}/tree_data.pkl" \
    --out_nodes      "${TMPDIR}/node_table.csv" \
    --out_edges      "${TMPDIR}/edge_table.csv" \
    --out_distances  "${TMPDIR}/pairwise_dist.csv"
echo "    OK"

# ── Step 3: Assign clones ─────────────────────────────────────────────────────
echo "[3/6] assign_clones.py"
python "${BIN}/assign_clones.py" \
    --tree_pkl             "${TMPDIR}/tree_data.pkl" \
    --node_table           "${TMPDIR}/node_table.csv" \
    --scunique_events      "${DATA}/scunique_events.tsv" \
    --cell_metadata        "${TMPDIR}/cell_metadata_validated.csv" \
    --bam_manifest         "${TMPDIR}/bam_manifest_validated.csv" \
    --strategy             internal_node \
    --min_cells_per_clone  2 \
    --min_branch_length    0.0 \
    --event_similarity_thr 0.5 \
    --small_clone_action   merge \
    --out_assignments      "${TMPDIR}/cell_clone_assignments.csv" \
    --out_summary          "${TMPDIR}/clone_summary.csv" \
    --out_events           "${TMPDIR}/clone_events.csv" \
    --out_tree             "${TMPDIR}/clone_tree_annotated.new"
echo "    OK — assignments:"
awk -F',' 'NR>1 {print "      " $1, "→", $2}' "${TMPDIR}/cell_clone_assignments.csv" | head -15

# ── Step 4: Plot clone tree ───────────────────────────────────────────────────
echo "[4/6] plot_clone_tree.py"
python "${BIN}/plot_clone_tree.py" \
    --tree_pkl    "${TMPDIR}/tree_data.pkl" \
    --assignments "${TMPDIR}/cell_clone_assignments.csv" \
    --out_pdf     "${TMPDIR}/clone_tree.pdf" \
    --out_png     "${TMPDIR}/clone_tree.png"
echo "    OK"

# ── Step 5: Generate report (with dummy inputs) ───────────────────────────────
echo "[5/6] generate_report.py"
# Create dummy consensus table
echo "chrom,pos,ref,alt" > "${TMPDIR}/consensus_table.csv"
echo "chr1,100000,A,G"  >> "${TMPDIR}/consensus_table.csv"
# Create dummy variant matrix
echo "chrom,pos,ref,alt" > "${TMPDIR}/variant_matrix.csv"

python "${BIN}/generate_report.py" \
    --validation_report  "${TMPDIR}/validation_report.json" \
    --clone_summary      "${TMPDIR}/clone_summary.csv" \
    --consensus_table    "${TMPDIR}/consensus_table.csv" \
    --cross_clone_matrix "${TMPDIR}/variant_matrix.csv" \
    --pipeline_version   "SMOKE_TEST" \
    --out_html           "${TMPDIR}/report.html" \
    --out_md             "${TMPDIR}/report.md"
echo "    OK"

# ── Step 6: Test all strategies ───────────────────────────────────────────────
echo "[6/6] Testing all clone strategies"
for strategy in distance event_profile hybrid; do
    extra_args=""
    if [ "$strategy" = "distance" ]; then
        extra_args="--distance_threshold 0.3"
    fi
    python "${BIN}/assign_clones.py" \
        --tree_pkl             "${TMPDIR}/tree_data.pkl" \
        --node_table           "${TMPDIR}/node_table.csv" \
        --scunique_events      "${DATA}/scunique_events.tsv" \
        --cell_metadata        "${TMPDIR}/cell_metadata_validated.csv" \
        --bam_manifest         "${TMPDIR}/bam_manifest_validated.csv" \
        --strategy             "${strategy}" \
        --min_cells_per_clone  2 \
        --min_branch_length    0.0 \
        --event_similarity_thr 0.5 \
        --small_clone_action   merge \
        ${extra_args} \
        --out_assignments      "${TMPDIR}/assignments_${strategy}.csv" \
        --out_summary          "${TMPDIR}/summary_${strategy}.csv" \
        --out_events           "${TMPDIR}/events_${strategy}.csv" \
        --out_tree             "${TMPDIR}/tree_${strategy}.new" \
        > /dev/null 2>&1
    n_clones=$(awk -F',' 'NR>1 {print $2}' "${TMPDIR}/assignments_${strategy}.csv" | sort -u | wc -l | tr -d ' ')
    echo "    ${strategy}: ${n_clones} clones"
done

echo ""
echo "=== All smoke tests PASSED ==="
echo "Outputs written to ${TMPDIR}"
ls -lh "${TMPDIR}"
