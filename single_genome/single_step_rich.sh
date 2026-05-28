#!/bin/bash
#SBATCH --job-name=single_rich
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=20G
#SBATCH --time=01:00:00
#SBATCH --mail-type=END,FAIL
#
# Rich-medium LOO for a single-genome model (reference, ARG-MAG or in-vitro).
# Depends on single_step.sh having produced
#   $OUTPUT_TSV_DIR/${ID}_auxotrophies.tsv
# (minimal-medium LOO output) for the Tier A vs Tier C cross-reference.

if [ -z "$SCRIPT_DIR" ]; then
    SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
fi
source "$SCRIPT_DIR/pipeline_single.config"

LEAVE_ONE_OUT_RICH_R="${LEAVE_ONE_OUT_RICH_R:-$PIPELINE_BASE/leave_one_out_rich.R}"

module load mamba/1.3.1
source activate gapseq

cd "$WORKDIR"
FILES=(*.faa)
CURRENT="${FILES[$SLURM_ARRAY_TASK_ID - 1]}"
ID=$(basename "$CURRENT" .faa)

MODEL_FILE="$WORKDIR/${ID}.RDS"
[ -s "$MODEL_FILE" ] || { echo "MISSING: $MODEL_FILE"; exit 1; }

export MODEL_FILE
export OUTPUT_PREFIX="$ID"
export ORGANISM_ID="$ID"
export NUTRIENTS_TSV
export OUTPUT_DIR="$OUTPUT_TSV_DIR"
export USE_FBA="${USE_FBA:-FALSE}"
export RICH_LB="${RICH_LB:--10}"

echo "========== Rich-medium LOO: $ID =========="
Rscript "$LEAVE_ONE_OUT_RICH_R"

MERGED="$OUTPUT_TSV_DIR/merged"
mkdir -p "$MERGED"
cd "$OUTPUT_TSV_DIR"
awk 'FNR==1 && NR!=1 {next} {print}' *_rich_LOO.tsv         | grep -v "^all_" > "$MERGED/all_rich_LOO.tsv"         2>/dev/null
awk 'FNR==1 && NR!=1 {next} {print}' *_rich_LOO_summary.tsv | grep -v "^all_" > "$MERGED/all_rich_LOO_summary.tsv" 2>/dev/null

echo ">>> Rich-LOO done for $ID."
