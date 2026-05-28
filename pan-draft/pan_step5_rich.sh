#!/bin/bash
#SBATCH --job-name=pan_step5_rich
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=20G
#SBATCH --time=01:00:00
#SBATCH --mail-type=END,FAIL
#
# Rich-medium LOO for a single pan model.
# Depends on pan_step4.sh having produced
#   $OUTPUT_TSV_DIR/${SGB}_pan_auxotrophies.tsv
# (the minimal-medium LOO output), which leave_one_out_rich.R uses
# for the Tier A vs Tier C cross-reference.

if [ -z "$SCRIPT_DIR" ]; then
    SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
fi
source "$SCRIPT_DIR/pipeline.config"

LEAVE_ONE_OUT_RICH_R="${LEAVE_ONE_OUT_RICH_R:-$PIPELINE_BASE/leave_one_out_rich.R}"

module load mamba/1.3.1
source activate gapseq

PAN="$PIPELINE_BASE/${SGB}_rep/pan_model"
mkdir -p "$OUTPUT_TSV_DIR"

[ -s "$PAN/panModel.RDS" ] || { echo "MISSING: $PAN/panModel.RDS"; exit 1; }

export MODEL_FILE="$PAN/panModel.RDS"
export OUTPUT_PREFIX="${SGB}_pan"
export ORGANISM_ID="${SGB}_pan"
export NUTRIENTS_TSV
export OUTPUT_DIR="$OUTPUT_TSV_DIR"
# One pan model needs FBA-mode (same one as pan_step4); set via env if needed:
#   USE_FBA=TRUE for MGYG000291777
export USE_FBA="${USE_FBA:-FALSE}"
export RICH_LB="${RICH_LB:--10}"

echo "========== Step 5 (rich-medium LOO) for $SGB =========="
Rscript "$LEAVE_ONE_OUT_RICH_R"

# Merge per-SGB rich LOO TSVs into merged/ subdir
MERGED="$OUTPUT_TSV_DIR/merged"
mkdir -p "$MERGED"
cd "$OUTPUT_TSV_DIR"
awk 'FNR==1 && NR!=1 {next} {print}' *_pan_rich_LOO.tsv         > "$MERGED/all_rich_LOO.tsv"         2>/dev/null
awk 'FNR==1 && NR!=1 {next} {print}' *_pan_rich_LOO_summary.tsv > "$MERGED/all_rich_LOO_summary.tsv" 2>/dev/null

echo ">>> Step 5 done for $SGB. Cumulative rich-LOO TSVs in $MERGED/"
