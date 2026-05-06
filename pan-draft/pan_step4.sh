#!/bin/bash
#SBATCH --job-name=pan_step4
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=20G
#SBATCH --time=02:00:00
#SBATCH --mail-type=END,FAIL

if [ -z "$SCRIPT_DIR" ]; then
    SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
fi
source "$SCRIPT_DIR/pipeline.config"

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

echo "========== Step 4: phenotyping for $SGB =========="
Rscript "$LEAVE_ONE_OUT_R"

# Auto-merge per-SGB TSVs into a separate merged/ subdir
MERGED="$OUTPUT_TSV_DIR/merged"
mkdir -p "$MERGED"
cd "$OUTPUT_TSV_DIR"

awk 'FNR==1 && NR!=1 {next} {print}' *_pan_auxotrophies.tsv     > "$MERGED/all_auxotrophies.tsv"     2>/dev/null
awk 'FNR==1 && NR!=1 {next} {print}' *_pan_substrates_full.tsv  > "$MERGED/all_substrates_full.tsv"  2>/dev/null
awk 'FNR==1 && NR!=1 {next} {print}' *_pan_substrates_top15.tsv > "$MERGED/all_substrates_top15.tsv" 2>/dev/null
awk 'FNR==1 && NR!=1 {next} {print}' *_pan_summary.tsv          > "$MERGED/all_summary.tsv"          2>/dev/null

echo ">>> Step 4 done for $SGB. Cumulative TSVs in $MERGED/"
