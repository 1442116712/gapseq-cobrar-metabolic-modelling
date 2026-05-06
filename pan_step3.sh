#!/bin/bash
#SBATCH --job-name=pan_step3
#SBATCH --nodes=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --time=00:15:00
#SBATCH --mail-type=END,FAIL

if [ -z "$SCRIPT_DIR" ]; then
    SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
fi
source "$SCRIPT_DIR/pipeline.config"

module load mamba/1.3.1
source activate gapseq

WORKDIR="$PIPELINE_BASE/${SGB}_rep"
mkdir -p "$OUTPUT_TSV_DIR"

echo "========== Step 3: parsing reports for $SGB =========="

# 3a. Aggregate per-MAG step1 stdouts
python3 "$PARSE_STAGE1_PY" \
    "${SGB}_pan" \
    "$WORKDIR/models/logs" \
    "$OUTPUT_TSV_DIR/all_stage1_summary.tsv"

# 3b. Parse step2 stdout
STEP2_LOG=$(ls -t "$WORKDIR/pan_model/logs/stdout_step2_"*.txt 2>/dev/null | head -1)
if [ -s "$STEP2_LOG" ]; then
    python3 "$PARSE_REPORT_PY" "${SGB}_pan" "$STEP2_LOG" "$OUTPUT_TSV_DIR/all_reconstruction.tsv"
else
    echo "WARN: no step2 stdout found in $WORKDIR/pan_model/logs/"
fi

echo ">>> Step 3 done for $SGB."
