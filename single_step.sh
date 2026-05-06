#!/bin/bash
#SBATCH --job-name=single
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=40G
#SBATCH --time=03:00:00
#SBATCH --mail-type=END,FAIL

if [ -z "$SCRIPT_DIR" ]; then
    SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
fi
source "$SCRIPT_DIR/pipeline_single.config"

module load mamba/1.3.1
source activate gapseq

cd "$WORKDIR"
FILES=(*.faa)
CURRENT="${FILES[$SLURM_ARRAY_TASK_ID - 1]}"
ID=$(basename "$CURRENT" .faa)

echo "========== Single full pipeline: $ID =========="

$GAPSEQ find -p all -t Bacteria -M prot -K 16 "$CURRENT"
$GAPSEQ find-transport -K 16 "$CURRENT"
$GAPSEQ draft \
    -r "${ID}-all-Reactions.tbl" \
    -t "${ID}-Transporter.tbl" \
    -p "${ID}-all-Pathways.tbl" \
    -b auto

$GAPSEQ medium -m "${ID}-draft.RDS" -p "${ID}-all-Pathways.tbl" -f .
$GAPSEQ fill -m "${ID}-draft.RDS" -n "${ID}-medium.csv" -b 100

LOG_FILE="$WORKDIR/logs/stdout_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.txt"
mkdir -p "$OUTPUT_TSV_DIR"
if [ -s "$LOG_FILE" ]; then
    python3 "$PARSE_REPORT_PY" "$ID" "$LOG_FILE" "$OUTPUT_TSV_DIR/all_reconstruction.tsv"
fi

export MODEL_FILE="$WORKDIR/${ID}.RDS"
export OUTPUT_PREFIX="$ID"
export ORGANISM_ID="$ID"
export NUTRIENTS_TSV
export OUTPUT_DIR="$OUTPUT_TSV_DIR"

if [ -s "$MODEL_FILE" ]; then
    Rscript "$LEAVE_ONE_OUT_R"
else
    echo "WARN: no $MODEL_FILE produced, skipping phenotyping"
fi

MERGED="$OUTPUT_TSV_DIR/merged"
mkdir -p "$MERGED"
cd "$OUTPUT_TSV_DIR"
for kind in auxotrophies substrates_full substrates_top15 summary; do
    files=$(ls *_${kind}.tsv 2>/dev/null | grep -v "^all_" | grep -v "^merged")
    [ -n "$files" ] && awk 'FNR==1 && NR!=1 {next} {print}' $files > "$MERGED/all_${kind}.tsv"
done

echo ">>> $ID done."
