#!/bin/bash
#SBATCH --job-name=pan_step1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=40G
#SBATCH --time=3:00:00
#SBATCH --mail-type=END,FAIL

if [ -z "$SCRIPT_DIR" ]; then
    SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
fi
source "$SCRIPT_DIR/pipeline.config"

module load mamba/1.3.1
source activate gapseq

WORKDIR="$PIPELINE_BASE/${SGB}_rep"
MODELS="$WORKDIR/models"

cd "$WORKDIR"
FILES=(*.faa)
CURRENT="${FILES[$SLURM_ARRAY_TASK_ID - 1]}"
ID=$(basename "$CURRENT" .faa)

echo "========== Step 1: per-MAG draft for $CURRENT (task $SLURM_ARRAY_TASK_ID) =========="

$GAPSEQ find -p all -t Bacteria -M prot -K 16 -f "$MODELS" "$CURRENT"
$GAPSEQ find-transport -K 16 -f "$MODELS" "$CURRENT"
$GAPSEQ draft -f "$MODELS" \
    -r "$MODELS/${ID}-all-Reactions.tbl" \
    -t "$MODELS/${ID}-Transporter.tbl" \
    -p "$MODELS/${ID}-all-Pathways.tbl" \
    -b auto

echo ">>> $ID draft done."
