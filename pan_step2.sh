#!/bin/bash
#SBATCH --job-name=pan_step2
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
PAN="$WORKDIR/pan_model"
mkdir -p "$PAN"

echo "========== Step 2: pan-Draft + medium + fill for $SGB =========="

# 1. Build pan-draft
$GAPSEQ pan -m "$MODELS/" -w "$MODELS/" -f "$PAN/"

# 2. Predict gap-fill medium
$GAPSEQ medium -m "$PAN/panModel-draft.RDS" -p "$PAN/panModel-tmp-Pathways.tbl" -f "$PAN/"

# 3. Gap-fill (fill ignores -f, run from PAN so output lands here)
cd "$PAN"
$GAPSEQ fill -m "panModel-draft.RDS" -n "panModel-medium.csv" -b 100

echo ">>> Pan-GEM for $SGB complete."
