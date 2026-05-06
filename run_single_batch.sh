#!/bin/bash
# Run single-genome pipeline on all .faa files in a directory.
# Usage: bash run_single_batch.sh [<faa_dir>]
# Default: $SINGLE_BASE/faa

set -e
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/pipeline_single.config"

WORKDIR=${1:-"$SINGLE_BASE/faa"}
[ -d "$WORKDIR" ] || { echo "Missing: $WORKDIR"; exit 1; }

cd "$WORKDIR"
N=$(ls *.faa 2>/dev/null | wc -l)
[ "$N" -eq 0 ] && { echo "No .faa in $WORKDIR"; exit 1; }
echo "Found $N .faa files in $WORKDIR"

LOG_DIR="$WORKDIR/logs"
mkdir -p "$LOG_DIR" "$OUTPUT_TSV_DIR/logs"

sbatch \
    --partition=$SLURM_PARTITION \
    --mail-user=$SLURM_MAIL \
    --array=1-$N \
    --output="$LOG_DIR/stdout_%A_%a.txt" \
    --error="$LOG_DIR/stderr_%A_%a.txt" \
    --export=ALL,WORKDIR=$WORKDIR,SCRIPT_DIR=$SCRIPT_DIR \
    "$SCRIPT_DIR/single_step.sh"
