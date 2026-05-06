#!/bin/bash
# Usage: bash run_pan_pipeline.sh <SGB_ID> [<SGB_ID> ...]

set -e
[ $# -eq 0 ] && { echo "Usage: $0 <SGB_ID> [<SGB_ID> ...]"; exit 1; }

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/pipeline.config"

extract_jobid() { grep -oP '\d+$' | tail -1; }

COMMON_ARGS="--partition=$SLURM_PARTITION --mail-user=$SLURM_MAIL"

for SGB in "$@"; do
    WORKDIR="$PIPELINE_BASE/${SGB}_rep"
    [ -d "$WORKDIR" ] || { echo "[$SGB] SKIP - missing $WORKDIR"; continue; }

    N=$(ls "$WORKDIR"/*.faa 2>/dev/null | wc -l)
    [ "$N" -eq 0 ] && { echo "[$SGB] SKIP - no .faa"; continue; }

    LOG1="$WORKDIR/models/logs"
    LOG2="$WORKDIR/pan_model/logs"
    LOG34="$OUTPUT_TSV_DIR/logs"
    mkdir -p "$LOG1" "$LOG2" "$LOG34"

    echo "========== [$SGB] N=$N MAGs =========="

    JOB1=$(sbatch $COMMON_ARGS \
        --array=1-$N \
        --output="$LOG1/stdout_step1_%A_%a.txt" \
        --error="$LOG1/stderr_step1_%A_%a.txt" \
        --export=ALL,SGB=$SGB,SCRIPT_DIR=$SCRIPT_DIR \
        "$SCRIPT_DIR/pan_step1.sh" | extract_jobid)
    echo "[$SGB] step1=$JOB1 (array 1-$N)"

    JOB2=$(sbatch $COMMON_ARGS \
        --dependency=afterok:$JOB1 \
        --output="$LOG2/stdout_step2_%j.txt" \
        --error="$LOG2/stderr_step2_%j.txt" \
        --export=ALL,SGB=$SGB,SCRIPT_DIR=$SCRIPT_DIR \
        "$SCRIPT_DIR/pan_step2.sh" | extract_jobid)
    echo "[$SGB] step2=$JOB2"

    JOB3=$(sbatch $COMMON_ARGS \
        --dependency=afterok:$JOB2 \
        --output="$LOG34/stdout_step3_${SGB}_%j.txt" \
        --error="$LOG34/stderr_step3_${SGB}_%j.txt" \
        --export=ALL,SGB=$SGB,SCRIPT_DIR=$SCRIPT_DIR \
        "$SCRIPT_DIR/pan_step3.sh" | extract_jobid)
    echo "[$SGB] step3=$JOB3"

    JOB4=$(sbatch $COMMON_ARGS \
        --dependency=afterok:$JOB3 \
        --output="$LOG34/stdout_step4_${SGB}_%j.txt" \
        --error="$LOG34/stderr_step4_${SGB}_%j.txt" \
        --export=ALL,SGB=$SGB,SCRIPT_DIR=$SCRIPT_DIR \
        "$SCRIPT_DIR/pan_step4.sh" | extract_jobid)
    echo "[$SGB] step4=$JOB4"
done

echo ""
