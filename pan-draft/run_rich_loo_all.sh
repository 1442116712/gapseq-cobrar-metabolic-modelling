#!/bin/bash
# Dispatch rich-medium LOO across all reconstructed models.
#
# Usage:
#   bash run_rich_loo_all.sh pan   <SGB_ID> [<SGB_ID> ...]   # pan models
#   bash run_rich_loo_all.sh single [<faa_dir>]              # all single-genome models
#
# Notes:
#   - Requires the minimal-medium LOO output to already exist
#     (i.e. pan_step4.sh or single_step.sh has been run for the model).
#   - Pan models with known pFBA degeneracy can be flagged via:
#       USE_FBA=TRUE bash run_rich_loo_all.sh pan MGYG000291777

set -e
MODE="${1:-}"; shift || true
[ -z "$MODE" ] && { echo "Usage: $0 pan <SGB_ID>...   |   $0 single [<faa_dir>]"; exit 1; }

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
COMMON="--partition=k2-hipri --mail-user=zwu18@qub.ac.uk --export=ALL,SCRIPT_DIR=$SCRIPT_DIR"

case "$MODE" in
  pan)
    [ $# -eq 0 ] && { echo "pan mode needs at least one SGB_ID"; exit 1; }
    source "$SCRIPT_DIR/pipeline.config"
    LOG="$OUTPUT_TSV_DIR/logs"; mkdir -p "$LOG"
    for SGB in "$@"; do
      WORK="$PIPELINE_BASE/${SGB}_rep/pan_model"
      [ -s "$WORK/panModel.RDS" ] || { echo "[$SGB] SKIP - missing $WORK/panModel.RDS"; continue; }
      JOB=$(sbatch $COMMON \
        --output="$LOG/stdout_rich_${SGB}_%j.txt" \
        --error="$LOG/stderr_rich_${SGB}_%j.txt" \
        --export=ALL,SGB=$SGB,SCRIPT_DIR=$SCRIPT_DIR,USE_FBA=${USE_FBA:-FALSE},RICH_LB=${RICH_LB:--10} \
        "$SCRIPT_DIR/pan_step5_rich.sh" | grep -oP '\d+$' | tail -1)
      echo "[$SGB] rich LOO submitted: $JOB"
    done
    ;;

  single)
    # locate the single-genome scripts dir (config + single_step_rich.sh).
    # Override via env var SINGLE_SCRIPT_DIR if your layout differs.
    if [ -z "${SINGLE_SCRIPT_DIR:-}" ]; then
        for cand in "$SCRIPT_DIR/../single_genome" "$SCRIPT_DIR/../single"; do
            [ -d "$cand" ] && SINGLE_SCRIPT_DIR="$(cd "$cand" && pwd)" && break
        done
    fi
    [ -z "${SINGLE_SCRIPT_DIR:-}" ] && { echo "ERROR: cannot find single-genome script dir; set SINGLE_SCRIPT_DIR"; exit 1; }
    [ -s "$SINGLE_SCRIPT_DIR/pipeline_single.config" ] || { echo "ERROR: no pipeline_single.config in $SINGLE_SCRIPT_DIR"; exit 1; }
    source "$SINGLE_SCRIPT_DIR/pipeline_single.config"
    WORKDIR="${1:-$SINGLE_BASE/faa}"
    [ -d "$WORKDIR" ] || { echo "Missing: $WORKDIR"; exit 1; }
    cd "$WORKDIR"
    N=$(ls *.faa 2>/dev/null | wc -l)
    [ "$N" -eq 0 ] && { echo "No .faa in $WORKDIR"; exit 1; }
    echo "Submitting rich-LOO array for $N single-genome models"
    LOG="$WORKDIR/logs"; mkdir -p "$LOG" "$OUTPUT_TSV_DIR/logs"
    sbatch $COMMON \
        --array=1-$N \
        --output="$LOG/stdout_rich_%A_%a.txt" \
        --error="$LOG/stderr_rich_%A_%a.txt" \
        --export=ALL,WORKDIR=$WORKDIR,SCRIPT_DIR=$SINGLE_SCRIPT_DIR,USE_FBA=${USE_FBA:-FALSE},RICH_LB=${RICH_LB:--10} \
        "$SINGLE_SCRIPT_DIR/single_step_rich.sh"
    ;;

  *)
    echo "Unknown MODE: $MODE  (expected 'pan' or 'single')"
    exit 1
    ;;
esac
