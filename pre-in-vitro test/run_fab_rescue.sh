#!/bin/bash
#SBATCH --job-name=fab_rescue
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=20G
#SBATCH --time=02:00:00
#SBATCH --partition=k2-medpri
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=zwu18@qub.ac.uk
#
# In-silico FAB rescue pre-experiment for DSM 11370.
#
# Before submitting:
#   1. scp the new FAB.csv to HPC, e.g.:
#        scp others/FAB.csv 40335635@kelvin.qub.ac.uk:/users/40335635/sharedscratch/EggNOG/manuscript_4/single_genome/FAB.csv
#   2. scp this script + fab_rescue_test.R to a working dir on HPC:
#        scp others/fab_rescue_test.R others/run_fab_rescue.sh 40335635@kelvin.qub.ac.uk:/users/40335635/sharedscratch/EggNOG/manuscript_4/single_genome/
#
# Submit on HPC:
#   cd /users/40335635/sharedscratch/EggNOG/manuscript_4/single_genome
#   sbatch run_fab_rescue.sh
#
# Override defaults if needed:
#   sbatch --export=ALL,USE_FBA=TRUE run_fab_rescue.sh

SINGLE_BASE="/users/40335635/sharedscratch/EggNOG/manuscript_4/single_genome"

export MODEL_FILE="${MODEL_FILE:-$SINGLE_BASE/faa/GCF_000426565.1.RDS}"
export FAB_CSV="${FAB_CSV:-$SINGLE_BASE/FAB.csv}"
export M2_CSV="${M2_CSV:-$SINGLE_BASE/M2_no_rumen.csv}"
export OUTPUT_DIR="${OUTPUT_DIR:-$SINGLE_BASE/output_tsv/fab_rescue}"
export USE_FBA="${USE_FBA:-FALSE}"
# nutrients.tsv lets the R script auto-discover cobalamin-family compounds
# by name (catches any cpd we didn't curate in the script)
export NUTRIENTS_TSV="${NUTRIENTS_TSV:-/users/40335635/sharedscratch/EggNOG/gapseq/gapseq-2.0.1/dat/nutrients.tsv}"

mkdir -p "$OUTPUT_DIR"

module load mamba/1.3.1
source activate gapseq

# SLURM copies the script to /var/spool/slurmd/... so BASH_SOURCE[0] no longer
# points to where you submitted from. Use SLURM_SUBMIT_DIR (set by sbatch to
# the directory you submitted from) instead.
SCRIPT_DIR="${SCRIPT_DIR:-${SLURM_SUBMIT_DIR:-$PWD}}"
RSCRIPT="$SCRIPT_DIR/fab_rescue_test.R"

if [ ! -s "$RSCRIPT" ]; then
    echo "ERROR: cannot find fab_rescue_test.R in $SCRIPT_DIR"
    echo "       set SCRIPT_DIR via: sbatch --export=ALL,SCRIPT_DIR=/path ..."
    exit 1
fi

Rscript "$RSCRIPT"

echo ">>> Done. Outputs in: $OUTPUT_DIR"
