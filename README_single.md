# Single-Genome Metabolic Modelling Pipeline

A SLURM-based pipeline for reconstructing genome-scale metabolic models from individual bacterial genomes (reference strains, isolates, or single MAGs) using `gapseq`, followed by high-throughput phenotyping (auxotrophy detection and substrate preference profiling).

This pipeline is the single-genome counterpart of the pan-Draft pipeline. The two are designed to produce directly comparable outputs; see also `README.md` for the pan-Draft pipeline.

## Overview

This pipeline takes a directory of `.faa` files (one per organism) and produces, for each:

1. A gap-filled, FBA-ready metabolic model (`<ID>.RDS`)
2. A profile of essential growth factors (auxotrophies)
3. A ranked list of preferred carbon/energy sources

All four `gapseq` stages (find, find-transport, draft, fill) plus phenotyping run as a single SLURM array job, with one task per input genome.

## Directory Layout

Scripts are kept alongside the pan-Draft scripts in the project root:

```
unculturable/
├── pipeline_single.config       # Single-genome pipeline configuration
├── single_step.sh               # Full per-genome pipeline (SLURM array task)
├── run_single_batch.sh          # Wrapper: submits the array job
├── leave_one_out.R              # Phenotyping engine (shared with pan pipeline)
└── parse_gapseq_report.py       # Log parser (shared with pan pipeline)
```

Inputs and outputs live in a separate workspace (configurable via `SINGLE_BASE`):

```
single_genome/
├── faa/                         # Place your input .faa files here
│   ├── E_coli.faa
│   ├── K_pneumoniae.faa
│   ├── B_fragilis.faa
│   ├── *-draft.RDS              # gapseq intermediates land here too
│   ├── *.RDS                    # Final per-genome models
│   └── logs/                    # SLURM stdout/stderr per array task
└── output_tsv/                  # Structured outputs
    ├── all_reconstruction.tsv
    ├── <ID>_*.tsv               # Per-genome phenotyping outputs
    └── merged/                  # Cumulative tables across all genomes
```

## Prerequisites

### Software

- `gapseq` v2.0.1 (https://github.com/jotech/gapseq)
- `R` ≥ 4.0 with the `cobrar` package
- `Python` ≥ 3.6
- SLURM workload manager

A conda/mamba environment named `gapseq` is assumed:

```bash
mamba create -n gapseq -c bioconda gapseq r-base
mamba activate gapseq
R -e 'install.packages("cobrar")'
```

### Data

- `nutrients.tsv` (provided with `gapseq`, used for compound name lookup)
- One `.faa` file per genome, all placed in a single input directory

### Naming convention

Each `.faa` filename (without extension) becomes the organism ID used throughout all outputs (TSV column values, phenotyping file prefixes, etc.). To avoid downstream parsing issues, **do not use dots in genome names**:

```bash
# Bad:  E.coli.faa  → ID becomes "E.coli", may break R/awk parsing
# Good: E_coli.faa  → ID becomes "E_coli"
```

## Configuration

Edit `pipeline_single.config` to match your environment:

```bash
SINGLE_BASE       # Root directory for inputs and outputs
GAPSEQ            # Path to gapseq executable
NUTRIENTS_TSV     # Path to gapseq's nutrients.tsv
SLURM_PARTITION   # Your SLURM partition name
SLURM_MAIL        # Your email for SLURM notifications
```

All other paths are derived automatically.

## Workflow

### 1. Place input genomes

```bash
ls $SINGLE_BASE/faa/*.faa
```

Each file should be a valid protein FASTA. Filenames define the organism IDs.

### 2. Run the pipeline

```bash
bash run_single_batch.sh [<faa_dir>]
```

The argument is optional; if omitted, the wrapper uses `$SINGLE_BASE/faa`.

This submits a single SLURM array job with one task per `.faa` file. Each task runs the full pipeline independently, so the wall-clock time equals that of the slowest genome.

| Step within each task | Time | Description |
|------------------------|------|-------------|
| `gapseq find`          | ~5–10 min | Pathway and reaction prediction |
| `gapseq find-transport`| ~1 min    | Transporter prediction |
| `gapseq draft`         | < 1 min   | Draft network assembly |
| `gapseq medium`        | < 1 min   | Predict minimal growth medium |
| `gapseq fill`          | ~2–5 min  | Gap-filling |
| Log parsing            | < 1 min   | Extract reconstruction stats |
| Phenotyping (R)        | ~3–10 min | LOO + substrate scan |
| Cumulative table merge | < 1 min   | Refresh `merged/` outputs |

### 3. Monitor

```bash
squeue -u $USER
ls $SINGLE_BASE/output_tsv/
```

## Outputs

### Per-genome outputs

Located in `$SINGLE_BASE/output_tsv/<ID>_*.tsv`:

| File | Content |
|------|---------|
| `<ID>_summary.tsv`          | One-line model summary (reactions, exchanges, auxotrophies, growth) |
| `<ID>_auxotrophies.tsv`     | Long-format list of essential growth factors |
| `<ID>_substrates_full.tsv`  | All ~200 candidate organic substrates with growth/flux |
| `<ID>_substrates_top15.tsv` | Top 15 substrates by net growth |

### Cumulative tables (one row per genome, appended)

| File | Content |
|------|---------|
| `output_tsv/all_reconstruction.tsv` | Gap-fill statistics from each genome's stdout (ORF coverage, added reactions, final growth, etc.) |

### Merged tables (auto-refreshed at end of each task)

| File | Content |
|------|---------|
| `output_tsv/merged/all_summary.tsv`          | All genomes concatenated |
| `output_tsv/merged/all_auxotrophies.tsv`     | All genomes concatenated |
| `output_tsv/merged/all_substrates_top15.tsv` | All genomes concatenated |
| `output_tsv/merged/all_substrates_full.tsv`  | All genomes concatenated |

## Methodology Notes

### Gap-filling medium

Each genome is gap-filled against a minimal medium predicted by `gapseq medium`, identical to the pan-Draft pipeline. This ensures direct comparability between single-genome and pan-Draft phenotyping outputs (auxotrophy counts, substrate preferences) when both pipelines are applied to the same project.

If the project requires comparison against rich-medium reconstructions (e.g. TSB), modify the `gapseq fill` line in `single_step.sh` to use `-n /path/to/TSBmed.csv` and document the choice in your Methods section.

### Phenotyping

See the pan-Draft pipeline README (`README.md`) for a full description of `leave_one_out.R`. The procedure is identical for single-genome and pan models.

### Combining outputs with the pan-Draft pipeline

When both pipelines have been run on the same project, the cumulative tables can be combined for cross-comparison:

```bash
mkdir -p combined_output
for kind in summary auxotrophies substrates_top15 substrates_full; do
    awk 'FNR==1 && NR!=1 {next} {print}' \
        unculturable/output_tsv/merged/all_${kind}.tsv \
        single_genome/output_tsv/merged/all_${kind}.tsv \
        > combined_output/all_${kind}.tsv
done
```

The `organism` column in each TSV uses the suffix `_pan` for pan-Draft models, while single-genome IDs carry no suffix. Distinguish groups in downstream analysis with a regex on this column.

## Resuming after failures

### A single task in the array failed

Identify the failed task from the SLURM logs, then resubmit just that index:

```bash
# Find which array index failed (e.g. task 3 in job 8612345)
ls -lh $SINGLE_BASE/faa/logs/stderr_8612345_*.txt | sort -k5 -n -r

# Resubmit only that index
sbatch --array=3 \
    --partition=$SLURM_PARTITION \
    --output=$SINGLE_BASE/faa/logs/stdout_%A_%a.txt \
    --error=$SINGLE_BASE/faa/logs/stderr_%A_%a.txt \
    --export=ALL,WORKDIR=$SINGLE_BASE/faa,SCRIPT_DIR=/path/to/unculturable \
    /path/to/unculturable/single_step.sh
```

### Re-running a single genome from scratch

```bash
ID=E_coli
cd $SINGLE_BASE/faa
rm -f ${ID}-* ${ID}.RDS ${ID}.xml

# Remove the genome's row from the cumulative reconstruction table
grep -v "^${ID}\b" $SINGLE_BASE/output_tsv/all_reconstruction.tsv > tmp \
    && mv tmp $SINGLE_BASE/output_tsv/all_reconstruction.tsv

# Remove the per-genome phenotyping TSVs
rm -f $SINGLE_BASE/output_tsv/${ID}_*.tsv

# Identify which array index this genome is at, then resubmit
ls $SINGLE_BASE/faa/*.faa | nl   # find the line number = array index
sbatch --array=<index> ...        # as above
```

## Citation

If you use this pipeline, please cite:

- gapseq: Zimmermann et al. (2021), *Genome Biology*, 22:81.
- cobrar: Waschina, Zimmermann & Froitzheim (2026), R package v0.2.3, https://github.com/Waschina/cobrar 

