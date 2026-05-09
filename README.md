# Pan-Draft Metabolic Modelling Pipeline for Cow-Rumen SGBs

A SLURM-based pipeline for reconstructing pan-genome-scale metabolic models from species-level genome bins (SGBs) using `gapseq` + pan-Draft, followed by high-throughput phenotyping (auxotrophy detection and substrate preference profiling) and a quality-control battery.

## Overview

This pipeline takes a set of metagenome-assembled genomes (MAGs) belonging to a single SGB (≥30 members recommended) and produces:

1. A gap-filled, FBA-ready pan-metabolic model (`panModel.RDS`)
2. A profile of essential growth factors (auxotrophies)
3. A ranked list of preferred carbon/energy sources
4. A QC report summarising thermodynamic and stoichiometric integrity

The pipeline is designed for batch processing of multiple SGBs on an HPC cluster and produces structured TSV outputs suitable for downstream visualisation.

## Directory Layout

Project root (the directory containing all the scripts):

```
unculturable/
├── pipeline.config             # Central configuration (edit paths here)
├── pan_step1.sh                # Per-MAG draft (gapseq find/draft)
├── pan_step2.sh                # Pan-Draft + medium prediction + gap-fill
├── pan_step3.sh                # Parse stage-1/2 logs into structured TSVs
├── pan_step4.sh                # FBA phenotyping on the pan model
├── run_pan_pipeline.sh         # Wrapper: submits all four steps with dependencies
├── download_mags.sh            # Fetch MAGs from MGnify catalogue
├── leave_one_out.R             # Phenotyping engine (pFBA mode, default)
├── leave_one_out_fba.R         # Phenotyping engine (FBA fallback mode)
├── run_fba_pan.sh              # Submitter for FBA-mode phenotyping
├── qc_battery.R                # QC battery for all pan-models (batch)
├── qc_battery_pilot.R          # QC battery on a single model (verbose, for diagnostics)
├── parse_gapseq_report.py      # Parse one gapseq fill log
├── parse_pan_stage1.py         # Aggregate per-MAG ORF coverage stats
└── output_tsv/                 # All structured outputs land here
    ├── all_reconstruction.tsv
    ├── all_stage1_summary.tsv
    ├── qc_battery_results.tsv
    ├── <SGB>_pan_*.tsv         # Per-SGB phenotyping outputs
    └── merged/                 # Cumulative tables for cross-SGB comparison
```

For each processed SGB, the following per-SGB directories are created:

```
<SGB>_rep/
├── *.faa                       # Per-MAG protein sequences (input)
├── models/                     # Stage-1 per-MAG drafts
│   ├── <MAG>-draft.RDS
│   ├── <MAG>-all-Pathways.tbl
│   └── logs/                   # Stage-1 SLURM stdout/stderr
└── pan_model/                  # Stage-2 pan-Draft outputs
    ├── panModel.RDS            # Final gap-filled pan model
    ├── panModel.xml            # SBML format
    ├── panModel-medium.csv     # Predicted minimal medium
    ├── pan-reactome_stat.tsv   # Core/shell/cloud reaction statistics
    ├── rxnXmod.tsv             # Reaction × MAG presence/absence matrix
    └── logs/                   # Stage-2 SLURM stdout/stderr
```

## Prerequisites

### Software

- `gapseq` v2.0.1 (https://github.com/jotech/gapseq)
- `prodigal` v2.6.3 (for re-predicting ORFs from MAGs)
- `R` ≥ 4.0 with the `cobrar` package
- `Python` ≥ 3.6
- SLURM workload manager (or adapt scripts to your scheduler)

A conda/mamba environment named `gapseq` is assumed:

```bash
mamba create -n gapseq -c bioconda gapseq prodigal r-base
mamba activate gapseq
R -e 'install.packages("cobrar")'
```

### Data

- `nutrients.tsv` (provided with `gapseq`, used for compound name lookup)
- MAG protein FASTA files, one `.faa` per MAG, all in a single directory

## Configuration

Edit `pipeline.config` to match your environment:

```bash
PIPELINE_BASE     # Root directory for all SGB working directories
GAPSEQ            # Path to gapseq executable
NUTRIENTS_TSV     # Path to gapseq's nutrients.tsv
SLURM_PARTITION   # Your SLURM partition name
SLURM_MAIL        # Your email for SLURM notifications
```

All other paths are derived automatically.

## Workflow

### 1. Acquire MAGs (optional helper)

For the MGnify cow-rumen catalogue:

```bash
bash download_mags.sh <SGB_ID>
```

This downloads all cluster member GFFs from the MGnify FTP, extracts the embedded nucleotide FASTA, and re-runs Prodigal to ensure consistent gene prediction across all MAGs.

Output: `<SGB>_rep/*.faa` (one per cluster member).

### 2. Run the full pipeline

For one or more SGBs:

```bash
bash run_pan_pipeline.sh <SGB_ID> [<SGB_ID> ...]
```

This submits four SLURM jobs per SGB, chained by `--dependency=afterok`:

| Stage | Job | Time | Description |
|-------|-----|------|-------------|
| 1 | Array (1..N MAGs) | ~30 min | gapseq find / find-transport / draft per MAG |
| 2 | Single | ~10 min | gapseq pan / medium / fill |
| 3 | Single | ~1 min  | Parse stdout logs into TSV |
| 4 | Single | ~5 min  | FBA phenotyping (LOO + substrate scan) |

Multiple SGBs run in parallel; within each SGB the four stages run sequentially.

### 3. Monitor

```bash
squeue -u $USER
ls output_tsv/
```

### 4. FBA fallback for numerically unstable pan-models

A small fraction of pan-models (typically those with > 1000 reactions and dense pan-reactome connectivity) trigger GLPK simplex degeneracy during the parsimony stage of pFBA. Symptoms include `step 4` running indefinitely or producing millions of solver warnings before terminating. For such models, switch from pFBA to standard FBA using a parallel R engine (`leave_one_out_fba.R`).

Identify candidates by inspecting `output_tsv/logs/stdout_step4_<SGB>_*.txt` for stalled progress, then submit:

```bash
sbatch --export=ALL,SGB=<SGB_ID>,SCRIPT_DIR=$SCRIPT_DIR run_fba_pan.sh
```

The FBA engine writes outputs into the same `output_tsv/` directory using the same naming convention (`<SGB>_pan_*.tsv`) so downstream merge logic and figure code do not need to distinguish between modes. The `Total_Flux` column is `NA` for FBA-mode runs because parsimonious flux is not computed.

For example invocations on three known unstable models in our dataset:

```bash
for SGB in MGYG000291777 MGYG000295175 MGYG000294398; do
    sbatch --export=ALL,SGB=$SGB,SCRIPT_DIR=$SCRIPT_DIR run_fba_pan.sh
done
```

When the merged tables are refreshed, the FBA-mode rows are concatenated alongside the pFBA-mode rows automatically.

### 5. Quality control

After phenotyping, run the QC battery on every reconstructed pan-model in the project:

```bash
Rscript qc_battery.R [<base_dir>] [<out_tsv>]
```

Defaults: `base_dir = $PIPELINE_BASE`, `out_tsv = $PIPELINE_BASE/output_tsv/qc_battery_results.tsv`. The script discovers all `<base_dir>/*_rep/pan_model/panModel.RDS` files and applies three independent tests per model. Sequential runtime is ~1-2 min per model on a typical login node.

For diagnosing a single model in detail (verbose output including sample imbalanced reactions and formula parsing statistics), use the pilot variant:

```bash
Rscript qc_battery_pilot.R <model.RDS>
```

The QC tests are described in detail in the **Methodology Notes** section below. Models passing all three tests are suitable for downstream analysis without further intervention.

## Outputs

### Per-SGB outputs

Located in `output_tsv/<SGB>_pan_*.tsv`:

| File | Content |
|------|---------|
| `<SGB>_pan_summary.tsv`          | One-line model summary (reactions, exchanges, auxotrophies, growth) |
| `<SGB>_pan_auxotrophies.tsv`     | Long-format list of essential growth factors |
| `<SGB>_pan_substrates_full.tsv`  | All ~200 candidate organic substrates with growth/flux |
| `<SGB>_pan_substrates_top15.tsv` | Top 15 substrates by net growth |

### Cumulative tables (one row per SGB)

| File | Content |
|------|---------|
| `output_tsv/all_reconstruction.tsv` | Gap-fill statistics from stage 2 stdout (added reactions, final growth, etc.) |
| `output_tsv/all_stage1_summary.tsv` | Per-MAG ORF coverage aggregated to SGB level (median, range, mean) |
| `output_tsv/qc_battery_results.tsv` | Per-model QC verdicts (ATP-from-water, mass balance, element balance) |

### Merged tables (auto-refreshed at end of stage 4)

| File | Content |
|------|---------|
| `output_tsv/merged/all_summary.tsv`          | All SGBs concatenated |
| `output_tsv/merged/all_auxotrophies.tsv`     | All SGBs concatenated |
| `output_tsv/merged/all_substrates_top15.tsv` | All SGBs concatenated |
| `output_tsv/merged/all_substrates_full.tsv`  | All SGBs concatenated |

## Methodology Notes

### Pan-reactome construction

Reactions present in fewer than 6% of cluster MAGs (default `--min.rxn.freq.in.mods 0.06`) are excluded from the pan-reactome to mitigate noise from MAG fragmentation. This threshold follows De Bernardini et al. (2024).

### Gap-filling medium

Pan-models are gap-filled against minimal media predicted by `gapseq medium`, not against rich media (e.g. TSB). This avoids retention of biosynthetic shortcuts that would otherwise be falsely detected as auxotrophies.

### Phenotyping

The `leave_one_out.R` script performs three procedures:

1. **Auxotrophy detection (LOO)**. Each currently-uptaken organic exchange is set to zero in turn; if predicted growth drops below 10⁻⁴ hr⁻¹, the compound is classified as essential.
2. **Trace-limited background**. All organic exchanges are constrained to −0.05 mmol gDW⁻¹ hr⁻¹; identified auxotrophies are opened to −0.1 mmol gDW⁻¹ hr⁻¹, which is sufficient to satisfy structural availability of essential cofactors without constituting a substantial carbon source.
3. **Substrate scan**. Each candidate organic exchange is opened to −10 in turn; growth and total flux are recorded by parsimonious FBA (pFBA, GLPK solver via `cobrar`).

Net growth is computed as growth on substrate minus trace background growth; substrates with net growth > 0.01 hr⁻¹ are reported as active.

### FBA fallback

`leave_one_out_fba.R` mirrors the three procedures above but uses standard FBA (`fba()`) in place of `pfba()` throughout, including in the substrate scan stage. Growth values are numerically equivalent to the pFBA stage-1 result; the parsimonious flux distribution is not computed, so the `Total_Flux` column is reported as `NA`. This branch is invoked only for models that fail to terminate under pFBA due to GLPK simplex degeneracy (typically 2-10% of SGBs in a typical cohort).

### Quality control

Three independent tests are applied per model. Thresholds and pass/warn/fail classification follow the conventions of the MEMOTE community standard (Lieven et al. 2020), adapted for a non-MILP environment.

1. **ATP synthesis from water**. All exchange reactions are closed (lb = 0) except H₂O (cpd00001) and H⁺ (cpd00067), which are opened to −1000 mmol gDW⁻¹ hr⁻¹. The objective is switched to the ATP maintenance reaction (rxn00062), and FBA is solved. A model passes when the maximum ATP flux falls below 1×10⁻⁶ mmol gDW⁻¹ hr⁻¹, indicating the absence of thermodynamically infeasible cycles capable of generating ATP without organic substrates.

2. **Total mass balance**. For each internal reaction (i.e. reactions other than exchanges, biomass, and demand or sink reactions), the column sum of the stoichiometric matrix is computed weighted by the molecular mass of each metabolite, derived from `met_attr$chemicalFormula`. Reactions involving any metabolite with a missing or unparseable formula are excluded. A reaction passes when the absolute weighted sum is below 0.001 Da. A model passes when ≥99% of checkable internal reactions are balanced.

3. **Per-element balance**. The same internal-reaction set is tested for elemental closure across C, H, N, O, P and S. For each element separately, the absolute element-weighted column sum must remain below 0.001 atom equivalents. A reaction passes when this condition holds simultaneously for all six elements; a model passes when ≥95% of checkable internal reactions are balanced.

The fourth test in the MEMOTE standard (loopless flux contribution) requires mixed-integer programming not natively supported by cobrar 0.2.x and is therefore omitted from this battery; the closely related and stricter ATP-from-water test serves as the primary check for thermodynamically infeasible cycles. A separate cobrapy-based loopless assessment can be added externally if required.

The QC implementation has been validated by manual inspection of a known balanced reaction (rxn00216, ATP:D-glucose 6-phosphotransferase) and by injection of an artificial imbalance that was correctly detected, confirming the parser and column-sum logic. When all reconstructions originate from gapseq and inherit ModelSEED-curated stoichiometry, mass and element balance pass rates of 100% are expected (Henry et al. 2010; Seaver et al. 2021).

### Known limitations

- LOO occasionally classifies internal metabolic intermediates (e.g. chorismate) as essential due to topological gaps in the pan-model. Manual filtering of the auxotrophy list against ModelSEED metabolite categories is recommended before downstream interpretation.
- `gapseq medium` predictions are minimal by design; auxotrophy counts on TSB-fill models will differ from minimal-fill models. Methodological choice must be consistent within a given comparison.
- FBA-mode phenotyping does not produce parsimonious flux distributions. Comparisons of `Total_Flux` between pFBA-mode and FBA-mode SGBs are therefore not meaningful and should be avoided.

## Resuming after failures

### A SLURM step is killed (timeout, OOM, etc.)

Identify the failed step from logs, then resubmit only that step plus downstream steps. For example, when stage 2 times out for `MGYG000291777`:

```bash
SGB=MGYG000291777
SCRIPT_DIR=/path/to/unculturable
LOG2=$SCRIPT_DIR/${SGB}_rep/pan_model/logs
LOG34=$SCRIPT_DIR/output_tsv/logs

# Check whether panModel.RDS was actually produced
ls $SCRIPT_DIR/${SGB}_rep/pan_model/panModel.RDS
tail -5 $LOG2/stdout_step2_*.txt

# If stage 2 is incomplete, resubmit with a longer time limit
JOB2=$(sbatch \
       --time=12:00:00 --partition=k2-medpri \
       --output=$LOG2/stdout_step2_%j.txt \
       --error=$LOG2/stderr_step2_%j.txt \
       --export=ALL,SGB=$SGB,SCRIPT_DIR=$SCRIPT_DIR \
       $SCRIPT_DIR/pan_step2.sh | grep -oP '\d+$')

# Chain stages 3 and 4
JOB3=$(sbatch \
       --dependency=afterok:$JOB2 \
       --output=$LOG34/stdout_step3_${SGB}_%j.txt \
       --error=$LOG34/stderr_step3_${SGB}_%j.txt \
       --export=ALL,SGB=$SGB,SCRIPT_DIR=$SCRIPT_DIR \
       $SCRIPT_DIR/pan_step3.sh | grep -oP '\d+$')

sbatch \
    --dependency=afterok:$JOB3 \
    --output=$LOG34/stdout_step4_${SGB}_%j.txt \
    --error=$LOG34/stderr_step4_${SGB}_%j.txt \
    --export=ALL,SGB=$SGB,SCRIPT_DIR=$SCRIPT_DIR \
    $SCRIPT_DIR/pan_step4.sh
```

### Stage 4 stalls under pFBA

If `stdout_step4_<SGB>_*.txt` shows no progress for hours and the model is not unusually large (under 1500 reactions), GLPK simplex degeneracy is the likely cause. Switch to FBA mode:

```bash
sbatch --export=ALL,SGB=<SGB_ID>,SCRIPT_DIR=$SCRIPT_DIR run_fba_pan.sh
```

This replaces stage 4 only; stages 1-3 are not re-run.

### Re-running an SGB from scratch

```bash
SGB=MGYG000291777
rm -rf $PIPELINE_BASE/${SGB}_rep/models/ $PIPELINE_BASE/${SGB}_rep/pan_model/

# Remove the SGB's row from cumulative tables
for f in output_tsv/all_reconstruction.tsv output_tsv/all_stage1_summary.tsv \
         output_tsv/qc_battery_results.tsv; do
    grep -v "${SGB}_pan" "$f" > "$f.tmp" && mv "$f.tmp" "$f"
done

bash run_pan_pipeline.sh $SGB
```

## Citation

If you use this pipeline, please cite:

- gapseq: Zimmermann et al. (2021), *Genome Biology*, 22:81.
- pan-Draft: De Bernardini et al. (2024), *Genome Biology*, 25(1):280.
- cobrar: Waschina, Zimmermann & Froitzheim (2026), R package v0.2.3, https://github.com/Waschina/cobrar
- MGnify catalogue: Gurbich et al. (2023), *J Mol Biol*, 435(14):168016.
- ModelSEED biochemistry: Henry et al. (2010), *Nat Biotechnol*, 28(9):977-982; Seaver et al. (2021), *Nucleic Acids Res*, 49(D1):D575-D588.
- MEMOTE QC standard: Lieven et al. (2020), *Nat Biotechnol*, 38(3):272-276.
