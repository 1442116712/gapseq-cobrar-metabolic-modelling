# GEM-Guided Cultivation Modelling — Code Pipeline

Two complementary genome-scale metabolic modelling (GEM) pipelines used in this manuscript, plus the shared validation and aggregation tooling that sits on top of them. Both pipelines turn protein FASTAs into FBA-ready models, run a phenotyping pass (auxotrophy detection + substrate preference profiling), and feed into the same downstream QC, rich-medium LOO and tier-classification framework.

The README is split into three parts:

- **Part I — Core modelling system** — what you need to actually build models and run the standard minimal-medium phenotyping. Cover the two operating modes (single-genome and pan-Draft) and the shared phenotyping engine.
- **Part II — Validation & confidence stratification** — QC battery, rich-medium LOO, 4-tier auxotrophy classification, and the predicted-media aggregator. These do not produce new models; they assess and stratify what the core pipeline outputs.
- **Part III — Operations** — monitoring, resuming after failures, re-running from scratch.

If you just want to use the pipeline, Part I + the operations section of Part III is the canonical path. Part II is for evaluating output quality.

---

## Repository layout

```
Manuscript4/
├── README.md                  # this file
│
├── single/                    # Single-genome pipeline (Part I.A)
│   ├── README_single.md       # pointer back to this README
│   ├── pipeline_single.config # config for the single-genome pipeline
│   ├── single_step.sh         # SLURM array task: gapseq + minimal-medium phenotyping
│   ├── single_step_rich.sh    # SLURM array task: rich-medium LOO (Part II.B)
│   └── run_single_batch.sh    # wrapper: submits the array job
│
├── pan-draft/                 # Pan-Draft pipeline (Part I.B) + shared tooling
│   ├── README.md              # pointer back to this README
│   ├── pipeline.config        # config for the pan-Draft pipeline
│   ├── pan_step1.sh           # stage 1 — per-MAG draft (gapseq find/draft)
│   ├── pan_step2.sh           # stage 2 — pan-Draft + medium prediction + gap-fill
│   ├── pan_step3.sh           # stage 3 — parse stdout logs into structured TSVs
│   ├── pan_step4.sh           # stage 4 — FBA phenotyping on minimal medium
│   ├── pan_step5_rich.sh      # stage 5 — rich-medium LOO (Part II.B)
│   ├── run_pan_pipeline.sh    # wrapper: submits stages 1-4 with dependencies
│   ├── run_rich_loo_all.sh    # dispatcher for stage 5 (pan and single modes)
│   ├── run_fba_pan.sh         # submitter for FBA-mode phenotyping (fallback)
│   ├── download_mags.sh       # fetch MAGs from MGnify catalogue
│   ├── extract_media.sh       # aggregate gapseq-predicted minimal media across all models (Part II.D)
│   │
│   ├── leave_one_out.R        # shared phenotyping engine (pFBA mode, default)
│   ├── leave_one_out_fba.R    # shared phenotyping engine (FBA fallback)
│   ├── leave_one_out_rich.R   # shared rich-medium LOO engine (Part II.B)
│   ├── qc_battery.R           # QC battery for pan-models (Part II.A)
│   ├── qc_battery_pilot.R     # verbose QC for a single model (diagnostics)
│   │
│   ├── parse_gapseq_report.py # parse one gapseq fill log
│   └── parse_pan_stage1.py    # aggregate per-MAG ORF coverage stats
│
└── R code/                    # Manuscript figures (downstream consumers of the modelling outputs)
    ├── Fig3_tier_classification.R
    └── ...
```

The R phenotyping engines, the rich-LOO engine and the predicted-media aggregator live alongside the pan-Draft scripts but are invoked by both pipelines.

---

## Prerequisites

### Software

- `gapseq` v2.0.1 (https://github.com/jotech/gapseq)
- `prodigal` v2.6.3 (only needed for the pan-Draft pipeline, to re-predict ORFs from MAG nucleotide FASTAs)
- `R` ≥ 4.0 with the `cobrar` package
- `Python` ≥ 3.6
- SLURM workload manager (or adapt the scripts to your scheduler)

A conda/mamba environment named `gapseq` is assumed:

```bash
mamba create -n gapseq -c bioconda gapseq prodigal r-base
mamba activate gapseq
R -e 'install.packages("cobrar")'
```

### Data

- `nutrients.tsv` (provided with `gapseq`, used for compound-name lookup in the phenotyping engines)
- For single-genome mode: one `.faa` per genome
- For pan-Draft mode: MAG protein FASTAs, one `.faa` per cluster member, grouped by SGB

### Naming convention

Each `.faa` filename (without extension) becomes the organism ID used throughout the outputs (TSV column values, phenotyping file prefixes, etc.). To avoid downstream parsing issues, **do not use dots in genome names**:

```bash
# Bad:  E.coli.faa  → ID becomes "E.coli", may break R/awk parsing
# Good: E_coli.faa  → ID becomes "E_coli"
```

### Configuration

Each pipeline has its own config file with paths to the executables, your SLURM partition, and the location of `nutrients.tsv`.

- `single/pipeline_single.config` — `SINGLE_BASE`, `GAPSEQ`, `NUTRIENTS_TSV`, `SLURM_PARTITION`, `SLURM_MAIL`
- `pan-draft/pipeline.config` — `PIPELINE_BASE`, `GAPSEQ`, `NUTRIENTS_TSV`, `SLURM_PARTITION`, `SLURM_MAIL`, plus `LEAVE_ONE_OUT_R` and `LEAVE_ONE_OUT_RICH_R` pointers

All other paths derive automatically.

---

# Part I — Core modelling system

This part covers everything needed to **build** models and run the canonical minimal-medium phenotyping. Both pipelines produce directly comparable outputs (same column schema, same phenotyping engine), so single-genome and pan-Draft TSVs can be merged for cross-comparison without further adjustment.

## I.A — Single-genome mode (`single/`)

Use this mode for **reference strains, isolates, ARG-MAGs and any single-genome input**, including the in-vitro panel models analysed in the manuscript (DSM 11370 / GCF_000426565.1 and the three SH strains).

### Inputs

- A directory of `.faa` files, one per organism, placed in `$SINGLE_BASE/faa/` (lab / in-vitro models) or `$SINGLE_BASE/faa_insilico/` (references, ARG-MAGs, SGB representatives — optional secondary input set).
- Each `.faa` should be a valid protein FASTA.

### Submission

```bash
bash run_single_batch.sh [<faa_dir>]
```

The argument is optional; if omitted, the wrapper uses `$SINGLE_BASE/faa`. The wrapper submits a single SLURM array job with one task per `.faa`. Each task runs the full pipeline independently.

### Per-task stages

| Step within each task | Time | Description |
|------------------------|------|-------------|
| `gapseq find`          | ~5–10 min | Pathway and reaction prediction (HMM-based) |
| `gapseq find-transport`| ~1 min    | Transporter prediction |
| `gapseq draft`         | < 1 min   | Draft network assembly |
| `gapseq medium`        | < 1 min   | Predict minimal growth medium → `<ID>-medium.csv` |
| `gapseq fill`          | ~2–5 min  | Gap-fill against the predicted minimal medium |
| Log parsing            | < 1 min   | Extract reconstruction stats via `parse_gapseq_report.py` |
| Phenotyping (R)        | ~3–10 min | Minimal-medium LOO + substrate scan via `leave_one_out.R` |
| Cumulative table merge | < 1 min   | Refresh `merged/` outputs |

Wall-clock equals the slowest genome in the array (typically 15–25 min).

### Per-genome outputs

Located in `$SINGLE_BASE/output_tsv/<ID>_*.tsv`:

| File | Content |
|------|---------|
| `<ID>_summary.tsv`            | One-line model summary (reactions, exchanges, auxotrophies, growth) |
| `<ID>_auxotrophies.tsv`       | Long-format list of essential growth factors (minimal-medium LOO) |
| `<ID>_substrates_full.tsv`    | All ~200 candidate organic substrates with growth/flux |
| `<ID>_substrates_top15.tsv`   | Top 15 substrates by net growth |

The gapseq-predicted minimal medium (`<ID>-medium.csv`) sits next to the model RDS in `$SINGLE_BASE/faa/` (not in `output_tsv/`).

### Cumulative & merged tables

| File | Content |
|------|---------|
| `output_tsv/all_reconstruction.tsv`              | Per-genome gap-fill statistics (one row per genome, appended) |
| `output_tsv/merged/all_summary.tsv`              | All genomes concatenated |
| `output_tsv/merged/all_auxotrophies.tsv`         | All genomes concatenated |
| `output_tsv/merged/all_substrates_top15.tsv`     | All genomes concatenated |
| `output_tsv/merged/all_substrates_full.tsv`      | All genomes concatenated |

## I.B — Pan-Draft mode (`pan-draft/`)

Use this mode for **MAG cohorts** belonging to a single species-level genome bin (SGB), ≥ 30 members recommended. Produces a gap-filled pan-metabolic model (`panModel.RDS`) that represents the SGB's combined reaction repertoire.

### Inputs

- One MAG protein FASTA per cluster member, all in `<SGB>_rep/*.faa`. Use the helper `download_mags.sh` for MGnify catalogues (see step 1 below).

### Workflow

Optionally fetch MAGs:

```bash
bash download_mags.sh <SGB_ID>
```

Then submit the full pipeline for one or more SGBs:

```bash
bash run_pan_pipeline.sh <SGB_ID> [<SGB_ID> ...]
```

This submits **four** SLURM jobs per SGB, chained by `--dependency=afterok`. Multiple SGBs run in parallel; within each SGB the four stages run sequentially.

| Stage | Job              | Time     | Description |
|-------|------------------|----------|-------------|
| 1     | Array (1..N MAGs)| ~30 min  | gapseq find / find-transport / draft per MAG |
| 2     | Single           | ~10 min  | gapseq pan / medium / fill → `panModel.RDS`, `panModel-medium.csv` |
| 3     | Single           | ~1 min   | Parse stdout logs into TSV via `parse_pan_stage1.py` and `parse_gapseq_report.py` |
| 4     | Single           | ~5 min   | FBA phenotyping on minimal medium (LOO + substrate scan) via `leave_one_out.R` |

Monitor:

```bash
squeue -u $USER
ls output_tsv/
```

### FBA fallback for numerically unstable pan-models

A small fraction of pan-models (typically those with > 1000 reactions and dense pan-reactome connectivity) trigger GLPK simplex degeneracy during the parsimony stage of pFBA. Symptoms include `step 4` running indefinitely or producing millions of solver warnings before terminating. For such models, switch from pFBA to standard FBA using the parallel R engine `leave_one_out_fba.R`:

```bash
sbatch --export=ALL,SGB=<SGB_ID>,SCRIPT_DIR=$SCRIPT_DIR run_fba_pan.sh
```

The FBA engine writes outputs into the same `output_tsv/` directory using the same naming convention (`<SGB>_pan_*.tsv`) so downstream merge logic and figure code do not need to distinguish between modes. The `Total_Flux` column is `NA` for FBA-mode runs because parsimonious flux is not computed.

Three SGBs in our dataset required this fallback:

```bash
for SGB in MGYG000291777 MGYG000295175 MGYG000294398; do
    sbatch --export=ALL,SGB=$SGB,SCRIPT_DIR=$SCRIPT_DIR run_fba_pan.sh
done
```

When the merged tables are refreshed, FBA-mode rows are concatenated alongside the pFBA-mode rows automatically.

### Per-SGB directory layout

For each processed SGB the pipeline creates:

```
<SGB>_rep/
├── *.faa                       # per-MAG protein sequences (input)
├── models/                     # stage-1 per-MAG drafts
│   ├── <MAG>-draft.RDS
│   ├── <MAG>-all-Pathways.tbl
│   └── logs/                   # stage-1 SLURM stdout/stderr
└── pan_model/                  # stage-2 pan-Draft outputs
    ├── panModel.RDS            # final gap-filled pan model
    ├── panModel.xml            # SBML format
    ├── panModel-medium.csv     # predicted minimal medium
    ├── pan-reactome_stat.tsv   # core/shell/cloud reaction statistics
    ├── rxnXmod.tsv             # reaction × MAG presence/absence matrix
    └── logs/                   # stage-2 SLURM stdout/stderr
```

### Per-SGB outputs

Located in `output_tsv/<SGB>_pan_*.tsv`:

| File | Content |
|------|---------|
| `<SGB>_pan_summary.tsv`             | One-line model summary (reactions, exchanges, auxotrophies, growth) |
| `<SGB>_pan_auxotrophies.tsv`        | Long-format list of essential growth factors (minimal-medium LOO) |
| `<SGB>_pan_substrates_full.tsv`     | All ~200 candidate organic substrates with growth/flux |
| `<SGB>_pan_substrates_top15.tsv`    | Top 15 substrates by net growth |

### Cumulative & merged tables

| File | Content |
|------|---------|
| `output_tsv/all_reconstruction.tsv`              | Gap-fill statistics from stage-2 stdout (added reactions, final growth, etc.) |
| `output_tsv/all_stage1_summary.tsv`              | Per-MAG ORF coverage aggregated to SGB level (median, range, mean) |
| `output_tsv/merged/all_summary.tsv`              | All SGBs concatenated |
| `output_tsv/merged/all_auxotrophies.tsv`         | All SGBs concatenated |
| `output_tsv/merged/all_substrates_top15.tsv`     | All SGBs concatenated |
| `output_tsv/merged/all_substrates_full.tsv`     | All SGBs concatenated |

## I.C — Shared phenotyping engine (`leave_one_out.R`)

Both pipelines invoke `leave_one_out.R` after gap-fill. The script runs three procedures in sequence:

1. **Auxotrophy detection (LOO)**. Each currently-uptaken organic exchange is closed in turn; if predicted growth drops below 10⁻⁴ h⁻¹, the compound is classified as essential.
2. **Trace-limited background**. All organic exchanges are constrained to −0.05 mmol gDW⁻¹ h⁻¹; identified auxotrophies are reopened to −0.1 mmol gDW⁻¹ h⁻¹ (sufficient to satisfy structural availability of essential cofactors without constituting a substantial carbon source).
3. **Substrate scan**. Each candidate organic exchange is opened to −10 in turn; growth and total flux are recorded by parsimonious FBA (pFBA, GLPK solver via `cobrar`).

Net growth = growth on substrate − trace background growth. Substrates with net growth > 0.01 h⁻¹ are reported as active.

`leave_one_out_fba.R` mirrors the three procedures above but uses standard FBA (`fba()`) instead of `pfba()`, including in the substrate scan. Growth values are numerically equivalent to the pFBA stage-1 result; the parsimonious flux distribution is not computed, so the `Total_Flux` column is reported as `NA`. This branch is invoked only for the three pan-models that fail to terminate under pFBA (see Part I.B).

### Universal inorganic exchanges

The phenotyping engine excludes a fixed set of eight inorganic exchanges from the candidate LOO set: H₂O (cpd00001), H⁺ (cpd00067), phosphate (cpd00009), ammonium (cpd00013), sulfate (cpd00048), K⁺ (cpd00205), Mg²⁺ (cpd00254), Na⁺ (cpd00971). These are treated as universally available and are not tested for auxotrophy.

### Combining single + pan outputs

When both pipelines have been run on the same project, cumulative tables can be combined for cross-comparison:

```bash
mkdir -p combined_output
for kind in summary auxotrophies substrates_top15 substrates_full; do
    awk 'FNR==1 && NR!=1 {next} {print}' \
        pan-draft/output_tsv/merged/all_${kind}.tsv \
        single/output_tsv/merged/all_${kind}.tsv \
        > combined_output/all_${kind}.tsv
done
```

The `organism` column carries the suffix `_pan` for pan-Draft models; single-genome IDs carry no suffix. Distinguish groups in downstream analysis with a regex on this column.

---

# Part II — Validation & confidence stratification

The validation layer **does not produce new models** — it evaluates and stratifies what the core pipeline already built. Use it once a cohort of models has been reconstructed; rerun it whenever the model set changes.

## II.A — QC battery (`qc_battery.R`)

Three independent tests are applied per model. Thresholds and pass/warn/fail follow the conventions of the MEMOTE community standard (Lieven et al. 2020), adapted for a non-MILP environment.

1. **ATP synthesis from water**. All exchange reactions are closed (lb = 0) except H₂O (cpd00001) and H⁺ (cpd00067), opened to −1000 mmol gDW⁻¹ h⁻¹. The objective is switched to the ATP maintenance reaction (rxn00062), and FBA is solved. A model passes when the maximum ATP flux falls below 1×10⁻⁶ mmol gDW⁻¹ h⁻¹ — i.e. no thermodynamically infeasible cycles generate ATP without organic substrates.

2. **Total mass balance**. For each internal reaction (i.e. non-exchange, non-biomass, non-demand/sink), the column sum of the stoichiometric matrix is computed weighted by molecular mass derived from `met_attr$chemicalFormula`. Reactions involving metabolites with missing or unparseable formulas are excluded. A reaction passes when the absolute weighted sum is below 0.001 Da; a model passes when ≥ 99 % of checkable internal reactions are balanced.

3. **Per-element balance**. The same internal-reaction set is tested for elemental closure across C, H, N, O, P, S. For each element separately, the absolute element-weighted column sum must remain below 0.001 atom equivalents. A reaction passes when this holds for all six elements; a model passes when ≥ 95 % of checkable internal reactions are balanced.

The fourth MEMOTE test (loopless flux contribution) requires MILP not natively supported by `cobrar` 0.2.x and is omitted; the stricter ATP-from-water test is the primary thermodynamic check.

### How to run

Batch mode (every reconstructed pan-model in the project):

```bash
Rscript qc_battery.R [<base_dir>] [<out_tsv>]
```

Defaults: `base_dir = $PIPELINE_BASE`, `out_tsv = $PIPELINE_BASE/output_tsv/qc_battery_results.tsv`. The script auto-discovers `<base_dir>/*_rep/pan_model/panModel.RDS`. Sequential runtime is ~1–2 min per model on a typical login node.

Pilot mode (verbose, single model, for diagnostics):

```bash
Rscript qc_battery_pilot.R <model.RDS>
```

Output: `output_tsv/qc_battery_results.tsv` (one row per model; ATP-from-water flux, mass-balance pass rate, element-balance pass rate, overall verdict).

## II.B — Rich-medium LOO (`leave_one_out_rich.R`)

After steps 4/5 have produced the canonical minimal-medium auxotrophy lists, run a second leave-one-out scan in which every organic exchange is opened to `RICH_LB` (default −10 mmol gDW⁻¹ h⁻¹) before each compound is dropped. The contrast between rich-LOO and minimal-LOO outcomes is the basis for the Tier A vs Tier C split.

### Pan-Draft submission

```bash
# Stable cohort (default pFBA)
bash run_rich_loo_all.sh pan <SGB_ID> [<SGB_ID> ...]

# FBA mode for the three unstable models in our dataset
USE_FBA=TRUE bash run_rich_loo_all.sh pan MGYG000291777 MGYG000294398 MGYG000295175
```

Each submission triggers `pan_step5_rich.sh`, a single SLURM job (~5–10 min per model) that:

1. Loads `panModel.RDS`
2. Opens every organic exchange to `RICH_LB`
3. Solves baseline growth on the rich medium and confirms `> 1e-4 h⁻¹`
4. Loops over every organic exchange, closing one at a time, recording rich-LOO growth
5. Reads `<SGB>_pan_auxotrophies.tsv` to cross-reference with the minimal LOO
6. Writes per-compound and per-model summary TSVs into `output_tsv/`
7. Appends to the cumulative merged tables

### Single-genome submission

```bash
# Single-genome models (references, ARG-MAGs, in-vitro)
bash run_rich_loo_all.sh single <faa_dir_full_path>

# Or directly via the array submitter:
N=$(ls $SINGLE_BASE/faa/*.faa | wc -l)
sbatch --array=1-$N \
    --partition=$SLURM_PARTITION \
    --output=$SINGLE_BASE/faa/logs/stdout_%A_%a.txt \
    --error=$SINGLE_BASE/faa/logs/stderr_%A_%a.txt \
    --export=ALL,WORKDIR=$SINGLE_BASE/faa,SCRIPT_DIR=$PWD \
    single_step_rich.sh
```

`single_step_rich.sh` is the per-genome counterpart of `pan_step5_rich.sh`. Same env vars (`USE_FBA`, `RICH_LB`, `MODEL_FILE`, `OUTPUT_PREFIX`, `OUTPUT_DIR`, `NUTRIENTS_TSV`), same outputs.

### Outputs (both modes)

| File | Content |
|------|---------|
| `<ID>_rich_LOO.tsv`           | Per-compound rich-LOO outcome with intermediate labels |
| `<ID>_rich_LOO_summary.tsv`   | Per-organism summary (counts, base rich growth, solver mode) |
| `merged/all_rich_LOO.tsv`     | All models concatenated |
| `merged/all_rich_LOO_summary.tsv` | All summaries concatenated |

### Solver fallback

The same three pan-models that fail pFBA in stage 4 also fail under pFBA in the rich-medium LOO. Submit them with `USE_FBA=TRUE`. Rich-LOO outputs use the same column schema regardless of solver mode (the `Total_Flux` column is not produced by stage 5, since the rich LOO only needs growth/no-growth).

### Intermediate labels emitted by the rich-LOO engine

`leave_one_out_rich.R` itself emits only the two-way cross of (rich-essential × min-essential):

| Intermediate label | rich-essential | min-essential | Notes |
|--------------------|----------------|----------------|-------|
| `A_hard`           | ✓ | ✓ | feeds Tier A |
| `C_network_coupled`| ✗ | ✓ | feeds Tier C |
| `X_rich_only`      | ✓ | ✗ | sanity flag — should be zero under FBA monotonicity |
| `not_essential`    | ✗ | ✗ | candidate for Tier B if also in the predicted medium |

Tier B (gap-fill bridged) and Tier D (anomaly) require the third axis — the gapseq-predicted minimal medium — and are assembled in Part II.C below using the aggregator in Part II.D.

## II.C — 4-tier classification (A / B / C / D)

The classification used throughout the manuscript (Figure 3; Supplementary Table S10a_Auxotrophy_tier) is built by cross-referencing **three** per-(model × compound) signals:

1. **rich-LOO essential** — from `<ID>_rich_LOO.tsv` (Part II.B)
2. **minimal-LOO essential** — from `<ID>_auxotrophies.tsv` (Part I.C)
3. **in-medium** — whether the compound appears in the gapseq-predicted minimal medium for that model, from `all_predicted_media.tsv` (Part II.D)

| Tier | in_medium | min-LOO essential | rich-LOO essential | Code label                | Interpretation |
|------|-----------|-------------------|--------------------|----------------------------|----------------|
| **A** | ✓ | ✓ | ✓ | `A_hard`                | Robust auxotroph: biosynthesis genuinely lacking in the model |
| **B** | ✓ | ✗ | ✗ | `B_gap_fill_rescued`    | Gap-fill bridged a synthesis route (compound not in the essential-factor list) |
| **C** | ✓ | ✓ | ✗ | `C_network_coupled`     | Condition-dependent: synthesis enabled when other organic precursors are simultaneously available |
| **D** | ✗ | ✓ | (any) | `D_anomaly`           | Anomaly: essential in minimal LOO but absent from the predicted medium (rare) |

The reverse case (`X_rich_only`: rich-essential but not min-essential) should be zero under FBA monotonicity and is excluded from the final 4-tier table; it is reported by the rich-LOO engine only as a sanity flag.

For downstream analysis (e.g. Figure 3), Tier A is the high-confidence essential-factor set; Tier C is flagged as condition-dependent and should be interpreted with caution when designing defined media.

> **Note.** The current build assembles the S10a_Auxotrophy_tier table by joining the three inputs in the supplementary spreadsheet build; no standalone R/Python script is checked in for the join. The inputs are sufficient to reproduce it: see the column schema in `Fig3_tier_classification.R` for the expected layout.

## II.D — Predicted-media aggregator (`extract_media.sh`)

`gapseq medium` (invoked inside both pipelines) produces, for every reconstructed model, a `*-medium.csv` describing the per-compound minimal medium that gap-fill targeted (cpd_id, compound name, max uptake flux). `extract_media.sh` walks every model directory in the project — pan models, single-genome references, ARG-MAGs and lab/in-vitro single-genome models — and concatenates these per-model CSVs into a single TSV with a `model_type` column for downstream filtering:

```bash
cd pan-draft/
bash extract_media.sh
```

### Inputs (hard-coded at the top of the script)

| Glob | `model_type` |
|------|--------------|
| `$PAN_BASE/MGYG*_rep/pan_model/panModel-medium.csv` | `pan` |
| `$SINGLE_BASE/faa_insilico/<ID>-medium.csv` where `<ID>` starts with `Ecoli_`, `Btheta_` or `Efaecalis_` | `reference` |
| `$SINGLE_BASE/faa_insilico/<ID>-medium.csv` where `<ID>` ∈ {`MGYG000292883`, `MGYG000295553`} | `arg_mag` |
| Other `$SINGLE_BASE/faa_insilico/MGYG*-medium.csv` | `sgb_rep` |
| `$SINGLE_BASE/faa/<ID>-medium.csv` | `lab` |

Edit `PAN_BASE` and `SINGLE_BASE` at the top of `extract_media.sh` if your directory layout differs from the one assumed in `pipeline.config` / `pipeline_single.config`.

### Output

```
$PAN_BASE/output_tsv/all_predicted_media.tsv
columns: model  model_type  cpd_id  compound_name  max_flux
```

This is the third input to the 4-tier classification described in Part II.C (the `in_medium` column for each (model × compound) pair is derived by joining `all_predicted_media.tsv` against the per-model rich + minimal LOO outputs).

---

# Part III — Operations

## Monitoring

```bash
squeue -u $USER
ls pan-draft/output_tsv/
ls single/output_tsv/
```

## Resuming after failures

### Pan-Draft stage 2 timeout

When stage 2 times out for an SGB:

```bash
SGB=MGYG000291777
SCRIPT_DIR=/path/to/pan-draft
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

### Pan-Draft stage 4 stalls under pFBA

If `stdout_step4_<SGB>_*.txt` shows no progress for hours and the model is not unusually large (under 1500 reactions), GLPK simplex degeneracy is the likely cause. Switch to FBA mode:

```bash
sbatch --export=ALL,SGB=<SGB_ID>,SCRIPT_DIR=$SCRIPT_DIR run_fba_pan.sh
```

This replaces stage 4 only; stages 1–3 are not re-run.

### Pan-Draft stage 5 stalls under pFBA

The same three models that fail pFBA at stage 4 also fail at stage 5. Force the FBA solver:

```bash
USE_FBA=TRUE bash run_rich_loo_all.sh pan <SGB_ID>
```

This replaces stage 5 only; stages 1–4 are not re-run. The cross-reference logic in `leave_one_out_rich.R` reads the existing minimal-medium `_auxotrophies.tsv`, so re-running stage 4 is unnecessary.

### Single-genome: one array task failed

```bash
# Find which array index failed (e.g. task 3 in job 8612345)
ls -lh $SINGLE_BASE/faa/logs/stderr_8612345_*.txt | sort -k5 -n -r

# Resubmit only that index
sbatch --array=3 \
    --partition=$SLURM_PARTITION \
    --output=$SINGLE_BASE/faa/logs/stdout_%A_%a.txt \
    --error=$SINGLE_BASE/faa/logs/stderr_%A_%a.txt \
    --export=ALL,WORKDIR=$SINGLE_BASE/faa,SCRIPT_DIR=/path/to/pan-draft \
    /path/to/single/single_step.sh
```

## Re-running

### Re-running a pan SGB from scratch

```bash
SGB=MGYG000291777
rm -rf $PIPELINE_BASE/${SGB}_rep/models/ $PIPELINE_BASE/${SGB}_rep/pan_model/

# Remove the SGB's row from cumulative tables
for f in output_tsv/all_reconstruction.tsv output_tsv/all_stage1_summary.tsv \
         output_tsv/qc_battery_results.tsv \
         output_tsv/merged/all_rich_LOO.tsv \
         output_tsv/merged/all_rich_LOO_summary.tsv; do
    grep -v "${SGB}_pan" "$f" > "$f.tmp" && mv "$f.tmp" "$f"
done

bash run_pan_pipeline.sh $SGB
bash run_rich_loo_all.sh pan $SGB     # (or USE_FBA=TRUE for the three unstable models)
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

After re-running, rerun `extract_media.sh` to refresh `all_predicted_media.tsv`.

---

# Methodology notes (consolidated)

## Pan-reactome construction

Reactions present in fewer than 6 % of cluster MAGs (default `--min.rxn.freq.in.mods 0.06`) are excluded from the pan-reactome to mitigate noise from MAG fragmentation. This threshold follows De Bernardini et al. (2024).

## Gap-filling medium

Both pipelines gap-fill against minimal media predicted by `gapseq medium`, not against rich media (e.g. TSB). This avoids retention of biosynthetic shortcuts that would otherwise be falsely detected as auxotrophies. If a project requires comparison against rich-medium reconstructions, modify the `gapseq fill` line in `single_step.sh` (or `pan_step2.sh`) to use `-n /path/to/TSBmed.csv` and document the choice in your Methods section.

## Solver methodology

- **pFBA (default)**: two-stage LP — maximise biomass, then minimise total flux among biomass-optimal solutions. Produces a sparse, near-unique flux distribution physiologically more interpretable than an arbitrary FBA optimum. Used everywhere by default.
- **FBA (fallback)**: single-stage LP — maximise biomass only. Used when pFBA hangs on GLPK simplex degeneracy. Growth values are numerically equivalent to the pFBA stage-1 result; `Total_Flux` column is `NA`.

Both engines use the GLPK solver via the `cobrar` R package. Growth-rate values (h⁻¹) returned by both engines are directly comparable: they are the optimised biomass-reaction flux from the same LP formulation.

## Rich-medium LOO methodology

`leave_one_out_rich.R` performs a parallel leave-one-out scan to distinguish robust auxotrophies from condition-dependent essentials. Every organic exchange is opened to `RICH_LB = −10 mmol gDW⁻¹ h⁻¹` before each compound is dropped, and the resulting growth is classified using the same 10⁻⁴ h⁻¹ threshold. The same eight inorganic exchanges (H₂O, H⁺, phosphate, ammonium, sulfate, K⁺, Mg²⁺, Na⁺) are excluded from the organic set, matching `leave_one_out.R`.

The reverse case (rich-essential but not min-essential) is expected to be zero under FBA monotonicity with respect to exchange bounds and is reported as `X_rich_only` only as a sanity flag — these rows are excluded from the final 4-tier table.

## QC validation note

The QC implementation has been validated by manual inspection of a known balanced reaction (rxn00216, ATP:D-glucose 6-phosphotransferase) and by injection of an artificial imbalance that was correctly detected, confirming the parser and column-sum logic. When all reconstructions originate from `gapseq` and inherit ModelSEED-curated stoichiometry, mass and element balance pass rates of 100 % are expected (Henry et al. 2010; Seaver et al. 2021).

## Known limitations

- LOO occasionally classifies internal metabolic intermediates (e.g. chorismate) as essential due to topological gaps. Manual filtering of the auxotrophy list against ModelSEED metabolite categories is recommended, or use of the Tier A subset (Part II.C) to focus on high-confidence calls.
- `gapseq medium` predictions are minimal by design; auxotrophy counts on TSB-fill models will differ from minimal-fill models. The methodological choice must be consistent within a given comparison.
- FBA-mode phenotyping does not produce parsimonious flux distributions. Comparisons of `Total_Flux` between pFBA-mode and FBA-mode rows are therefore not meaningful and should be avoided.
- The rich-medium LOO does not eliminate gap-fill artefacts: a compound that gap-fill silently bridged appears as Tier B (not in the essential-factor list) regardless of whether the bridging reactions correspond to real biology. The tier framework is a confidence stratifier, not a gap-fill validator. Comparison with cultured reference models (where reconstruction completeness is much higher) remains the most reliable way to flag MAG-completeness artefacts in the pan-Draft set.

---

# Citation

If you use this pipeline, please cite:

- gapseq: Zimmermann et al. (2021), *Genome Biology*, 22:81.
- pan-Draft: De Bernardini et al. (2024), *Genome Biology*, 25(1):280.
- cobrar: Waschina, Zimmermann & Froitzheim (2026), R package v0.2.3, https://github.com/Waschina/cobrar
- MGnify catalogue: Gurbich et al. (2023), *J Mol Biol*, 435(14):168016.
- ModelSEED biochemistry: Henry et al. (2010), *Nat Biotechnol*, 28(9):977-982; Seaver et al. (2021), *Nucleic Acids Res*, 49(D1):D575-D588.
- MEMOTE QC standard: Lieven et al. (2020), *Nat Biotechnol*, 38(3):272-276.
