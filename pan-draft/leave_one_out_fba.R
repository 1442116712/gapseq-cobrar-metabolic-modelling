#!/usr/bin/env Rscript
# leave_one_out_fba.R
#
# A drop-in variant of leave_one_out.R that replaces pfba() with fba() to
# bypass the pFBA stage-2 (parsimony) step that triggers GLPK numerical
# instability on this specific pan-model (MGYG000291777 minimal-medium fill).
#
# Inputs (env vars, identical to leave_one_out.R):
#   MODEL_FILE      -- path to panModel.RDS
#   OUTPUT_PREFIX   -- file-name prefix for outputs
#   ORGANISM_ID     -- organism column value in TSVs
#   NUTRIENTS_TSV   -- (optional) path to gapseq nutrients.tsv for compound names
#   OUTPUT_DIR      -- output directory (default ".")
#
# Output: same 4 TSVs as leave_one_out.R, but Total_Flux column is NA
# (fba() does not minimise total flux).

suppressPackageStartupMessages(library(cobrar))

# ---------- Inputs ----------
model_file    <- Sys.getenv("MODEL_FILE")
out_prefix    <- Sys.getenv("OUTPUT_PREFIX")
organism_id   <- Sys.getenv("ORGANISM_ID")
nutrients_tsv <- Sys.getenv("NUTRIENTS_TSV")
output_dir    <- Sys.getenv("OUTPUT_DIR")

if (model_file == "")  stop("MODEL_FILE not set.")
if (out_prefix == "")  out_prefix  <- sub("\\.RDS$", "", basename(model_file))
if (organism_id == "") organism_id <- out_prefix
if (output_dir == "")  output_dir  <- "."

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
out_path <- function(suffix) file.path(output_dir, paste0(out_prefix, suffix))

mod <- readRDS(model_file)

# Compound-name lookup
cpd_lookup <- function(rxn_id) NA_character_
if (nzchar(nutrients_tsv) && file.exists(nutrients_tsv)) {
  nutr <- read.delim(nutrients_tsv, sep = "\t", stringsAsFactors = FALSE)
  cpd_lookup <- function(rxn_id) {
    cpd <- sub("EX_", "", sub("_e0$", "", rxn_id))
    nm  <- nutr$name[match(cpd, nutr$id)]
    if (length(nm) == 0 || is.na(nm)) NA_character_ else nm
  }
}

cat("\n=========================================================\n")
cat(sprintf(">>> Phenotyping (FBA mode, no parsimony) -- %s\n", organism_id))
cat("=========================================================\n")

# ---------- Identify exchanges ----------
all_rxns <- mod@react_id
ex_rxns  <- all_rxns[grep("^EX_", all_rxns)]

universal_inorganics <- c(
  "EX_cpd00001_e0", "EX_cpd00067_e0", "EX_cpd00009_e0", "EX_cpd00013_e0",
  "EX_cpd00048_e0", "EX_cpd00205_e0", "EX_cpd00254_e0", "EX_cpd00971_e0"
)
organics <- setdiff(ex_rxns, universal_inorganics)

# ---------- Module A: LOO auxotrophy ----------
cat(">>> [1/3] LOO auxotrophy detection (FBA)...\n")
base_growth <- fba(mod)@obj
dynamic_auxotrophies <- c()
open_organics <- organics[mod@lowbnd[match(organics, mod@react_id)] < 0]

for (rxn in open_organics) {
  mod_dropout <- changeBounds(mod, react = rxn, lb = 0)
  tg <- tryCatch(fba(mod_dropout)@obj, error = function(e) 0)
  if (length(tg) == 0 || is.na(tg)) tg <- 0
  if (tg < 1e-4) dynamic_auxotrophies <- c(dynamic_auxotrophies, rxn)
}
cat(sprintf("    -> %d essential factors found.\n", length(dynamic_auxotrophies)))

aux_df <- data.frame(
  organism    = organism_id,
  reaction_id = dynamic_auxotrophies,
  cpd_id      = sub("EX_", "", sub("_e0$", "", dynamic_auxotrophies)),
  compound    = vapply(dynamic_auxotrophies, cpd_lookup, character(1)),
  stringsAsFactors = FALSE
)
write.table(aux_df, out_path("_auxotrophies.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# ---------- Module B: trace medium ----------
cat(">>> [2/3] Trace-limited background (FBA)...\n")
mod_trace <- mod
for (rxn in organics) {
  if (mod_trace@lowbnd[match(rxn, mod_trace@react_id)] < -0.05) {
    mod_trace <- changeBounds(mod_trace, react = rxn, lb = -0.05)
  }
}
if (length(dynamic_auxotrophies) > 0) {
  mod_trace <- changeBounds(mod_trace, react = dynamic_auxotrophies, lb = -0.1)
}
bg_growth <- tryCatch(fba(mod_trace)@obj, error = function(e) 0)
cat(sprintf("    -> Trace bg growth: %.4f hr^-1\n", bg_growth))

# ---------- Module C: substrate scan ----------
cat(">>> [3/3] Substrate scan (FBA, no parsimony)...\n")
candidate_sources <- setdiff(organics, dynamic_auxotrophies)

results_df <- data.frame(
  Reaction_ID = character(),
  Growth_Rate = numeric(),
  Total_Flux  = numeric(),
  stringsAsFactors = FALSE
)

for (rxn_id in candidate_sources) {
  mod_test <- changeBounds(mod_trace, react = rxn_id, lb = -10)
  test_res <- tryCatch(
    suppressWarnings(fba(mod_test)),
    error = function(e) NULL
  )
  if (is.null(test_res)) {
    gr <- 0
  } else {
    gr <- if (.hasSlot(test_res, "obj")) test_res@obj else 0
    if (length(gr) == 0 || is.na(gr)) gr <- 0
  }
  # fba() does not provide total flux; leave NA for cross-script consistency
  results_df <- rbind(results_df, data.frame(
    Reaction_ID = rxn_id,
    Growth_Rate = round(as.numeric(gr), 6),
    Total_Flux  = NA_real_
  ))
}

results_df$Net_Growth <- round(results_df$Growth_Rate - bg_growth, 6)
results_df$organism <- organism_id
results_df$cpd_id   <- sub("EX_", "", sub("_e0$", "", results_df$Reaction_ID))
results_df$compound <- vapply(results_df$Reaction_ID, cpd_lookup, character(1))
results_df <- results_df[, c("organism", "Reaction_ID", "cpd_id", "compound",
                             "Growth_Rate", "Total_Flux", "Net_Growth")]

write.table(results_df[order(-results_df$Net_Growth), ],
            out_path("_substrates_full.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

active <- results_df[results_df$Net_Growth > 0.01, ]
active <- active[order(-active$Net_Growth), ]
write.table(head(active, 15),
            out_path("_substrates_top15.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

summary_df <- data.frame(
  organism            = organism_id,
  n_reactions         = length(all_rxns),
  n_exchanges         = length(ex_rxns),
  n_organics          = length(organics),
  n_auxotrophies      = length(dynamic_auxotrophies),
  base_growth         = round(base_growth, 4),
  trace_bg_growth     = round(bg_growth, 4),
  n_active_substrates = nrow(active)
)
write.table(summary_df, out_path("_summary.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

cat(sprintf("\n>>> Output written to: %s/%s_*.tsv\n",
            output_dir, out_prefix))
cat(sprintf(">>> Note: Total_Flux column is NA (fba mode, no parsimony).\n"))
