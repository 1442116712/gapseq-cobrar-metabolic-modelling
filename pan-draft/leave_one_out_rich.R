#!/usr/bin/env Rscript
# leave_one_out_rich.R
#
# Rich-medium LOO companion to leave_one_out.R / leave_one_out_fba.R.
# Distinguishes "hard" auxotrophies (still essential when every organic
# exchange is opened to RICH_LB) from "network-coupled" essentials (only
# essential under the gap-filled minimal medium of the original LOO).
#
# Inputs (env vars, matching the original scripts):
#   MODEL_FILE      -- path to panModel.RDS
#   OUTPUT_PREFIX   -- file-name prefix for outputs
#   ORGANISM_ID     -- organism column value in TSVs
#   NUTRIENTS_TSV   -- (optional) gapseq nutrients.tsv for compound names
#   OUTPUT_DIR      -- output directory (default ".")
#   USE_FBA         -- "TRUE" to use fba() instead of pfba() (default FALSE);
#                      set this for the one model that triggered GLPK simplex
#                      degeneracy in pFBA mode (MGYG000291777)
#   RICH_LB         -- rich-medium lower bound for every organic exchange
#                      (default -10, matching the substrate-scan bound)
#
# Outputs:
#   <OUTPUT_PREFIX>_rich_LOO.tsv          -- per-compound rich-medium LOO
#                                            with tier classification
#   <OUTPUT_PREFIX>_rich_LOO_summary.tsv  -- per-organism summary
#
# Tier classification (rich-medium LOO x minimal-medium LOO):
#   A_hard            essential in BOTH rich and minimal medium
#                     -> robust auxotroph (Tier A)
#   C_network_coupled essential only in minimal medium, not in rich
#                     -> condition-dependent (Tier C)
#   X_rich_only       essential only in rich medium (shouldn't normally
#                     happen under FBA, flagged for inspection)
#   not_essential     never essential

suppressPackageStartupMessages(library(cobrar))

# ---------- Inputs ----------
model_file    <- Sys.getenv("MODEL_FILE")
out_prefix    <- Sys.getenv("OUTPUT_PREFIX")
organism_id   <- Sys.getenv("ORGANISM_ID")
nutrients_tsv <- Sys.getenv("NUTRIENTS_TSV")
output_dir    <- Sys.getenv("OUTPUT_DIR")
use_fba_env   <- toupper(Sys.getenv("USE_FBA", "FALSE"))
use_fba       <- use_fba_env %in% c("TRUE", "T", "1", "YES")
rich_lb_env   <- Sys.getenv("RICH_LB", "-10")
RICH_LB       <- as.numeric(rich_lb_env)

if (model_file == "")  stop("MODEL_FILE not set.")
if (out_prefix == "")  out_prefix  <- sub("\\.RDS$", "", basename(model_file))
if (organism_id == "") organism_id <- out_prefix
if (output_dir == "")  output_dir  <- "."

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
out_path <- function(suffix) file.path(output_dir, paste0(out_prefix, suffix))

mod <- readRDS(model_file)
solve_growth <- if (use_fba) fba else pfba

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
cat(sprintf(">>> Rich-medium LOO -- %s\n", organism_id))
cat(sprintf(">>> Solver = %s, RICH_LB = %.2f mmol gDW^-1 h^-1\n",
            if (use_fba) "FBA" else "pFBA", RICH_LB))
cat("=========================================================\n")

# ---------- Identify exchanges ----------
all_rxns <- mod@react_id
ex_rxns  <- all_rxns[grep("^EX_", all_rxns)]
universal_inorganics <- c(
  "EX_cpd00001_e0", "EX_cpd00067_e0", "EX_cpd00009_e0", "EX_cpd00013_e0",
  "EX_cpd00048_e0", "EX_cpd00205_e0", "EX_cpd00254_e0", "EX_cpd00971_e0"
)
organics <- setdiff(ex_rxns, universal_inorganics)
cat(sprintf(">>> %d organic exchanges to test\n", length(organics)))

# ---------- Build rich medium: open every organic exchange to RICH_LB ----------
mod_rich <- mod
for (rxn in organics) {
  cur <- mod_rich@lowbnd[match(rxn, mod_rich@react_id)]
  if (is.na(cur) || cur > RICH_LB) {
    mod_rich <- changeBounds(mod_rich, react = rxn, lb = RICH_LB)
  }
}

base_growth_rich <- tryCatch(solve_growth(mod_rich)@obj, error = function(e) 0)
if (length(base_growth_rich) == 0 || is.na(base_growth_rich)) base_growth_rich <- 0
cat(sprintf(">>> Base growth in rich medium: %.4f hr^-1\n", base_growth_rich))

if (base_growth_rich < 1e-4) {
  cat("!!! WARNING: rich-medium base growth below 1e-4; LOO results may be unreliable.\n")
}

# ---------- LOO under rich medium ----------
cat(">>> Running rich-medium LOO across all organic exchanges...\n")
loo_growth <- numeric(length(organics))
for (i in seq_along(organics)) {
  rxn <- organics[i]
  mod_drop <- changeBounds(mod_rich, react = rxn, lb = 0)
  tg <- tryCatch(solve_growth(mod_drop)@obj, error = function(e) 0)
  if (length(tg) == 0 || is.na(tg)) tg <- 0
  loo_growth[i] <- as.numeric(tg)
}
rich_essential <- loo_growth < 1e-4
n_rich_ess <- sum(rich_essential)
cat(sprintf("    -> %d compounds essential in rich medium\n", n_rich_ess))

# ---------- Cross-reference with minimal-medium LOO ----------
min_loo_file <- out_path("_auxotrophies.tsv")
min_essentials <- character(0)
if (file.exists(min_loo_file)) {
  min_df <- read.delim(min_loo_file, sep = "\t", stringsAsFactors = FALSE)
  min_essentials <- min_df$reaction_id
  cat(sprintf(">>> Cross-ref with minimal-medium LOO (%s): %d essentials\n",
              basename(min_loo_file), length(min_essentials)))
} else {
  cat(sprintf(">>> No minimal-medium LOO file found at %s; Tier A/C split skipped\n",
              min_loo_file))
}

min_essential <- organics %in% min_essentials

# ---------- Tier classification ----------
tier <- character(length(organics))
tier[ rich_essential &  min_essential] <- "A_hard"
tier[!rich_essential &  min_essential] <- "C_network_coupled"
tier[ rich_essential & !min_essential] <- "X_rich_only"
tier[!rich_essential & !min_essential] <- "not_essential"

cpd_id <- sub("EX_", "", sub("_e0$", "", organics))
compound <- vapply(organics, cpd_lookup, character(1))

loo_records <- data.frame(
  organism       = organism_id,
  reaction_id    = organics,
  cpd_id         = cpd_id,
  compound       = compound,
  growth_without = round(loo_growth, 6),
  rich_essential = rich_essential,
  min_essential  = min_essential,
  tier           = tier,
  stringsAsFactors = FALSE
)
loo_records <- loo_records[order(loo_records$tier, loo_records$cpd_id), ]

write.table(loo_records, out_path("_rich_LOO.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# ---------- Summary ----------
summary_df <- data.frame(
  organism            = organism_id,
  base_growth_rich    = round(base_growth_rich, 4),
  n_organics_tested   = length(organics),
  n_min_essential     = length(min_essentials),
  n_rich_essential    = n_rich_ess,
  n_tier_A_hard       = sum(tier == "A_hard"),
  n_tier_C_network    = sum(tier == "C_network_coupled"),
  n_tier_X_richonly   = sum(tier == "X_rich_only"),
  solver              = if (use_fba) "FBA" else "pFBA",
  rich_lb             = RICH_LB,
  stringsAsFactors    = FALSE
)
write.table(summary_df, out_path("_rich_LOO_summary.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

cat(sprintf("\n>>> Output written to: %s/%s_rich_LOO.tsv\n",
            output_dir, out_prefix))
cat(sprintf(">>> Tier A (hard) = %d  |  Tier C (network-coupled) = %d  |  Tier X = %d\n",
            summary_df$n_tier_A_hard,
            summary_df$n_tier_C_network,
            summary_df$n_tier_X_richonly))
