#!/usr/bin/env Rscript
# fab_rescue_test.R
#
# In-silico pre-experiment for the FAB rescue assay (Manuscript 4, Exp 2).
#
# Workflow (all on the DSM 11370 / GCF_000426565.1 single-genome model):
#   [1] Build a FAB-only medium from FAB.csv (close all organic exchanges,
#       reopen those listed in FAB.csv to lb = -maxFlux).
#       Expectation: growth ~ 0 (mirrors wet-lab observation).
#   [2] Diagnose: cross-reference DSM 11370's 10 published auxotrophies
#       (S11b) against FAB.csv. List per-compound availability and the
#       wet-lab SOP cocktail status.
#   [3] Add the 7-component SOP cocktail at SOP-derived flux (max(FAB,
#       cocktail) per exchange). Expectation: growth recovers.
#   [4] Leave-one-out within the cocktail: drop each component back to
#       FAB-baseline, re-test growth. Identifies which cocktail
#       components are individually necessary for rescue.
#   [5] Saturating-cocktail sanity check (all 7 at lb = -10): confirms
#       the model CAN grow given ample supply — i.e. rules out an
#       unrelated auxotrophy.
#   [6] Positive control on the gapseq-default minimal medium (if a
#       *-medium.csv is colocated with the RDS).
#
# Inputs (env vars):
#   MODEL_FILE      -- /path/to/GCF_000426565.1.RDS
#   FAB_CSV         -- /path/to/FAB.csv
#   OUTPUT_DIR      -- output directory (default ".")
#   USE_FBA         -- "TRUE" to use fba() instead of pfba() (default TRUE,
#                      because pFBA can hang on simplex degeneracy and we only
#                      care about growth/no-growth, not flux distribution)
#
# Outputs (in OUTPUT_DIR):
#   DSM11370_fab_rescue_summary.tsv  -- one row per condition
#   DSM11370_fab_diagnostics.tsv     -- per-auxotrophy availability + SOP

suppressPackageStartupMessages(library(cobrar))

# ---------------------- Configuration ----------------------
model_file <- Sys.getenv("MODEL_FILE",
  "/users/40335635/sharedscratch/EggNOG/manuscript_4/single_genome/faa/GCF_000426565.1.RDS")
fab_csv    <- Sys.getenv("FAB_CSV",
  "/users/40335635/sharedscratch/EggNOG/manuscript_4/single_genome/FAB.csv")
output_dir <- Sys.getenv("OUTPUT_DIR", ".")
# Default to FBA (faster, no simplex-degeneracy risk; we only need growth y/n).
use_fba    <- toupper(Sys.getenv("USE_FBA", "TRUE")) %in% c("TRUE","T","1","YES")

# Timing helper -- prints elapsed wall-clock per step for diagnostics
t0 <- Sys.time()
tic <- function(msg) {
  cat(sprintf("[%6.1fs] %s\n",
              as.numeric(Sys.time() - t0, units = "secs"), msg))
  flush.console()
}

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
out_path <- function(suffix) file.path(output_dir, paste0("DSM11370_", suffix))
solve_growth <- if (use_fba) fba else pfba

cat("=========================================================\n")
cat(sprintf(">>> FAB rescue in-silico pre-experiment\n"))
cat(sprintf(">>> Solver: %s\n", if (use_fba) "FBA" else "pFBA"))
cat(sprintf(">>> Model:  %s\n", model_file))
cat(sprintf(">>> FAB:    %s\n", fab_csv))
cat(sprintf(">>> Out:    %s\n", output_dir))
cat("=========================================================\n")

if (!file.exists(model_file)) stop("MODEL_FILE not found.")
if (!file.exists(fab_csv))    stop("FAB_CSV not found.")

tic("Loading model RDS")
mod <- readRDS(model_file)
tic(sprintf("Loaded: %d reactions, %d metabolites",
            length(mod@react_id), length(mod@met_id)))

# ---------------------- Helpers ----------------------
ex_of <- function(cpd) sprintf("EX_%s_e0", cpd)
has_ex <- function(m, cpd) ex_of(cpd) %in% m@react_id
get_lb <- function(m, ex) m@lowbnd[match(ex, m@react_id)]

growth_of <- function(m) {
  s <- tryCatch(solve_growth(m)@obj, error = function(e) 0)
  if (length(s) == 0 || is.na(s)) 0 else as.numeric(s)
}

# Batched medium setter: one changeBounds call total, not one per compound.
# Looping changeBounds N times can be ~N x slower than a single vectorised call
# because cobrar re-validates the model after each modification.
set_medium <- function(m, medium_df, close_first = TRUE) {
  all_ex <- m@react_id[grep("^EX_", m@react_id)]
  # Build target lb per exchange: 0 by default (closed), negative for medium
  target_lb <- setNames(rep(0, length(all_ex)), all_ex)
  ex_ids <- ex_of(medium_df$compounds)
  flx    <- as.numeric(medium_df$maxFlux)
  keep   <- ex_ids %in% all_ex
  target_lb[ex_ids[keep]] <- -flx[keep]
  if (close_first) {
    m <- changeBounds(m, react = all_ex, lb = unname(target_lb))
  } else {
    # Only modify exchanges that are in the medium; leave others untouched.
    m <- changeBounds(m, react = ex_ids[keep], lb = -flx[keep])
  }
  m
}

# ---------------------- 1. FAB-only ----------------------
fab <- read.csv(fab_csv, stringsAsFactors = FALSE)
tic(sprintf("[Step 1] Building FAB-only medium (%d compounds)", nrow(fab)))
mod_fab <- set_medium(mod, fab, close_first = TRUE)
tic("[Step 1] Medium set, solving baseline growth")
g_fab <- growth_of(mod_fab)
tic(sprintf("[Step 1] FAB-only growth = %.4f h^-1 %s", g_fab,
            if (g_fab < 1e-4) "[no growth, matches wet-lab]"
            else              "[WARN: model grows, expectation was no growth]"))

# ---------------------- 2. Diagnose auxotrophies vs FAB ----------------------
# DSM 11370's 10 published auxotrophies (S11b), with SOP-cocktail status.
# SOP final concentrations from Invitro_supplement_protocol.docx.
auxotrophies <- data.frame(
  cpd_id    = c("cpd00030","cpd00034","cpd00149","cpd10515","cpd10516",
                "cpd00063","cpd00058","cpd00099","cpd00239","cpd00305"),
  name      = c("Mn2+","Zn2+","Co2+","Fe2+","Fe3+",
                "Ca2+","Cu2+","Cl-","H2S","Thiamin"),
  sop_added = c(TRUE,TRUE,TRUE,TRUE,FALSE,
                TRUE,TRUE,FALSE,FALSE,TRUE),
  sop_uM    = c(10,10,10,25,0,
                500,2,0,0,1),
  sop_note  = c("","","","FeSO4 only", "not added (FeSO4 -> Fe2+ only)",
                "","","FAB matrix (NaCl)","FAB matrix (thioglyc/Cys)",""),
  stringsAsFactors = FALSE
)
auxotrophies$in_fab <- vapply(auxotrophies$cpd_id,
  function(c) c %in% fab$compounds, logical(1))
auxotrophies$fab_maxFlux <- vapply(auxotrophies$cpd_id, function(c) {
  r <- fab[fab$compounds == c, , drop = FALSE]
  if (nrow(r) > 0) as.numeric(r$maxFlux[1]) else 0
}, numeric(1))
auxotrophies$model_has_ex <- vapply(auxotrophies$cpd_id,
  function(c) has_ex(mod, c), logical(1))

tic("[Step 2] Auxotrophy availability in FAB:")
print(auxotrophies[, c("cpd_id","name","in_fab","fab_maxFlux",
                       "sop_added","sop_uM")], row.names = FALSE)
write.table(auxotrophies, out_path("fab_diagnostics.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# ---------------------- 3. FAB + SOP cocktail ----------------------
# Cocktail entries (SOP). flux = 2 x mM (assumes OD600 ~1, gDW/L ~0.5,
# 1 h consumption convention).
cocktail <- data.frame(
  cpd_id = c("cpd00030","cpd00034","cpd00149","cpd00063",
             "cpd00058","cpd10515","cpd00305"),
  name   = c("Mn2+","Zn2+","Co2+","Ca2+","Cu2+","Fe2+","Thiamin"),
  uM     = c(10,10,10,500,2,25,1),
  stringsAsFactors = FALSE
)
cocktail$flux <- 2 * cocktail$uM / 1000   # uM -> mM -> mmol/gDW/h

tic("[Step 3] Adding SOP cocktail (batched)")
# Compute new lb for cocktail exchanges = max(|current|, cocktail flux); apply in one call
cocktail_ex <- ex_of(cocktail$cpd_id)
in_model <- cocktail_ex %in% mod_fab@react_id
cur_lb   <- get_lb(mod_fab, cocktail_ex)
new_lb   <- -pmax(abs(cur_lb), cocktail$flux)
mod_fc <- changeBounds(mod_fab,
                       react = cocktail_ex[in_model],
                       lb    = new_lb[in_model])
tic("[Step 3] Solving FAB+cocktail growth")
g_fc <- growth_of(mod_fc)
tic(sprintf("[Step 3] FAB + SOP cocktail growth = %.4f h^-1 %s", g_fc,
            if (g_fc > 1e-4) "[rescued]" else "[WARN: still no growth]"))

# ---------------------- 4. Cocktail leave-one-out ----------------------
tic("[Step 4] Cocktail leave-one-out (one solve per component)")
loo_rows <- list()
for (i in seq_len(nrow(cocktail))) {
  step_start <- Sys.time()
  ex <- ex_of(cocktail$cpd_id[i])
  if (ex %in% mod_fc@react_id) {
    fab_flux <- auxotrophies$fab_maxFlux[
      match(cocktail$cpd_id[i], auxotrophies$cpd_id)]
    if (is.na(fab_flux)) fab_flux <- 0
    mod_loo <- changeBounds(mod_fc, react = ex, lb = -fab_flux)
    g_loo <- growth_of(mod_loo)
  } else {
    g_loo <- NA_real_
  }
  dt <- as.numeric(Sys.time() - step_start, units = "secs")
  loo_rows[[i]] <- data.frame(
    condition   = sprintf("FAB+cocktail (drop %s)", cocktail$name[i]),
    growth_rate = g_loo,
    stringsAsFactors = FALSE
  )
  tic(sprintf("    drop %-10s  growth = %.4f  (%.1fs)  %s",
              cocktail$name[i], g_loo, dt,
              if (!is.na(g_loo) && g_loo < 1e-4) "[essential to rescue]"
              else                                "[dispensable]"))
}

# ---------------------- 5. Saturating cocktail sanity check ----------------------
tic("[Step 5] Saturating cocktail (lb=-10), batched")
mod_sat <- changeBounds(mod_fab,
                        react = cocktail_ex[in_model],
                        lb    = rep(-10, sum(in_model)))
g_sat <- growth_of(mod_sat)
tic(sprintf("[Step 5] Saturating cocktail growth = %.4f h^-1", g_sat))

# ---------------------- 6. Positive control on gapseq min medium ----------------------
min_csv <- sub("\\.RDS$", "-medium.csv", model_file)
g_min <- NA_real_
if (file.exists(min_csv)) {
  tic("[Step 6] Positive control on gapseq min medium")
  min_df <- read.csv(min_csv, stringsAsFactors = FALSE)
  mod_min <- set_medium(mod, min_df, close_first = TRUE)
  g_min <- growth_of(mod_min)
  tic(sprintf("[Step 6] gapseq min medium growth = %.4f h^-1", g_min))
} else {
  tic(sprintf("[Step 6] No %s — skipping positive control", basename(min_csv)))
}

# ---------------------- 7. Cobalamin-family sweep (FAA / blood rescue hypothesis) ----------------------
# Wet-lab observation: DSM 11370 does NOT grow on FAB but DOES grow on FAA
# (FAB + 5% horse blood). Hypothesis: blood provides cobalamin (some form
# of vitamin B12) directly, bypassing the de novo B12 synthesis pathway
# that requires Co2+. Since the model may import cobalamin in any of several
# pathway intermediates (CN-Cbl, AdoCbl, cobinamide, etc.), we auto-discover
# which cobalamin-family exchanges exist and sweep each.
#
# Approximate concentration mapping (flux 2 = 1 mM convention):
#   flux 1e-6 -> ~0.5 nM   flux 1e-4 -> ~50 nM (physiological blood B12)
#   flux 1e-2 -> ~5 uM     flux 1    -> ~500 uM (saturating)
tic("[Step 7] Cobalamin-family exchange scan (with model self-inspection)")

# Pattern for B12 family names (excludes "Cobalt" metal)
pat <- "(?i)(cobalamin|cobinamide|cobyrate|cobyric|cob\\(i+\\)|adenosylcob|corrinoid)"

# ---- Approach A: model self-inspection ----
# Dump ALL exchange reactions to file so you can inspect what the model
# actually exposes. Use the model's own metabolite naming if available.
all_ex_rxns <- mod_fab@react_id[grep("^EX_", mod_fab@react_id)]
all_ex_cpds <- sub("^EX_", "", sub("_e0$", "", all_ex_rxns))

# Try to get metabolite names from the model itself (multiple possible slots)
get_met_names <- function(m, cpd_ids) {
  # Try met_attr$name, met_name, or fall back to NA
  if ("met_attr" %in% slotNames(m) && "name" %in% names(m@met_attr)) {
    mid <- paste0(cpd_ids, "_e0")
    return(m@met_attr$name[match(mid, m@met_id)])
  }
  if ("met_name" %in% slotNames(m)) {
    mid <- paste0(cpd_ids, "_e0")
    return(m@met_name[match(mid, m@met_id)])
  }
  rep(NA_character_, length(cpd_ids))
}
ex_names_model <- get_met_names(mod_fab, all_ex_cpds)

# Also resolve names via nutrients.tsv (gapseq dict)
nutrients_tsv <- Sys.getenv("NUTRIENTS_TSV", "")
ex_names_nutr <- rep(NA_character_, length(all_ex_cpds))
if (nzchar(nutrients_tsv) && file.exists(nutrients_tsv)) {
  nutr <- read.delim(nutrients_tsv, sep = "\t", stringsAsFactors = FALSE)
  ex_names_nutr <- nutr$name[match(all_ex_cpds, nutr$id)]
}

# Combine: prefer model's own name, fall back to nutrients.tsv
ex_names <- ifelse(!is.na(ex_names_model) & ex_names_model != "",
                   ex_names_model, ex_names_nutr)

all_ex_df <- data.frame(
  reaction_id = all_ex_rxns,
  cpd_id      = all_ex_cpds,
  name        = ex_names,
  lb          = mod_fab@lowbnd[match(all_ex_rxns, mod_fab@react_id)],
  ub          = mod_fab@uppbnd[match(all_ex_rxns, mod_fab@react_id)],
  stringsAsFactors = FALSE
)
write.table(all_ex_df, out_path("model_all_exchanges.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
tic(sprintf("[Step 7] Model has %d exchange reactions in total -> %s",
            nrow(all_ex_df), out_path("model_all_exchanges.tsv")))

# ---- Approach B: search by name in the dumped exchange list ----
hits_by_name <- all_ex_df[!is.na(all_ex_df$name) &
                          grepl(pat, all_ex_df$name, perl = TRUE), ]
if (nrow(hits_by_name) > 0) {
  tic(sprintf("[Step 7] By NAME search in model exchanges: %d cobalamin-family hits",
              nrow(hits_by_name)))
  print(hits_by_name[, c("reaction_id","cpd_id","name","lb")], row.names = FALSE)
} else {
  tic("[Step 7] By NAME search in model exchanges: 0 cobalamin-family hits")
}

# ---- Approach C: search by curated cpd IDs (sanity backup) ----
cobl_db <- c(
  cpd00166 = "Adenosylcobalamin (gapseq canonical B12)",
  cpd00423 = "Cob(II)alamin (B12r)",
  cpd00635 = "Cob(I)alamin (B12s)",
  cpd03422 = "Cobinamide (B12 precursor)",
  cpd03424 = "Cobalamin (legacy ID)"
)
hits_by_id <- cobl_db[names(cobl_db) %in% all_ex_df$cpd_id]
if (length(hits_by_id) > 0) {
  tic(sprintf("[Step 7] By curated cpd ID search: %d hits  (%s)",
              length(hits_by_id), paste(names(hits_by_id), collapse = ", ")))
}

# ---- Also search ALL reactions (not just exchanges) for cobalamin keywords ----
# This tells us whether the model HAS B12 synthesis/use reactions even without
# a transporter.
all_rxn_count <- length(mod_fab@react_id)
# Search reaction names if available (try multiple possible slot locations)
rxn_name_slot <- rep(NA_character_, all_rxn_count)
if ("react_name" %in% slotNames(mod_fab)) {
  rxn_name_slot <- mod_fab@react_name
} else if ("react_attr" %in% slotNames(mod_fab) &&
           "name" %in% names(mod_fab@react_attr)) {
  rxn_name_slot <- mod_fab@react_attr$name
}
b12_rxn_hits <- mod_fab@react_id[!is.na(rxn_name_slot) &
                                 grepl(pat, rxn_name_slot, perl = TRUE)]
if (length(b12_rxn_hits) > 0) {
  tic(sprintf("[Step 7] %d reactions in model match B12 keyword (transport/synthesis/use):",
              length(b12_rxn_hits)))
  for (rxn in head(b12_rxn_hits, 20)) {
    nm <- rxn_name_slot[match(rxn, mod_fab@react_id)]
    cat(sprintf("       %s : %s\n", rxn, nm))
  }
  if (length(b12_rxn_hits) > 20) cat(sprintf("       ... and %d more\n",
                                              length(b12_rxn_hits) - 20))
}

# Use the union of name-based and ID-based hits as the set to test
cobl_present <- setNames(
  c(hits_by_name$name, hits_by_id),
  c(hits_by_name$cpd_id, names(hits_by_id))
)
cobl_present <- cobl_present[!duplicated(names(cobl_present))]
tic(sprintf("[Step 7] -> %d UNIQUE cobalamin-family exchanges to sweep",
            length(cobl_present)))

# Also report what's in the gapseq min medium (the model grows on min medium,
# so this shows where the model's B12/cofactor pool comes from in the positive
# control)
if (exists("min_df")) {
  # Check both curated IDs and nutrients.tsv name-pattern matches
  cobl_candidates <- unique(c(names(cobl_db),
                              if (exists("nutr")) {
                                nutr$id[grep(pat, nutr$name, perl = TRUE)]
                              } else character(0)))
  cobl_in_min <- intersect(min_df$compounds, cobl_candidates)
  if (length(cobl_in_min) > 0) {
    tic(sprintf("    gapseq min medium contains: %s",
                paste(cobl_in_min, collapse = ", ")))
  } else {
    tic("    gapseq min medium has NO cobalamin compound — model must synthesize B12 de novo from Co2+")
  }
}

b12_rows <- list()
if (length(cobl_present) == 0) {
  tic("[Step 7] No cobalamin exchanges in model — cannot directly test blood/B12 rescue.")
  tic("         Interpretation: gapseq did not annotate a B12 transporter for this organism.")
  tic("         The wet-lab FAA rescue (via blood B12) cannot be reproduced in silico without")
  tic("         manually adding a transporter. The Co2+ auxotrophy nevertheless captures the")
  tic("         underlying biology: blood->B12 supply bypasses the same de novo synthesis step")
  tic("         that requires Co2+ in the model.")
} else {
  b12_levels <- c(1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1, 10)
  for (cpd in names(cobl_present)) {
    ex <- ex_of(cpd)
    cat(sprintf("    -- %s (%s) --\n", cpd, cobl_present[[cpd]]))
    for (lvl in b12_levels) {
      mod_b12 <- changeBounds(mod_fab, react = ex, lb = -lvl)
      g <- growth_of(mod_b12)
      b12_rows[[length(b12_rows) + 1]] <- data.frame(
        cpd_id      = cpd,
        compound    = cobl_present[[cpd]],
        flux        = lvl,
        uM_approx   = lvl * 500,
        growth_rate = g,
        stringsAsFactors = FALSE
      )
      tic(sprintf("       flux %.0e (~%.3g uM)  growth = %.4f  %s",
                  lvl, lvl * 500, g,
                  if (g > 1e-4) "[rescued]" else ""))
    }
  }
}
if (length(b12_rows) > 0) {
  b12_results <- do.call(rbind, b12_rows)
  write.table(b12_results, out_path("fab_b12_sweep.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)
}

# ---------------------- 8. Co2+ titration on FAB-only (FAA hypothesis) ----------------------
# Question: is the trace cobalt in 5% horse blood (~0.5-1 nM Co, mostly bound in
# B12) by itself sufficient to rescue growth in silico?
# If the minimum rescue flux is <= ~1e-6 (corresponding to blood Co level), FAA
# rescue can be explained by direct Co2+ supply — no B12 transporter needed.
# If the minimum rescue flux is >> 1e-6, blood Co alone is too dilute; FAA
# rescue must involve B12 import (which this model lacks).
tic("[Step 8] Co2+ titration on FAB-only (testing blood-Co vs blood-B12 mechanism)")
co_ex <- ex_of("cpd00149")
co_rows <- list()
if (!(co_ex %in% mod_fab@react_id)) {
  tic("[Step 8] Model has no Co2+ exchange — should not happen, skipping")
} else {
  # Concentration mapping (flux 2 = 1 mM convention):
  #   1e-8  ~ 5 pM       (sub-physiological)
  #   1e-7  ~ 50 pM      (sub-physiological)
  #   1e-6  ~ 0.5 nM     (matches 5% horse blood Co content)
  #   1e-5  ~ 5 nM
  #   1e-4  ~ 50 nM
  #   1e-3  ~ 0.5 uM
  #   1e-2  ~ 5 uM       (below SOP)
  #   1e-1  ~ 50 uM
  #   1     ~ 500 uM     (saturating)
  co_levels <- c(0, 1e-8, 1e-7, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1, 10)
  for (lvl in co_levels) {
    step_start <- Sys.time()
    mod_co <- changeBounds(mod_fab, react = co_ex, lb = -lvl)
    g_co <- growth_of(mod_co)
    dt <- as.numeric(Sys.time() - step_start, units = "secs")
    co_rows[[length(co_rows) + 1]] <- data.frame(
      co_flux     = lvl,
      co_nM_approx = lvl * 500e3,   # flux 1 = 0.5 mM = 500 000 nM
      growth_rate = g_co,
      rescued     = g_co > 1e-4,
      stringsAsFactors = FALSE
    )
    tic(sprintf("       Co2+ flux = %.0e  (~%.3g nM)  growth = %.4f  (%.1fs)  %s",
                lvl, lvl * 500e3, g_co, dt,
                if (g_co > 1e-4) "[rescued]" else "[no growth]"))
  }
  co_results <- do.call(rbind, co_rows)
  # Find minimum rescue flux
  rescued <- co_results[co_results$rescued, ]
  if (nrow(rescued) > 0) {
    min_rescue_flux <- min(rescued$co_flux)
    min_rescue_nM   <- min_rescue_flux * 500e3
    blood_nM_max    <- 1
    tic(sprintf("[Step 8] Minimum rescue flux = %.0e (~%.3g nM Co2+)",
                min_rescue_flux, min_rescue_nM))
    if (min_rescue_nM <= blood_nM_max) {
      tic("[Step 8] CONCLUSION: blood Co2+ alone is sufficient to rescue.")
      tic("         The FAA rescue can be explained by direct Co2+ supply from blood,")
      tic("         without invoking an external B12 transporter.")
    } else {
      tic(sprintf("[Step 8] CONCLUSION: blood Co2+ (~1 nM) is %.1fx below the in-silico threshold.",
                  min_rescue_nM / blood_nM_max))
      tic("         FAA rescue cannot be explained by Co2+ alone — must involve B12 import")
      tic("         (which this model lacks an exchange for).")
    }
  } else {
    tic("[Step 8] No Co2+ flux tested rescued growth — unexpected, please inspect.")
  }
  write.table(co_results, out_path("fab_co_titration.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)
}

# ---------------------- Write summary ----------------------
summary_df <- rbind(
  data.frame(condition = "gapseq min medium (positive control)",
             growth_rate = g_min, stringsAsFactors = FALSE),
  data.frame(condition = "FAB only (negative expectation)",
             growth_rate = g_fab, stringsAsFactors = FALSE),
  data.frame(condition = "FAB + full SOP cocktail",
             growth_rate = g_fc,  stringsAsFactors = FALSE),
  data.frame(condition = "FAB + saturating cocktail (lb=-10)",
             growth_rate = g_sat, stringsAsFactors = FALSE),
  do.call(rbind, loo_rows)
)
# Append Co2+ titration rows to summary
if (length(co_rows) > 0) {
  co_for_summary <- do.call(rbind, lapply(co_rows, function(r) {
    data.frame(
      condition   = sprintf("FAB + Co2+ only (flux=%.0e, ~%.3g nM)",
                            r$co_flux, r$co_nM_approx),
      growth_rate = r$growth_rate,
      stringsAsFactors = FALSE
    )
  }))
  summary_df <- rbind(summary_df, co_for_summary)
}
# Append B12 sweep rows to summary
if (length(b12_rows) > 0) {
  b12_for_summary <- do.call(rbind, lapply(b12_rows, function(r) {
    data.frame(
      condition   = sprintf("FAB + B12 (flux=%.0e, ~%.3g uM)",
                            r$b12_flux, r$b12_uM_approx),
      growth_rate = r$growth_rate,
      stringsAsFactors = FALSE
    )
  }))
  summary_df <- rbind(summary_df, b12_for_summary)
}
summary_df$grew <- !is.na(summary_df$growth_rate) &
                   summary_df$growth_rate > 1e-4
write.table(summary_df, out_path("fab_rescue_summary.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

cat("\n=========================================================\n")
cat(">>> Summary\n")
cat("=========================================================\n")
print(summary_df, row.names = FALSE)
cat(sprintf("\n>>> Wrote: %s\n", out_path("fab_rescue_summary.tsv")))
cat(sprintf(">>> Wrote: %s\n", out_path("fab_diagnostics.tsv")))
if (length(b12_rows) > 0) {
  cat(sprintf(">>> Wrote: %s\n", out_path("fab_b12_sweep.tsv")))
}
if (length(co_rows) > 0) {
  cat(sprintf(">>> Wrote: %s\n", out_path("fab_co_titration.tsv")))
}
