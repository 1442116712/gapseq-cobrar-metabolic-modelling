#!/usr/bin/env Rscript
# qc_battery.R (v2 - if/else syntax fixed)
# Batch QC battery for pan-Draft metabolic models.
# Four checks:
#   (1) ATP-from-water       : TIC detection
#   (2) Loopless contribution: PLACEHOLDER (cobrar 0.2.x lacks loopless FBA)
#   (3) Total mass balance
#   (4) Per-element balance (C/H/N/O/P/S)
#
# Usage:
#   Rscript qc_battery.R [<base_dir>] [<out_tsv>]

suppressPackageStartupMessages({
  library(cobrar)
  library(Matrix)
})

ATOMIC_WEIGHTS <- c(
  H = 1.008,  C = 12.011, N = 14.007, O = 15.999, P = 30.974, S = 32.06,
  Na = 22.99, K = 39.098, Mg = 24.305, Ca = 40.078, Cl = 35.45,
  Fe = 55.845, Mn = 54.938, Zn = 65.38, Cu = 63.546, Co = 58.933,
  Ni = 58.693, Mo = 95.95, Se = 78.971, F = 18.998, Br = 79.904,
  I = 126.904, B = 10.811, As = 74.922
)
KEY_ELEMENTS <- c("C", "H", "N", "O", "P", "S")

parse_formula_to_counts <- function(f) {
  if (is.null(f) || is.na(f) || f == "" || f == "null" || f == "*" || f == "R") {
    return(NULL)
  }
  matches <- regmatches(f, gregexpr("([A-Z][a-z]?)([0-9]*)", f))[[1]]
  matches <- matches[matches != ""]
  if (length(matches) == 0) return(NULL)
  counts <- list()
  for (m in matches) {
    el      <- gsub("[0-9]", "", m)
    cnt_str <- gsub("[A-Za-z]", "", m)
    cnt     <- if (cnt_str == "") 1 else suppressWarnings(as.numeric(cnt_str))
    if (is.na(cnt)) next
    counts[[el]] <- (if (is.null(counts[[el]])) 0 else counts[[el]]) + cnt
  }
  if (length(counts) == 0) NULL else unlist(counts)
}

build_element_matrix <- function(formulas, element_set) {
  n <- length(formulas)
  E <- matrix(0, nrow = n, ncol = length(element_set),
              dimnames = list(NULL, element_set))
  has_formula <- logical(n)
  for (i in seq_len(n)) {
    counts <- parse_formula_to_counts(formulas[i])
    if (is.null(counts)) next
    has_formula[i] <- TRUE
    for (el in names(counts)) {
      if (el %in% element_set) E[i, el] <- counts[[el]]
    }
  }
  list(E = E, has_formula = has_formula)
}

internal_mask <- function(mod) {
  rxn_ids <- mod@react_id
  is_ex   <- grepl("^EX_", rxn_ids)
  is_bio  <- grepl("biomass|BIOMASS", rxn_ids, ignore.case = TRUE) |
             grepl("^bio[0-9]*$",     rxn_ids)
  is_obj  <- if (length(mod@obj_coef) == length(rxn_ids)) mod@obj_coef != 0
             else rep(FALSE, length(rxn_ids))
  is_dm   <- grepl("^DM_|^SK_|^sink_", rxn_ids)
  !(is_ex | is_bio | is_obj | is_dm)
}

classify_lt <- function(value, pass_th, warn_th) {
  if (is.na(value)) {
    "NA"
  } else if (value < pass_th) {
    "PASS"
  } else if (value < warn_th) {
    "WARN"
  } else {
    "FAIL"
  }
}

classify_gt <- function(value, pass_th, warn_th) {
  if (is.na(value)) {
    "NA"
  } else if (value >= pass_th) {
    "PASS"
  } else if (value >= warn_th) {
    "WARN"
  } else {
    "FAIL"
  }
}

check_atp_from_water <- function(mod) {
  ex_rxns <- mod@react_id[grep("^EX_", mod@react_id)]
  m <- changeBounds(mod, react = ex_rxns, lb = 0)
  for (cpd in c("EX_cpd00001_e0", "EX_cpd00067_e0")) {
    if (cpd %in% m@react_id) m <- changeBounds(m, react = cpd, lb = -1000)
  }
  atpm_idx <- which(m@react_id %in% c("rxn00062", "rxn00062_c0", "ATPM"))
  if (length(atpm_idx) == 0) atpm_idx <- grep("^rxn00062", m@react_id)
  if (length(atpm_idx) == 0) {
    return(list(value = NA_real_, status = "ATPM_not_found"))
  }
  atpm_id <- m@react_id[atpm_idx[1]]
  m@obj_coef <- rep(0, length(m@react_id))
  m@obj_coef[atpm_idx[1]] <- 1
  sol <- tryCatch(fba(m)@obj, error = function(e) NA_real_)
  if (length(sol) == 0 || is.na(sol)) {
    return(list(value = NA_real_, status = "solver_error"))
  }
  list(value = as.numeric(sol), status = paste0("ok:", atpm_id))
}

check_loopless <- function(mod) {
  list(base = NA_real_, loopless = NA_real_, penalty_pct = NA_real_,
       status = "cobrar_no_loopless")
}

check_total_mass_balance <- function(mod, elem_data) {
  mw_per_element <- sapply(colnames(elem_data$E), function(el) {
    if (el %in% names(ATOMIC_WEIGHTS)) ATOMIC_WEIGHTS[el] else 0
  })
  mw_per_met <- as.numeric(elem_data$E %*% mw_per_element)
  is_internal <- internal_mask(mod)
  S_int <- mod@S[, is_internal, drop = FALSE]
  no_formula  <- !elem_data$has_formula
  uncheckable <- as.numeric(t(abs(S_int)) %*% as.numeric(no_formula)) > 0
  mass_sum    <- as.numeric(t(S_int) %*% mw_per_met)
  n_total      <- ncol(S_int)
  n_checkable  <- sum(!uncheckable)
  n_consistent <- sum(!uncheckable & abs(mass_sum) < 0.001)
  pct          <- if (n_checkable > 0) n_consistent / n_checkable * 100 else NA_real_
  list(consistent_pct = pct, n_consistent = n_consistent,
       n_checkable = n_checkable, n_total = n_total)
}

check_element_balance <- function(mod, elem_data) {
  is_internal <- internal_mask(mod)
  S_int <- mod@S[, is_internal, drop = FALSE]
  no_formula  <- !elem_data$has_formula
  uncheckable <- as.numeric(t(abs(S_int)) %*% as.numeric(no_formula)) > 0
  E_subset       <- elem_data$E[, KEY_ELEMENTS, drop = FALSE]
  balance_matrix <- as.matrix(t(S_int) %*% E_subset)
  is_balanced    <- apply(abs(balance_matrix), 1, function(x) all(x < 0.001))
  n_total     <- ncol(S_int)
  n_checkable <- sum(!uncheckable)
  n_balanced  <- sum(!uncheckable & is_balanced)
  pct         <- if (n_checkable > 0) n_balanced / n_checkable * 100 else NA_real_
  list(balanced_pct = pct, n_balanced = n_balanced,
       n_checkable = n_checkable, n_total = n_total)
}

qc_one_model <- function(model_path, organism_id) {
  cat(sprintf("\n========== %s ==========\n", organism_id))
  mod <- readRDS(model_path)
  formulas <- NULL
  if ("met_attr" %in% slotNames(mod)) {
    if ("chemicalFormula" %in% names(mod@met_attr)) {
      formulas <- mod@met_attr$chemicalFormula
    } else if ("formula" %in% names(mod@met_attr)) {
      formulas <- mod@met_attr$formula
    }
  }
  if (is.null(formulas)) formulas <- rep(NA_character_, length(mod@met_id))
  all_elements <- unique(c(KEY_ELEMENTS, names(ATOMIC_WEIGHTS)))
  elem_data    <- build_element_matrix(formulas, all_elements)
  cat(sprintf("    formulas: %d / %d\n",
              sum(elem_data$has_formula), length(formulas)))

  c1 <- check_atp_from_water(mod)
  cat(sprintf("    [1] ATP-from-water: %g (%s)\n", c1$value, c1$status))

  c2 <- check_loopless(mod)
  cat(sprintf("    [2] Loopless: NA (%s)\n", c2$status))

  c3 <- check_total_mass_balance(mod, elem_data)
  cat(sprintf("    [3] Mass balance: %.2f%% (%d/%d)\n",
              c3$consistent_pct, c3$n_consistent, c3$n_checkable))

  c4 <- check_element_balance(mod, elem_data)
  cat(sprintf("    [4] Element balance: %.2f%% (%d/%d)\n",
              c4$balanced_pct, c4$n_balanced, c4$n_checkable))

  s1 <- classify_lt(c1$value, 1e-6, 0.01)
  s2 <- "NA"
  s3 <- classify_gt(c3$consistent_pct, 99, 95)
  s4 <- classify_gt(c4$balanced_pct, 95, 80)

  per_check <- c(s1, s2, s3, s4)
  overall <- if (any(per_check == "FAIL")) {
    paste0("FAIL_", paste(which(per_check == "FAIL"), collapse=","))
  } else if (any(per_check == "WARN")) {
    paste0("WARN_", paste(which(per_check == "WARN"), collapse=","))
  } else {
    "PASS"
  }

  data.frame(
    organism              = organism_id,
    atp_from_water        = c1$value,
    atp_status            = c1$status,
    atp_class             = s1,
    loopless_penalty_pct  = c2$penalty_pct,
    loopless_status       = c2$status,
    loopless_class        = s2,
    stoich_consistent_pct = c3$consistent_pct,
    stoich_n_checkable    = c3$n_checkable,
    stoich_class          = s3,
    mass_balanced_pct     = c4$balanced_pct,
    mass_n_checkable      = c4$n_checkable,
    mass_class            = s4,
    overall               = overall,
    stringsAsFactors      = FALSE
  )
}

args     <- commandArgs(trailingOnly = TRUE)
base_dir <- if (length(args) >= 1) args[1] else
            "/users/40335635/sharedscratch/EggNOG/manuscript_4/unculturable"
out_tsv  <- if (length(args) >= 2) args[2] else
            file.path(base_dir, "output_tsv", "qc_battery_results.tsv")

cat(sprintf("Base directory : %s\n", base_dir))
cat(sprintf("Output TSV     : %s\n", out_tsv))

all_files <- list.files(base_dir, pattern = "panModel\\.RDS$",
                        recursive = TRUE, full.names = TRUE)
all_files <- all_files[grep("/pan_model/panModel\\.RDS$", all_files)]
cat(sprintf("Found %d pan models\n", length(all_files)))

results_list <- list()
for (f in all_files) {
  parts   <- strsplit(f, "/")[[1]]
  rep_dir <- parts[length(parts) - 2]
  org_id  <- sub("_rep$", "_pan", rep_dir)
  res <- tryCatch(qc_one_model(f, org_id),
                  error = function(e) {
                    cat(sprintf("ERROR for %s: %s\n", org_id, conditionMessage(e)))
                    NULL
                  })
  if (!is.null(res)) results_list[[length(results_list) + 1]] <- res
  if (length(results_list) > 0) {
    tmp <- do.call(rbind, results_list)
    write.table(tmp, out_tsv, sep = "\t", quote = FALSE, row.names = FALSE)
  }
}

results <- do.call(rbind, results_list)
cat(sprintf("\n>>> Wrote %d rows to %s\n", nrow(results), out_tsv))
