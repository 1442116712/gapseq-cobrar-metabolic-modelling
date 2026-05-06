# =============================================================
# Comparison plot - MGYG000291777 pan-Draft
#   TSB-medium gap-fill versus minimal-medium gap-fill
#
# Reads pandraft_all_reconstruction sheet and contrasts the two
# rows: MGYG000291777_pan (TSB) and MGYG000291777_pan_minimal.
#
# Two panels:
#   A | Growth rate at each gap-fill stage (Stage 1, Stage 3, Final)
#   B | Top secreted products (parsed from top_products field,
#       compounds in uptake_at_limit are excluded, zero-flux entries
#       are dropped)
#
# Output: Fig_TSB_vs_minimal_MGYG000291777.png (10 x 4.5 in, 1000 dpi)
# =============================================================

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(stringr)
  library(ggplot2)
  library(patchwork)
})

# ---------- Configuration ----------
xlsx_path <- "C:/Users/CFL/OneDrive - Queen's University Belfast/PhD_onedrive/Manuscript4/M4_Supplement.xlsx"
out_dir   <- "C:/Users/CFL/OneDrive - Queen's University Belfast/PhD_onedrive/Manuscript4/figures"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

TOP_N        <- 10           # Top secreted products to display
TSB_COL      <- "#2E86AB"    # TSB-fill (main figure version)
MIN_COL      <- "#95A5A6"    # Minimal-fill (failed control)
TSB_LABEL    <- "TSB-fill"
MIN_LABEL    <- "Minimal-fill"

# ---------- Read data ----------
recon <- read_excel(xlsx_path, sheet = "pandraft_all_reconstruction")
tsb   <- recon %>% filter(organism == "MGYG000291777_pan")
mn    <- recon %>% filter(organism == "MGYG000291777_pan_minimal")

if (nrow(tsb) != 1 || nrow(mn) != 1) {
  stop("Expected exactly one row each for MGYG000291777_pan and ",
       "MGYG000291777_pan_minimal in pandraft_all_reconstruction.")
}

# ---------- Helper: parse "compound:flux, compound:flux" strings ----------
parse_compounds <- function(s) {
  if (is.na(s) || nchar(trimws(s)) == 0) {
    return(tibble(compound = character(), flux = double()))
  }
  m <- str_match_all(s, "([A-Za-z][A-Za-z0-9+\\-]*):([0-9]+\\.?[0-9]*)")[[1]]
  if (nrow(m) == 0) {
    return(tibble(compound = character(), flux = double()))
  }
  tibble(compound = m[, 2], flux = as.numeric(m[, 3])) %>%
    group_by(compound) %>%
    summarise(flux = max(flux, na.rm = TRUE), .groups = "drop")
}

# ---------- Panel A: stage-wise growth ----------
growth_df <- tibble(
  stage   = factor(c("Stage 1", "Stage 3", "Final"),
                   levels = c("Stage 1", "Stage 3", "Final")),
  TSB     = c(tsb$s1_growth, tsb$s3_growth, tsb$final_growth),
  Minimal = c(mn$s1_growth,  mn$s3_growth,  mn$final_growth)
) %>%
  pivot_longer(c(TSB, Minimal),
               names_to = "medium", values_to = "growth") %>%
  mutate(medium = factor(medium, levels = c("TSB", "Minimal")))

p_a <- ggplot(growth_df, aes(x = stage, y = growth, fill = medium)) +
  geom_col(position = position_dodge(width = 0.75),
           width = 0.65, colour = "grey20", linewidth = 0.3) +
  geom_text(aes(label = sprintf("%.2f", growth)),
            position = position_dodge(width = 0.75),
            vjust = -0.4, size = 3.4, fontface = "bold") +
  scale_fill_manual(values = c(TSB = TSB_COL, Minimal = MIN_COL),
                    labels = c(TSB = TSB_LABEL, Minimal = MIN_LABEL),
                    name = NULL) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.18))) +
  labs(title = "A. Growth recovery across gap-fill stages",
       x = NULL,
       y = expression(bold(paste("Growth rate (h"^-1, ")")))) +
  theme_classic(base_size = 12) +
  theme(
    plot.title       = element_text(face = "bold", size = 13,
                                    hjust = 0, margin = margin(b = 8)),
    axis.text        = element_text(face = "bold", colour = "grey15"),
    axis.title.y     = element_text(face = "bold", size = 11),
    legend.position  = "top",
    legend.text      = element_text(face = "bold", size = 10),
    legend.key.size  = unit(0.5, "cm")
  )

# ---------- Panel B: top secreted products ----------
tsb_uptakes <- parse_compounds(tsb$uptake_at_limit)$compound
mn_uptakes  <- parse_compounds(mn$uptake_at_limit)$compound

tsb_prod <- parse_compounds(tsb$top_products) %>%
  filter(!compound %in% tsb_uptakes, flux > 0) %>%
  mutate(medium = "TSB")
mn_prod  <- parse_compounds(mn$top_products) %>%
  filter(!compound %in% mn_uptakes, flux > 0) %>%
  mutate(medium = "Minimal")

prod_all <- bind_rows(tsb_prod, mn_prod)

top_compounds <- prod_all %>%
  group_by(compound) %>%
  summarise(max_flux = max(flux), .groups = "drop") %>%
  arrange(desc(max_flux)) %>%
  slice_head(n = TOP_N) %>%
  pull(compound)

prod_plot <- expand.grid(compound = top_compounds,
                         medium   = c("TSB", "Minimal"),
                         stringsAsFactors = FALSE) %>%
  left_join(prod_all, by = c("compound", "medium")) %>%
  mutate(flux     = ifelse(is.na(flux), 0, flux),
         compound = factor(compound, levels = rev(top_compounds)),
         medium   = factor(medium, levels = c("TSB", "Minimal")))

p_b <- ggplot(prod_plot, aes(x = flux, y = compound, fill = medium)) +
  geom_col(position = position_dodge(width = 0.7),
           width = 0.65, colour = "grey20", linewidth = 0.3) +
  geom_text(data = prod_plot %>% filter(flux > 0),
            aes(label = sprintf("%.2f", flux)),
            position = position_dodge(width = 0.7),
            hjust = -0.15, size = 3.0, fontface = "bold") +
  scale_fill_manual(values = c(TSB = TSB_COL, Minimal = MIN_COL),
                    labels = c(TSB = TSB_LABEL, Minimal = MIN_LABEL),
                    name = NULL) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.18))) +
  labs(title = "B. Top secreted products",
       x = expression(bold(paste("Flux (mmol gDW"^-1, " h"^-1, ")"))),
       y = NULL) +
  theme_classic(base_size = 12) +
  theme(
    plot.title       = element_text(face = "bold", size = 13,
                                    hjust = 0, margin = margin(b = 8)),
    axis.text.y      = element_text(face = "bold", colour = "grey15", size = 10),
    axis.text.x      = element_text(face = "bold", colour = "grey15"),
    axis.title.x     = element_text(face = "bold", size = 11),
    legend.position  = "top",
    legend.text      = element_text(face = "bold", size = 10),
    legend.key.size  = unit(0.5, "cm")
  )

# ---------- Combine and save ----------
fig <- p_a + p_b +
  plot_layout(ncol = 2, widths = c(1, 1.4), guides = "collect") &
  theme(legend.position = "top")

fig <- fig +
  plot_annotation(
    title = "MGYG000291777 pan-Draft: TSB versus minimal medium gap-fill",
    theme = theme(plot.title = element_text(face = "bold", size = 14,
                                            hjust = 0.5))
  )

out_path <- file.path(out_dir, "Fig_TSB_vs_minimal_MGYG000291777.png")
ggsave(out_path, fig,
       width = 11, height = 4.8, units = "in",
       dpi = 1000, bg = "white")

cat(sprintf("Saved: %s\n", out_path))

# ---------- Diagnostics ----------
cat("\n=== Stage-wise growth ===\n")
print(growth_df)

cat("\n=== Top secreted products (filtered) ===\n")
print(prod_plot %>%
        pivot_wider(names_from = medium, values_from = flux,
                    values_fill = 0) %>%
        arrange(match(compound, top_compounds)))

cat("\n=== Excluded uptake compounds ===\n")
cat(sprintf("  TSB uptake_at_limit:     %s\n",
            paste(tsb_uptakes, collapse = ", ")))
cat(sprintf("  Minimal uptake_at_limit: %s\n",
            paste(mn_uptakes, collapse = ", ")))

