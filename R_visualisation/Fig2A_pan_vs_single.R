# =============================================================
# Figure 2: Pan-model vs single-rep comparison across SGBs
#
# Four panels arranged in a 2x2 grid, one per metric:
#   A. Reactions (n)
#   B. Auxotrophic factors (n)
#   C. Active substrates (n)
#   D. Trace background growth (h^-1)
#
# Each panel: horizontal grouped bars (single rep MAG vs pan model)
# for the SGBs, with bar values on each bar and pan/single delta% at
# the right edge.
#
# MGYG000291777 is excluded from the main figure: its minimal-medium
# pan-fill failed to converge (trace bg growth = 1.24 h^-1 and only
# n_active_substrates = 2 in the current supplementary table). The
# manuscript narrative uses MGYG000291777_pan_TSB instead, but those
# summary metrics are not in the current Excel; flip the flag below to
# include the minimal-medium row for diagnostic plotting.
# =============================================================

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
})

# ---------- Configuration ----------
xlsx_path <- "C:/Users/CFL/OneDrive - Queen's University Belfast/PhD_onedrive/Manuscript4/M4_Supplement.xlsx"
out_dir   <- "C:/Users/CFL/OneDrive - Queen's University Belfast/PhD_onedrive/Manuscript4/figures"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# Colours
col_single <- "#95A5A6"   # neutral grey, baseline single-rep
col_pan    <- "#2E86AB"   # teal blue, pan model
col_delta  <- "#34495E"   # dark slate, neutral delta% text

# ---------- SGB list ----------
all_sgbs <- c("MGYG000290784","MGYG000291361","MGYG000291777","MGYG000292637",
              "MGYG000293427","MGYG000294127","MGYG000295164","MGYG000295308","MGYG000295316")

EXCLUDE_PROBLEMATIC <- TRUE
sgbs <- if (EXCLUDE_PROBLEMATIC) setdiff(all_sgbs, "MGYG000291777") else all_sgbs

# Phylum-grouped top-to-bottom display order (Bacteroidota first, then Firmicutes)
display_order <- c(
  "MGYG000290784",  # Bacteroidota - UBA4334
  "MGYG000292637",  # Bacteroidota - RC9
  "MGYG000291361",  # Firmicutes  - UBA1777
  "MGYG000291777",  # Firmicutes  - RUG420   (excluded by default)
  "MGYG000293427",  # Firmicutes  - CAG-791
  "MGYG000294127",  # Firmicutes  - CAG-791
  "MGYG000295164",  # Firmicutes  - Ruminococcus_E
  "MGYG000295308",  # Firmicutes  - CAG-791
  "MGYG000295316"   # Firmicutes  - Ruminococcus_E
)

# Display labels with genus in parentheses for biological context
display_labels <- c(
  MGYG000290784 = "MGYG000290784 (UBA4334)",
  MGYG000292637 = "MGYG000292637 (RC9)",
  MGYG000291361 = "MGYG000291361 (UBA1777)",
  MGYG000291777 = "MGYG000291777 (RUG420)",
  MGYG000293427 = "MGYG000293427 (CAG-791)",
  MGYG000294127 = "MGYG000294127 (CAG-791)",
  MGYG000295164 = "MGYG000295164 (Ruminococcus_E)",
  MGYG000295308 = "MGYG000295308 (CAG-791)",
  MGYG000295316 = "MGYG000295316 (Ruminococcus_E)"
)

# ---------- Read data ----------
single_summ <- read_excel(xlsx_path, sheet = "single_all_summary")
pan_summ    <- read_excel(xlsx_path, sheet = "pandraft_all_summary")

# ---------- Reshape ----------
single_dat <- single_summ %>%
  filter(organism %in% sgbs) %>%
  mutate(sgb = organism, model = "Single rep MAG") %>%
  select(sgb, model, n_reactions, n_auxotrophies,
         n_active_substrates, trace_bg_growth)

pan_dat <- pan_summ %>%
  mutate(sgb = sub("_pan$", "", organism)) %>%
  filter(sgb %in% sgbs) %>%
  mutate(model = "Pan model") %>%
  select(sgb, model, n_reactions, n_auxotrophies,
         n_active_substrates, trace_bg_growth)

# Final ordering used by the figure
ordered_sgbs <- intersect(display_order, sgbs)

df <- bind_rows(single_dat, pan_dat) %>%
  mutate(
    model = factor(model, levels = c("Single rep MAG", "Pan model")),
    label = factor(display_labels[sgb], levels = rev(display_labels[ordered_sgbs]))
  )

# Long format
df_long <- df %>%
  pivot_longer(c(n_reactions, n_auxotrophies, n_active_substrates, trace_bg_growth),
               names_to = "metric", values_to = "value")

# Per-metric deltas: pan vs single, percent
df_delta <- df_long %>%
  pivot_wider(id_cols = c(sgb, label, metric),
              names_from = model, values_from = value) %>%
  mutate(delta_pct   = (`Pan model` - `Single rep MAG`) / `Single rep MAG` * 100,
         delta_label = sprintf("%+.0f%%", delta_pct))

# ---------- Panel maker ----------
make_panel <- function(metric_col, x_label, value_fmt = "%.0f") {
  pdat <- df_long  %>% filter(metric == metric_col)
  ddat <- df_delta %>% filter(metric == metric_col)
  
  pmax <- max(pdat$value, na.rm = TRUE)
  x_delta <- pmax * 1.18
  
  ggplot(pdat, aes(x = value, y = label, fill = model)) +
    geom_col(position = position_dodge(width = 0.75), width = 0.65,
             colour = "grey20", linewidth = 0.3) +
    geom_text(aes(label = sprintf(value_fmt, value)),
              position = position_dodge(width = 0.75),
              hjust = -0.15, size = 2.5, colour = "grey20") +
    geom_text(data = ddat,
              aes(x = x_delta, y = label, label = delta_label),
              size = 2.8, fontface = "bold", colour = col_delta,
              inherit.aes = FALSE) +
    scale_fill_manual(values = c("Single rep MAG" = col_single,
                                  "Pan model"      = col_pan),
                      name = NULL) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.32)),
                       limits = c(0, pmax * 1.40)) +
    labs(x = x_label, y = NULL) +
    theme_classic(base_size = 10) +
    theme(
      axis.text.y = element_text(size = 9, colour = "grey20"),
      plot.margin = margin(5, 5, 5, 5)
    )
}

# ---------- Build 4 panels ----------
p_a <- make_panel("n_reactions",         "Reactions (n)",          "%d")
p_b <- make_panel("n_auxotrophies",      "Auxotrophic factors (n)", "%d")
p_c <- make_panel("n_active_substrates", "Active substrates (n)",   "%d")
p_d <- make_panel("trace_bg_growth",
                  expression(paste("Trace background growth (h"^-1, ")")),
                  "%.3f")

# ---------- Combine with shared legend and panel labels ----------
fig2 <- (p_a + p_b) / (p_c + p_d) +
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = "A") &
  theme(legend.position = "bottom",
        legend.text     = element_text(size = 9),
        legend.key.size = unit(0.4, "cm"),
        plot.tag        = element_text(face = "bold", size = 12))

# ---------- Save ----------
out_path <- file.path(out_dir, "Fig2_pan_vs_single.png")
ggsave(out_path, fig2,
       width = 13, height = 8.5, units = "in",
       dpi = 1000, bg = "white")

cat(sprintf("Saved: %s   (%d SGBs included)\n", out_path, length(sgbs)))

# ---------- Summary print of deltas ----------
cat("\n=== Pan vs Single delta summary ===\n")
df_delta %>%
  select(sgb, metric, `Single rep MAG`, `Pan model`, delta_pct) %>%
  arrange(metric, sgb) %>%
  print(n = Inf)

