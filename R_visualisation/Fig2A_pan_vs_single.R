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
# Statistical comparison: paired Wilcoxon signed-rank test (single
# vs pan within each SGB; n_pairs = 9). p-values are additionally
# adjusted with the Benjamini-Hochberg procedure across the four
# metrics; both raw p and BH q are written to the CSV companion file
# and printed in each panel header. Non-parametric is used because
# of small n and the heavy-tailed distribution of trace_bg_growth.
#
# All 9 SGBs are now included. The minimal-medium pan-fill of
# MGYG000291777, which previously failed to converge, has been
# rerun and produces normal summary metrics in the current Excel.
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
col_single <- "#95A5A6"
col_pan    <- "#2E86AB"
col_delta  <- "#34495E"
col_p_sig  <- "#1A5490"
col_p_ns   <- "#7F8C8D"

# ---------- SGB list (all 9, including 291777) ----------
sgbs <- c("MGYG000290784","MGYG000291361","MGYG000291777","MGYG000292637",
          "MGYG000293427","MGYG000294127","MGYG000295164","MGYG000295308","MGYG000295316")

display_order <- c(
  "MGYG000290784",  # Bacteroidota - UBA4334
  "MGYG000292637",  # Bacteroidota - RC9
  "MGYG000291361",  # Firmicutes  - UBA1777
  "MGYG000291777",  # Firmicutes  - RUG420
  "MGYG000293427",  # Firmicutes  - CAG-791
  "MGYG000294127",  # Firmicutes  - CAG-791
  "MGYG000295164",  # Firmicutes  - Ruminococcus_E
  "MGYG000295308",  # Firmicutes  - CAG-791
  "MGYG000295316"   # Firmicutes  - Ruminococcus_E
)

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

ordered_sgbs <- intersect(display_order, sgbs)

df <- bind_rows(single_dat, pan_dat) %>%
  mutate(
    model = factor(model, levels = c("Single rep MAG", "Pan model")),
    label = factor(display_labels[sgb], levels = rev(display_labels[ordered_sgbs]))
  )

df_long <- df %>%
  pivot_longer(c(n_reactions, n_auxotrophies, n_active_substrates, trace_bg_growth),
               names_to = "metric", values_to = "value")

df_delta <- df_long %>%
  pivot_wider(id_cols = c(sgb, label, metric),
              names_from = model, values_from = value) %>%
  mutate(delta_pct   = (`Pan model` - `Single rep MAG`) / `Single rep MAG` * 100,
         delta_label = sprintf("%+.0f%%", delta_pct))

# ---------- Paired Wilcoxon + BH correction ----------
metrics_list <- c("n_reactions", "n_auxotrophies",
                  "n_active_substrates", "trace_bg_growth")

stat_tab <- data.frame(metric = metrics_list,
                       n_pairs = NA_integer_,
                       median_pan_minus_single = NA_real_,
                       median_pct_change = NA_real_,
                       W = NA_real_, p_value = NA_real_,
                       stringsAsFactors = FALSE)

for (i in seq_along(metrics_list)) {
  m <- metrics_list[i]
  paired <- df_delta %>% filter(metric == m)
  x <- paired$`Single rep MAG`
  y <- paired$`Pan model`
  wt <- suppressWarnings(wilcox.test(y, x, paired = TRUE, exact = FALSE))
  stat_tab$n_pairs[i]                 <- length(x)
  stat_tab$median_pan_minus_single[i] <- median(y - x)
  stat_tab$median_pct_change[i]       <- median((y - x) / x * 100)
  stat_tab$W[i]                       <- unname(wt$statistic)
  stat_tab$p_value[i]                 <- wt$p.value
}

# Benjamini-Hochberg adjustment across the 4 metrics
stat_tab$q_value <- p.adjust(stat_tab$p_value, method = "BH")

# Display formatter
fmt_p <- function(p, q) {
  star <- if      (q < 0.001) "***"
          else if (q < 0.01)  "**"
          else if (q < 0.05)  "*"
          else                 "n.s."
  pstr <- if (p < 0.001) "p < 0.001" else sprintf("p = %.3f", p)
  qstr <- if (q < 0.001) "q < 0.001" else sprintf("q = %.3f", q)
  sprintf("%s, %s %s", pstr, qstr, star)
}
stat_tab$display <- mapply(fmt_p, stat_tab$p_value, stat_tab$q_value)
stat_tab$is_sig  <- stat_tab$q_value < 0.05

cat("\n=== Paired Wilcoxon signed-rank tests (single vs pan, n=9) ===\n")
print(stat_tab, row.names = FALSE, digits = 4)

# ---------- Panel maker ----------
make_panel <- function(metric_col, x_label, value_fmt = "%.0f") {
  pdat <- df_long  %>% filter(metric == metric_col)
  ddat <- df_delta %>% filter(metric == metric_col)
  s    <- stat_tab[stat_tab$metric == metric_col, ]
  
  pmax <- max(pdat$value, na.rm = TRUE)
  x_delta <- pmax * 1.18
  
  p_text   <- sprintf("Wilcoxon (paired): %s", s$display)
  p_colour <- if (s$is_sig) col_p_sig else col_p_ns
  
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
    annotate("text",
             x = pmax * 1.40, y = length(ordered_sgbs) + 0.7,
             label = p_text,
             hjust = 1, vjust = 0,
             size = 3.0, fontface = "bold",
             colour = p_colour) +
    scale_fill_manual(values = c("Single rep MAG" = col_single,
                                  "Pan model"      = col_pan),
                      name = NULL) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.32)),
                       limits = c(0, pmax * 1.40)) +
    coord_cartesian(clip = "off") +
    labs(x = x_label, y = NULL) +
    theme_classic(base_size = 10) +
    theme(
      axis.text.y = element_text(size = 9, colour = "grey20"),
      plot.margin = margin(18, 5, 5, 5)
    )
}

# ---------- Build 4 panels ----------
p_a <- make_panel("n_reactions",         "Reactions (n)",          "%d")
p_b <- make_panel("n_auxotrophies",      "Auxotrophic factors (n)", "%d")
p_c <- make_panel("n_active_substrates", "Active substrates (n)",   "%d")
p_d <- make_panel("trace_bg_growth",
                  expression(paste("Trace background growth (h"^-1, ")")),
                  "%.3f")

# ---------- Combine ----------
fig2 <- (p_a + p_b) / (p_c + p_d) +
  plot_layout(guides = "collect") +
  plot_annotation(
    tag_levels = "A",
    caption    = sprintf("Paired Wilcoxon signed-rank test, n = %d SGB pairs. Raw p and Benjamini-Hochberg adjusted q reported per panel; significance refers to q < 0.05. * q < 0.05, ** q < 0.01, *** q < 0.001.",
                         stat_tab$n_pairs[1])
  ) &
  theme(legend.position = "bottom",
        legend.text     = element_text(size = 9),
        legend.key.size = unit(0.4, "cm"),
        plot.tag        = element_text(face = "bold", size = 12),
        plot.caption    = element_text(size = 8.5, colour = "grey30",
                                       hjust = 0, margin = margin(t = 8)))

# ---------- Save ----------
out_path <- file.path(out_dir, "Fig2_pan_vs_single.png")
ggsave(out_path, fig2,
       width = 13, height = 9.0, units = "in",
       dpi = 1000, bg = "white")

cat(sprintf("\nSaved: %s   (%d SGBs included)\n", out_path, length(sgbs)))

# Save stat table
stat_out <- file.path(out_dir, "Fig2_paired_wilcoxon_results.csv")
write.csv(stat_tab[, c("metric","n_pairs","median_pan_minus_single",
                       "median_pct_change","W","p_value","q_value","display")],
          stat_out, row.names = FALSE)
cat(sprintf("Saved: %s\n", stat_out))

# ---------- Per-SGB delta print ----------
cat("\n=== Per-SGB Pan vs Single deltas ===\n")
df_delta %>%
  select(sgb, metric, `Single rep MAG`, `Pan model`, delta_pct) %>%
  arrange(metric, sgb) %>%
  print(n = Inf)
