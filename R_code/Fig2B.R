# =============================================================
# Figure 1B (8-row version)
# 6 references (PCG/Rumen paired by phylum) + 2 ARG-MAGs
#
# Two stacked rows, four metrics arranged as paired diverging bars:
#   Row 1: Reactions (L)            <- |organism| -> Active substrates (R)
#   Row 2: Auxotrophic factors (L)  <- |organism| -> Trace bg growth   (R)
# =============================================================

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(cowplot)
})

# ---------- Configuration ----------
xlsx_path <- "C:/Users/CFL/OneDrive - Queen's University Belfast/PhD_onedrive/Manuscript4/M4_Supplement.xlsx"
out_dir   <- "C:/Users/CFL/OneDrive - Queen's University Belfast/PhD_onedrive/Manuscript4/figures"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# Group palette (locked across the whole manuscript)
grp_pal <- c(
  "Gapseq Reference"   = "#D4A017",   # gold
  "Rumen Reference" = "#2E7D32",   # green
  "MAG Target"      = "#B22222"    # red
)

# ---------- Read data ----------
summ <- read_excel(xlsx_path, sheet = "single_all_summary")

# Top-to-bottom order: paired by phylum
# (Pseudomonadota -> Bacteroidota -> Bacillota -> MAGs)
target_ids <- c(
  "Ecoli_K12_MG1655",
  "Ecoli_PA3",
  "Btheta_VPI5482",
  "Btheta_KPPR3",
  "Efaecalis_ATCC19433",
  "Efaecalis_68A",
  "MGYG000292883",
  "MGYG000295553"
)

display_labels <- c(
  Ecoli_K12_MG1655    = "E. coli K-12 MG1655",
  Ecoli_PA3           = "E. coli PA3",
  Btheta_VPI5482      = "B. thetaiotaomicron VPI-5482",
  Btheta_KPPR3        = "B. thetaiotaomicron KPPR3",
  Efaecalis_ATCC19433 = "E. faecalis ATCC 19433",
  Efaecalis_68A       = "E. faecalis 68A",
  MGYG000292883       = "MGYG000292883",
  MGYG000295553       = "MGYG000295553"
)

group_map <- c(
  Ecoli_K12_MG1655    = "Gapseq Reference",
  Ecoli_PA3           = "Rumen Reference",
  Btheta_VPI5482      = "Gapseq Reference",
  Btheta_KPPR3        = "Rumen Reference",
  Efaecalis_ATCC19433 = "Gapseq Reference",
  Efaecalis_68A       = "Rumen Reference",
  MGYG000292883       = "MAG Target",
  MGYG000295553       = "MAG Target"
)

# Tolerate missing rows (e.g. when SLURM array still running)
missing_ids <- setdiff(target_ids, summ$organism)
if (length(missing_ids) > 0) {
  warning("Not yet present in single_all_summary (figure will skip them):\n  ",
          paste(missing_ids, collapse = "\n  "))
}
present_ids <- intersect(target_ids, summ$organism)

dat <- summ %>%
  filter(organism %in% present_ids) %>%
  arrange(match(organism, target_ids)) %>%
  mutate(
    label = factor(display_labels[organism],
                   levels = rev(display_labels[present_ids])),
    group = factor(group_map[organism], levels = names(grp_pal)),
    across(c(n_reactions, n_auxotrophies, n_active_substrates,
             trace_bg_growth), as.numeric)
  )

# ---------- Helper: left-extending panel ----------
make_left_panel <- function(df, var, axis_label, label_fmt = "%d") {
  vmax <- max(df[[var]])
  ggplot(df, aes(x = -.data[[var]], y = label, fill = group)) +
    geom_col(width = 0.65, colour = "grey20", linewidth = 0.3) +
    geom_text(aes(label = sprintf(label_fmt, .data[[var]])),
              hjust = 1.15, size = 3.0, colour = "grey20") +
    scale_x_continuous(
      expand = expansion(mult = c(0.22, 0)),
      limits = c(-vmax * 1.30, 0),
      labels = function(x) abs(x),
      breaks = scales::pretty_breaks(n = 4)
    ) +
    scale_fill_manual(values = grp_pal, drop = FALSE) +
    labs(x = axis_label, y = NULL) +
    theme_classic(base_size = 11) +
    theme(
      axis.text.y     = element_blank(),
      axis.ticks.y    = element_blank(),
      axis.line.y     = element_blank(),
      legend.position = "none",
      plot.margin     = margin(5, 0, 5, 5)
    )
}

# ---------- Helper: right-extending panel ----------
make_right_panel <- function(df, var, axis_label, label_fmt = "%d") {
  vmax <- max(df[[var]])
  ggplot(df, aes(x = .data[[var]], y = label, fill = group)) +
    geom_col(width = 0.65, colour = "grey20", linewidth = 0.3) +
    geom_text(aes(label = sprintf(label_fmt, .data[[var]])),
              hjust = -0.15, size = 3.0, colour = "grey20") +
    scale_x_continuous(
      expand = expansion(mult = c(0, 0.22)),
      limits = c(0, vmax * 1.30),
      breaks = scales::pretty_breaks(n = 4)
    ) +
    scale_fill_manual(values = grp_pal, drop = FALSE) +
    labs(x = axis_label, y = NULL) +
    theme_classic(base_size = 11) +
    theme(
      axis.text.y     = element_blank(),
      axis.ticks.y    = element_blank(),
      axis.line.y     = element_blank(),
      legend.position = "none",
      plot.margin     = margin(5, 5, 5, 0)
    )
}

# ---------- Helper: centre label panel ----------
make_centre_panel <- function(df) {
  ggplot(df, aes(x = 0, y = label, label = label)) +
    geom_text(size = 3.2, colour = "grey15") +
    scale_x_continuous(limits = c(-1, 1), expand = c(0, 0)) +
    theme_void() +
    theme(plot.margin = margin(5, 0, 5, 0))
}

# ---------- Build 4 metric panels + 2 centre panels ----------
p_top_left  <- make_left_panel (dat, "n_reactions",         "Reactions",          "%d")
p_top_right <- make_right_panel(dat, "n_active_substrates", "Active substrates",  "%d")
p_bot_left  <- make_left_panel (dat, "n_auxotrophies",      "Auxotrophic factors", "%d")
p_bot_right <- make_right_panel(dat, "trace_bg_growth",
                                expression(paste("Trace background growth (h"^-1, ")")),
                                "%.3f")

p_centre_top <- make_centre_panel(dat)
p_centre_bot <- make_centre_panel(dat)

# ---------- Assemble two stacked rows ----------
row_top <- p_top_left + p_centre_top + p_top_right +
            plot_layout(widths = c(1, 0.85, 1))
row_bot <- p_bot_left + p_centre_bot + p_bot_right +
            plot_layout(widths = c(1, 0.85, 1))

# ---------- Shared group legend ----------
legend_dummy <- ggplot(dat, aes(x = 1, y = label, fill = group)) +
  geom_col() +
  scale_fill_manual(values = grp_pal, name = NULL, drop = FALSE) +
  theme_void() +
  theme(legend.position = "bottom",
        legend.key.size = unit(0.4, "cm"),
        legend.text     = element_text(size = 9))

shared_legend <- cowplot::get_legend(legend_dummy)

# ---------- Combine ----------
fig1b <- row_top / row_bot / patchwork::wrap_elements(shared_legend) +
  plot_layout(heights = c(1, 1, 0.06))

# ---------- Save ----------
out_path <- file.path(out_dir, "Fig2B.png")
ggsave(out_path, fig1b,
       width = 9.5, height = 7.8, units = "in",
       dpi = 1000, bg = "white")

cat("Saved:", out_path, "\n")

