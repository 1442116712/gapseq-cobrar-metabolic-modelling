# =============================================================
# Figure 1A (8-row version)
# 6 references (PCG/Rumen paired by phylum) + 2 ARG-MAGs
# Diverging bar layout:
#   Left   = ORF coverage (% extending leftwards)
#   Centre = organism label
#   Right  = Stage 1 added reactions, stacked core / non-core
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

# Group palette (3 groups, locked across all manuscript figures)
grp_pal <- c(
  "Gapseq Reference"   = "#D4A017",   # gold
  "Rumen Reference" = "#2E7D32",   # green
  "MAG Target"      = "#B22222"    # red
)

# Core / non-core stack palette (kept from 5-row version)
type_pal <- c(
  "Core (genomic evidence)" = "#5B9BD5",
  "Non-core (topological)"  = "#C9D7E8"
)

# ---------- Read data ----------
recon <- read_excel(xlsx_path, sheet = "single_all_reconstruction")

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
missing_ids <- setdiff(target_ids, recon$organism)
if (length(missing_ids) > 0) {
  warning("Not yet present in single_all_reconstruction (figure will skip them):\n  ",
          paste(missing_ids, collapse = "\n  "))
}
present_ids <- intersect(target_ids, recon$organism)

dat <- recon %>%
  filter(organism %in% present_ids) %>%
  arrange(match(organism, target_ids)) %>%
  mutate(
    label        = factor(display_labels[organism],
                          levels = rev(display_labels[present_ids])),
    group        = factor(group_map[organism], levels = names(grp_pal)),
    s1_core_n    = as.numeric(s1_core),
    s1_added_n   = as.numeric(s1_added),
    s1_noncore_n = s1_added_n - s1_core_n,
    orf_cov_n    = as.numeric(orf_coverage_pct)
  )

# ---------- Left panel: ORF coverage extending leftwards ----------
p_left <- ggplot(dat, aes(x = -orf_cov_n, y = label, fill = group)) +
  geom_col(width = 0.65, colour = "grey20", linewidth = 0.3) +
  geom_text(aes(label = sprintf("%.1f%%", orf_cov_n)),
            hjust = 1.15, size = 3.0, colour = "grey20") +
  scale_x_continuous(
    expand = expansion(mult = c(0.18, 0)),
    limits = c(-max(dat$orf_cov_n) * 1.25, 0),
    labels = function(x) abs(x),
    breaks = scales::pretty_breaks(n = 4)
  ) +
  scale_fill_manual(values = grp_pal, name = NULL, drop = FALSE) +
  labs(x = "ORF coverage (%)", y = NULL) +
  theme_classic(base_size = 11) +
  theme(
    axis.text.y  = element_blank(),
    axis.ticks.y = element_blank(),
    axis.line.y  = element_blank(),
    plot.margin  = margin(5, 0, 5, 5)
  )

# ---------- Centre panel: organism labels ----------
p_centre <- ggplot(dat, aes(x = 0, y = label, label = label)) +
  geom_text(size = 3.2, colour = "grey15") +
  scale_x_continuous(limits = c(-1, 1), expand = c(0, 0)) +
  theme_void() +
  theme(plot.margin = margin(5, 0, 5, 0))

# ---------- Right panel: Stage 1 added reactions, stacked ----------
dat_long <- dat %>%
  select(label, s1_core_n, s1_noncore_n) %>%
  pivot_longer(cols = c(s1_core_n, s1_noncore_n),
               names_to = "rxn_type", values_to = "n") %>%
  mutate(rxn_type = factor(
    ifelse(rxn_type == "s1_core_n",
           "Core (genomic evidence)",
           "Non-core (topological)"),
    levels = names(type_pal)
  ))

totals <- dat %>% transmute(label, total = s1_added_n)

p_right <- ggplot(dat_long, aes(x = n, y = label, fill = rxn_type)) +
  geom_col(width = 0.65, colour = "grey20", linewidth = 0.3) +
  geom_text(data = totals, aes(x = total, y = label, label = total),
            inherit.aes = FALSE, hjust = -0.25, size = 3.0,
            colour = "grey20") +
  scale_x_continuous(expand = expansion(mult = c(0, 0.15)),
                     limits = c(0, max(totals$total) * 1.18)) +
  scale_fill_manual(values = type_pal, name = NULL) +
  labs(x = "Stage 1 added reactions", y = NULL) +
  theme_classic(base_size = 11) +
  theme(
    axis.text.y  = element_blank(),
    axis.ticks.y = element_blank(),
    axis.line.y  = element_blank(),
    plot.margin  = margin(5, 5, 5, 0)
  )

# ---------- Combine with collected legends ----------
fig1a <- p_left + p_centre + p_right +
  plot_layout(widths = c(1, 0.75, 1.2), guides = "collect") &
  theme(legend.position  = "bottom",
        legend.box       = "horizontal",
        legend.spacing.x = unit(0.6, "cm"),
        legend.text      = element_text(size = 9),
        legend.key.size  = unit(0.4, "cm"))

# ---------- Save ----------
out_path <- file.path(out_dir, "Fig2A.png")
ggsave(out_path, fig1a,
       width = 9.0, height = 4.6, units = "in",
       dpi = 1000, bg = "white")

cat("Saved:", out_path, "\n")

