# =============================================================================
# Fig_42_4metrics.R
# Manuscript 4 - Reconstruction quality + phenotypic capacity across all
# 42 metabolic models, with three role blocks visually separated:
#   Top:    6 reference single-genome models (3 paired clinical/rumen)
#   Mid:    2 uncultured ARG-MAG single-genome models
#   Bot:    34 pan-Draft species-representative reconstructions
#
# Two fill scales (via ggnewscale):
#   - Single models coloured by group (Gapseq/Rumen ref, ARG-MAG)
#   - Pan models    coloured by class (v4 palette, matches Fig 4)
#
# Four panels, equal width:
#   A. ORF coverage (%)   - pan value = mean ORF across contributing MAGs
#   B. Reactions (n)
#   C. Active substrates (n)
#   D. Auxotrophic factors (n)
#
# Y-axis blocks separated by spacer rows; a narrow brace panel on the
# very left labels the three blocks.
# =============================================================================

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(patchwork)
  library(ggnewscale)   # for two fill scales in one plot
  library(cowplot)      # for shared legends
})

# ---- Inputs ------------------------------------------------------------------
xlsx_path <- "C:/Users/CFL/OneDrive - Queen's University Belfast/PhD_onedrive/Manuscript4/M4_Supplement.xlsx"
out_dir   <- "C:/Users/CFL/OneDrive - Queen's University Belfast/PhD_onedrive/Manuscript4/figures"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ---- Read data ---------------------------------------------------------------
pan_stage1 <- read_excel(xlsx_path, sheet = "pandraft_all_stage1_summary")
sin_recon  <- read_excel(xlsx_path, sheet = "single_all_reconstruction")
pan_summ   <- read_excel(xlsx_path, sheet = "pandraft_all_summary")
sin_summ   <- read_excel(xlsx_path, sheet = "single_all_summary")
tg         <- read_excel(xlsx_path, sheet = "target-SGBs")

num_force <- function(df, cols) {
  for (c in cols) {
    if (c %in% colnames(df)) {
      df[[c]] <- suppressWarnings(as.numeric(df[[c]]))
    }
  }
  df
}
sin_recon  <- num_force(sin_recon,  c("orf_coverage_pct"))
pan_stage1 <- num_force(pan_stage1, c("orf_cov_mean","orf_cov_median","orf_cov_max"))
pan_summ   <- num_force(pan_summ,   c("n_reactions","n_active_substrates",
                                       "n_auxotrophies","trace_bg_growth"))
sin_summ   <- num_force(sin_summ,   c("n_reactions","n_active_substrates",
                                       "n_auxotrophies","trace_bg_growth"))

# Pan ORF coverage: mean across MAGs contributing to each pan-Draft
pan_orf <- pan_stage1 %>%
  transmute(organism = pan_id,
            orf_coverage_pct = orf_cov_mean)

# ---- Merge ORF + summary per organism ---------------------------------------
pan_full <- pan_orf %>%
  inner_join(pan_summ %>%
               select(organism, n_reactions, n_active_substrates,
                      n_auxotrophies),
             by = "organism")

# Singles: only keep the 8 we care about
keep_single <- c("Ecoli_K12_MG1655","Ecoli_PA3",
                 "Btheta_VPI5482","Btheta_KPPR3",
                 "Efaecalis_ATCC19433","Efaecalis_68A",
                 "MGYG000292883","MGYG000295553")
sin_recon <- sin_recon %>% filter(organism %in% keep_single)
sin_summ  <- sin_summ  %>% filter(organism %in% keep_single)

sin_full <- sin_recon %>%
  select(organism, orf_coverage_pct) %>%
  inner_join(sin_summ %>%
               select(organism, n_reactions, n_active_substrates,
                      n_auxotrophies),
             by = "organism")

# ---- Metadata: class + role + group -----------------------------------------
sp_tax <- tg %>%
  group_by(Species_rep) %>%
  summarise(Class = first(Class), .groups = "drop") %>%
  rename(sgb_id = Species_rep)

pan_meta <- tibble(organism = pan_full$organism) %>%
  mutate(sgb_id = sub("_pan$", "", organism)) %>%
  left_join(sp_tax, by = "sgb_id") %>%
  mutate(role  = "pan",
         group = NA_character_) %>%      # pan uses Class for fill, not group
  select(organism, Class, role, group)

sin_meta <- tribble(
  ~organism,             ~Class,                ~role,    ~group,
  "Btheta_VPI5482",      "Bacteroidia",         "ref",    "Gapseq Reference",
  "Btheta_KPPR3",        "Bacteroidia",         "ref",    "Rumen Reference",
  "Efaecalis_ATCC19433", "Bacilli",             "ref",    "Gapseq Reference",
  "Efaecalis_68A",       "Bacilli",             "ref",    "Rumen Reference",
  "Ecoli_K12_MG1655",    "Gammaproteobacteria", "ref",    "Gapseq Reference",
  "Ecoli_PA3",           "Gammaproteobacteria", "ref",    "Rumen Reference",
  "MGYG000292883",       "Lentisphaeria",       "argmag", "ARG-MAG",
  "MGYG000295553",       "Kiritimatiellia",     "argmag", "ARG-MAG")

all_dat <- bind_rows(pan_full, sin_full) %>%
  left_join(bind_rows(pan_meta, sin_meta), by = "organism")

# ---- Y-axis order: ref top, ARG mid, pan bottom, with spacers ----------------
ref_order <- c("Btheta_VPI5482","Btheta_KPPR3",
               "Efaecalis_ATCC19433","Efaecalis_68A",
               "Ecoli_K12_MG1655","Ecoli_PA3")
arg_order <- c("MGYG000292883","MGYG000295553")

class_order_pan <- c("Bacteroidia", "Clostridia", "Bacilli",
                     "Negativicutes", "Coriobacteriia",
                     "Methanobacteria")

pan_block <- all_dat %>%
  filter(role == "pan") %>%
  mutate(class_rank = match(Class, class_order_pan)) %>%
  arrange(class_rank, organism)
pan_order <- pan_block$organism

# Display labels - ARG-MAGs carry full class name in brackets
display_lookup <- c(
  Btheta_VPI5482      = "B. theta VPI-5482",
  Btheta_KPPR3        = "B. theta KPPR-3",
  Efaecalis_ATCC19433 = "E. faecalis ATCC 19433",
  Efaecalis_68A       = "E. faecalis 68A",
  Ecoli_K12_MG1655    = "E. coli K-12 MG1655",
  Ecoli_PA3           = "E. coli PA-3",
  MGYG000292883       = "MGYG000292883 (Lentisphaeria)",
  MGYG000295553       = "MGYG000295553 (Kiritimatiellia)"
)
all_dat <- all_dat %>%
  mutate(display = ifelse(organism %in% names(display_lookup),
                          display_lookup[organism],
                          sub("_pan$", "", organism)))

# Build ordered display vector + spacer labels (unique whitespace strings)
SPACER_1 <- " "    # one space
SPACER_2 <- "  "   # two spaces (must differ as factor levels)

ord_orgs <- c(ref_order, "_SPACER1_", arg_order, "_SPACER2_", pan_order)

ord_display <- vapply(ord_orgs, function(o) {
  if (o == "_SPACER1_") return(SPACER_1)
  if (o == "_SPACER2_") return(SPACER_2)
  if (o %in% names(display_lookup)) return(display_lookup[[o]])
  return(sub("_pan$", "", o))
}, character(1))

# Inject spacer rows
spacer_row <- function(label) {
  tibble(organism = paste0("_SPACER_", label),
         orf_coverage_pct = NA_real_,
         n_reactions = NA_real_,
         n_active_substrates = NA_real_,
         n_auxotrophies = NA_real_,
         Class = NA_character_,
         role = "spacer",
         group = NA_character_,
         display = label)
}

all_dat <- bind_rows(
  all_dat %>% filter(organism %in% ref_order),
  spacer_row(SPACER_1),
  all_dat %>% filter(organism %in% arg_order),
  spacer_row(SPACER_2),
  all_dat %>% filter(organism %in% pan_order)
) %>%
  mutate(display = factor(display, levels = rev(ord_display)))

# ---- Palettes ---------------------------------------------------------------
# Pan: class palette, only the classes present in the 34 pan-Draft set
pan_class_palette <- c(
  "Bacteroidia"     = "#4F8FA8",
  "Clostridia"      = "#E07A1F",
  "Bacilli"         = "#F2CC4F",
  "Negativicutes"   = "#7A8550",
  "Coriobacteriia"  = "#9B7BAE",
  "Methanobacteria" = "#A04A85"
)

# Single: role/group palette matching Fig 1A/1B
single_group_palette <- c(
  "Gapseq Reference" = "#D4A017",   # gold
  "Rumen Reference"  = "#2E7D32",   # green
  "ARG-MAG"          = "#B22222"    # red
)

# ---- Panel helper ------------------------------------------------------------
make_panel <- function(df, value_col, x_label, value_fmt = "%.0f",
                       show_y = FALSE) {
  df_plot   <- df %>% filter(!is.na(.data[[value_col]]))
  pan_sub   <- df_plot %>% filter(role == "pan")
  sing_sub  <- df_plot %>% filter(role %in% c("ref","argmag"))
  vmax      <- max(df_plot[[value_col]], na.rm = TRUE)

  p <- ggplot() +
    # Pan layer (Class fill)
    geom_col(data = pan_sub,
             aes(x = .data[[value_col]], y = display, fill = Class),
             width = 0.7, colour = "grey25", linewidth = 0.25) +
    scale_fill_manual(values = pan_class_palette,
                      name = "Pan-Draft (class)",
                      drop = FALSE,
                      na.translate = FALSE) +
    new_scale_fill() +
    # Single layer (group fill) - drawn on a fresh scale
    geom_col(data = sing_sub,
             aes(x = .data[[value_col]], y = display, fill = group),
             width = 0.7, colour = "grey25", linewidth = 0.25) +
    scale_fill_manual(values = single_group_palette,
                      name = "Single genome",
                      drop = FALSE,
                      na.translate = FALSE) +
    # Numeric labels (both layers)
    geom_text(data = df_plot,
              aes(x = .data[[value_col]], y = display,
                  label = sprintf(value_fmt, .data[[value_col]])),
              hjust = -0.15, size = 2.6, colour = "grey20") +
    scale_x_continuous(expand = expansion(mult = c(0, 0.18)),
                       limits = c(0, vmax * 1.20)) +
    scale_y_discrete(drop = FALSE) +
    labs(x = x_label, y = NULL) +
    theme_classic(base_size = 9) +
    theme(
      axis.text.x        = element_text(size = 8),
      axis.title.x       = element_text(size = 9, face = "bold"),
      legend.position    = "none",
      panel.grid.major.x = element_line(linewidth = 0.2, colour = "#eee"),
      plot.margin        = margin(5, 5, 5, 5)
    )

  if (show_y) {
    p <- p + theme(axis.text.y = element_text(size = 7.2, colour = "grey15"))
  } else {
    p <- p + theme(axis.text.y  = element_blank(),
                   axis.ticks.y = element_blank())
  }
  p
}

# ---- Brace panel: vertical lines + labels on the very left ------------------
# y positions in the discrete axis (levels reversed: bottom = pan, top = ref):
#   pan_1 ... pan_34 -> y = 1 .. 34
#   spacer_2         -> y = 35
#   arg_2, arg_1     -> y = 36, 37  (MGYG000295553, MGYG000292883)
#   spacer_1         -> y = 38
#   ref_6 ... ref_1  -> y = 39 .. 44

brace_x_line  <- 0.7
brace_x_text  <- 0.4

brace_panel <- ggplot() +
  # Reference brace
  geom_segment(aes(x = brace_x_line, xend = brace_x_line,
                   y = 39, yend = 44),
               linewidth = 0.7, colour = "grey25") +
  geom_segment(aes(x = brace_x_line, xend = brace_x_line + 0.15,
                   y = 39, yend = 39),
               linewidth = 0.7, colour = "grey25") +
  geom_segment(aes(x = brace_x_line, xend = brace_x_line + 0.15,
                   y = 44, yend = 44),
               linewidth = 0.7, colour = "grey25") +
  annotate("text", x = brace_x_text, y = 41.5,
           label = "Reference\n(single)", fontface = "bold",
           size = 3.0, angle = 90, lineheight = 0.9) +
  # ARG-MAG brace
  geom_segment(aes(x = brace_x_line, xend = brace_x_line,
                   y = 36, yend = 37),
               linewidth = 0.7, colour = "#B22222") +
  geom_segment(aes(x = brace_x_line, xend = brace_x_line + 0.15,
                   y = 36, yend = 36),
               linewidth = 0.7, colour = "#B22222") +
  geom_segment(aes(x = brace_x_line, xend = brace_x_line + 0.15,
                   y = 37, yend = 37),
               linewidth = 0.7, colour = "#B22222") +
  annotate("text", x = brace_x_text, y = 36.5,
           label = "ARG-MAG\n(single)", fontface = "bold",
           size = 3.0, angle = 90, colour = "#B22222", lineheight = 0.9) +
  # Pan-Draft brace
  geom_segment(aes(x = brace_x_line, xend = brace_x_line,
                   y = 1, yend = 34),
               linewidth = 0.7, colour = "grey25") +
  geom_segment(aes(x = brace_x_line, xend = brace_x_line + 0.15,
                   y = 1, yend = 1),
               linewidth = 0.7, colour = "grey25") +
  geom_segment(aes(x = brace_x_line, xend = brace_x_line + 0.15,
                   y = 34, yend = 34),
               linewidth = 0.7, colour = "grey25") +
  annotate("text", x = brace_x_text, y = 17.5,
           label = "Pan-Draft\n(species-rep)", fontface = "bold",
           size = 3.0, angle = 90, lineheight = 0.9) +
  scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0.5, 44.5), expand = c(0, 0)) +
  theme_void() +
  theme(plot.margin = margin(5, 0, 5, 5))

# ---- 4 metric panels (equal width) -----------------------------------------
p_a <- make_panel(all_dat, "orf_coverage_pct",    "ORF coverage (%)",        "%.1f", show_y = TRUE)
p_b <- make_panel(all_dat, "n_reactions",         "Reactions (n)",           "%d")
p_c <- make_panel(all_dat, "n_active_substrates", "Active substrates (n)",   "%d")
p_d <- make_panel(all_dat, "n_auxotrophies",      "Auxotrophic factors (n)", "%d")

# ---- Build separate legends -------------------------------------------------
# Pan-Draft class legend
pan_legend_dummy <- ggplot(all_dat %>% filter(role == "pan"),
                            aes(x = 1, y = display, fill = Class)) +
  geom_col() +
  scale_fill_manual(values = pan_class_palette,
                    name = "Pan-Draft (class)",
                    drop = FALSE, na.translate = FALSE) +
  guides(fill = guide_legend(nrow = 1)) +
  theme_void() +
  theme(legend.position  = "bottom",
        legend.key.size  = unit(0.4, "cm"),
        legend.text      = element_text(size = 8.5),
        legend.title     = element_text(size = 9, face = "bold"))

# Single-genome group legend
single_legend_dummy <- ggplot(all_dat %>% filter(role %in% c("ref","argmag")),
                               aes(x = 1, y = display, fill = group)) +
  geom_col() +
  scale_fill_manual(values = single_group_palette,
                    name = "Single genome",
                    drop = FALSE, na.translate = FALSE) +
  guides(fill = guide_legend(nrow = 1)) +
  theme_void() +
  theme(legend.position  = "bottom",
        legend.key.size  = unit(0.4, "cm"),
        legend.text      = element_text(size = 8.5),
        legend.title     = element_text(size = 9, face = "bold"))

pan_legend_grob    <- cowplot::get_legend(pan_legend_dummy)
single_legend_grob <- cowplot::get_legend(single_legend_dummy)

legend_row <- cowplot::plot_grid(single_legend_grob, pan_legend_grob,
                                  ncol = 2, rel_widths = c(1, 2))

# ---- Assemble main figure ---------------------------------------------------
panels <- brace_panel + p_a + p_b + p_c + p_d +
  plot_layout(widths = c(0.14, 1, 1, 1, 1)) +
  plot_annotation(tag_levels = list(c("", "A", "B", "C", "D"))) &
  theme(plot.tag = element_text(face = "bold", size = 13))

final_fig <- panels / patchwork::wrap_elements(legend_row) +
  plot_layout(heights = c(1, 0.04))

# ---- Save ------------------------------------------------------------------
out_path <- file.path(out_dir, "Fig_42_4metrics.png")
ggsave(out_path, final_fig,
       width = 20, height = 10, units = "in",
       dpi = 1000, bg = "white")

cat(sprintf("Saved: %s\n", out_path))
cat(sprintf("Total rows shown: %d (42 organisms + 2 spacers)\n",
            nrow(all_dat)))
cat("\nORF coverage (pan) source: pandraft_all_stage1_summary$orf_cov_mean\n")
cat("ORF coverage range across 42 organisms:\n")
print(summary(all_dat$orf_coverage_pct))
