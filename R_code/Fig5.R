# =============================================================
# Fig5.R  (originally Fig_SH_cards_and_family_comparison.R)
# Manuscript 4 — composite figure combining, on one page:
#   Panel A (top)    : 3-card in-vitro phenotype summary for the
#                      three SH-panel strains (SH_0542, SH_0717,
#                      SH_0760). Design mirrors
#                      Fig6_invitro_3cards_1row.R.
#   Panel B (bottom) : SH_0717 vs Ruminococcaceae only — predicted
#                      essential factors (I) and top preferred
#                      substrates (II). SH_0542 and SH_0760 now use a
#                      matched-SGB one-vs-one comparison (Supplementary
#                      Table S11e), so their family-context panels have
#                      been removed from this figure.
#
# Output: figures/Fig5.png
# =============================================================

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
  library(cowplot)
  library(ggtext)
  library(stringr)
  library(forcats)
})

# ---------- Paths ----------
base_dir  <- "C:/Users/CFL/OneDrive - Queen's University Belfast/PhD_onedrive/Manuscript4_gapseq"
xlsx_card <- file.path(base_dir, "others/M4_Supplement.xlsx")      # in-vitro card source
xlsx_main <- file.path(base_dir, "M4_Supplementary_Tables.xlsx")   # S10/S11 sheets
map_tsv   <- file.path(base_dir, "pan-draft-SH/sheet1_sgb_to_mag.tsv")
out_dir   <- file.path(base_dir, "figures")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# =============================================================
# ------------------- PANEL A: 3 cards ------------------------
# =============================================================
test_palette <- c(
  "SH_0760" = "#7A5C9E",
  "SH_0717" = "#2E8B6E",
  "SH_0542" = "#C46A39"
)
substrate_col <- "#5B9BD5"
product_col   <- "#C9622B"
aux_pal <- c(
  "Amino acid"       = "#4C9F70",
  "Vitamin/cofactor" = "#9B59B6",
  "Metal"            = "#D4A017",
  "Other"            = "#7F8C8D"
)
SLOT_MAX          <- 10
TRACE_FLUX_CUTOFF <- 0.005
SUB_THRESHOLD     <- 0.01
TITLE_SIZE_PT     <- 18   # was 13; increased for stacked-card legibility

shorten_compound <- function(x) {
  case_when(
    grepl("^N-Acetyl-beta-D-glucosaminyl-1,6", x) ~ "GlcNAc\u03b21-6",
    grepl("^starch \\(n=27", x)                    ~ "Starch (n=27)",
    grepl("^starch \\(n=19", x)                    ~ "Starch (n=19)",
    grepl("^Inulin",         x)                    ~ "Inulin",
    x == "N-Acetyl-D-glucosamine"                  ~ "GlcNAc",
    x == "N-acetylneuraminate"                     ~ "Neu5Ac",
    x == "D-Glucosamine"                           ~ "D-Glucosamine",
    x == "D-Cellobiose"                            ~ "D-Cellobiose",
    x == "Maltoheptaose"                           ~ "Maltoheptaose",
    TRUE ~ x
  )
}

aux_class <- function(x) {
  case_when(
    x %in% c("L-Glutamate","L-Methionine","L-Cysteine","L-Aspartate",
             "L-Isoleucine","L-Leucine","L-Threonine","L-Valine",
             "L-Histidine","L-Phenylalanine","L-Tyrosine","L-Lysine",
             "L-Tryptophan","L-Arginine","L-Proline","L-Serine",
             "L-Alanine","Glycine","L-Asparagine","L-Glutamine") ~ "Amino acid",
    x %in% c("Heme","Folate","Thiamin","Pyridoxal","Riboflavin",
             "Pantothenic acid","Biotin","Cobalamin","Vitamin B12",
             "Niacin","Pyridoxol","NAD","NADP","FAD") ~ "Vitamin/cofactor",
    x %in% c("Zn2+","Cobalt","Fe3+","Ca2+","Cu2+","Fe2+","Mn2+",
             "Mg","Mg2+","Ni2+","Nickel","K+","Mo","Molybdate") ~ "Metal",
    TRUE ~ "Other"
  )
}

recon <- read_excel(xlsx_card, sheet = "in-vitro", skip = 0, n_max = 4)
aux   <- read_excel(xlsx_card, sheet = "in-vitro", skip = 6, n_max = 44) |>
         filter(!is.na(compound)) |>
         filter(organism %in% names(test_palette))
sub   <- read_excel(xlsx_card, sheet = "in-vitro", skip = 52, n_max = 60) |>
         filter(!is.na(compound)) |>
         mutate(Net_Growth = suppressWarnings(as.numeric(Net_Growth))) |>
         filter(Net_Growth >= SUB_THRESHOLD) |>
         filter(organism %in% names(test_palette))

count_real_substrates <- function(org_id) {
  n <- sub |> filter(organism == org_id) |> nrow()
  min(n, SLOT_MAX)
}
count_real_metabolites <- function(org_id) {
  n <- recon |>
    filter(organism == org_id) |>
    select(top_products) |>
    separate_rows(top_products, sep = ",\\s*") |>
    separate(top_products, into = c("metabolite", "flux"),
             sep = ":", convert = FALSE) |>
    mutate(flux = suppressWarnings(as.numeric(flux))) |>
    filter(!is.na(flux), flux >= TRACE_FLUX_CUTOFF) |>
    nrow()
  min(n, SLOT_MAX)
}

build_card <- function(org_id, display_label, title_colour,
                       sub_slots, met_slots) {

  aux_m <- aux |>
    filter(organism == org_id) |>
    select(organism, compound) |>
    mutate(compound = shorten_compound(compound),
           category = factor(aux_class(compound), levels = names(aux_pal))) |>
    arrange(category, compound) |>
    group_by(category) |>
    mutate(col_idx = row_number()) |>
    ungroup() |>
    mutate(category = factor(category, levels = rev(names(aux_pal))))

  if (nrow(aux_m) == 0) {
    p_aux <- ggplot() +
      annotate("text", x = 0, y = 0,
               label = "Essential factors (n = 0)",
               size = 5, fontface = "bold", colour = "grey20") +
      theme_void() +
      theme(plot.title = element_text(size = 14, face = "bold",
                                      margin = margin(b = 6)),
            plot.margin = margin(t = 4, r = 4, b = 4, l = 4))
  } else {
    max_col <- max(aux_m$col_idx)
    p_aux <- ggplot(aux_m,
                    aes(x = col_idx, y = category,
                        fill = category, label = compound)) +
      geom_label(colour = "white", size = 4.0, fontface = "bold",
                 label.padding = unit(0.28, "lines"),
                 label.r       = unit(0.12, "lines")) +
      scale_fill_manual(values = aux_pal, drop = FALSE) +
      scale_x_continuous(limits = c(0.4, max_col + 0.6), expand = c(0, 0)) +
      scale_y_discrete(drop = FALSE) +
      labs(title = paste0("Essential factors (n = ", nrow(aux_m), ")"),
           x = NULL, y = NULL) +
      theme_classic(base_size = 14) +
      theme(
        axis.text.x  = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line.x  = element_blank(),
        axis.text.y  = element_text(size = 13, colour = "grey20"),
        axis.ticks.y = element_blank(),
        axis.line.y  = element_blank(),
        legend.position = "none",
        plot.title   = element_text(size = 14, face = "bold",
                                    margin = margin(b = 6)),
        plot.margin  = margin(t = 4, r = 4, b = 4, l = 4)
      )
  }

  sub_m <- sub |>
    filter(organism == org_id) |>
    select(organism, compound, Net_Growth) |>
    mutate(compound = shorten_compound(compound)) |>
    arrange(desc(Net_Growth)) |>
    slice_head(n = sub_slots)
  n_sub_real <- nrow(sub_m)
  if (n_sub_real < sub_slots) {
    sub_m <- bind_rows(
      sub_m,
      tibble(organism = org_id,
             compound = paste0("__pad_", seq_len(sub_slots - n_sub_real)),
             Net_Growth = NA_real_))
  }
  sub_m <- sub_m |> mutate(compound = factor(compound, levels = rev(compound)))

  p_sub <- ggplot(sub_m, aes(x = Net_Growth, y = compound)) +
    geom_col(fill = substrate_col, colour = "grey20",
             linewidth = 0.4, width = 0.7, na.rm = TRUE) +
    geom_text(data = . %>% filter(!is.na(Net_Growth)),
              aes(label = sprintf("%.3f", Net_Growth)),
              hjust = -0.15, size = 4.0, colour = "grey20") +  # was 2.7
    scale_x_continuous(expand = expansion(mult = c(0, 0.22))) +
    scale_y_discrete(labels = function(x) ifelse(grepl("^__pad_", x), "", x),
                     drop = FALSE) +
    labs(title = paste0("Preferred substrates (top ", n_sub_real, ")"),
         x = expression(paste("Net growth (h"^-1, ")")),
         y = NULL) +
    theme_classic(base_size = 14) +  # was 10
    theme(axis.text.y = element_text(size = 12, colour = "grey20"),  # was 8.5
          axis.text.x = element_text(size = 11),
          axis.title.x = element_text(size = 12),
          plot.title  = element_text(size = 14, face = "bold",         # was 10.5
                                     margin = margin(b = 6)))

  met_m <- recon |>
    filter(organism == org_id) |>
    select(organism, top_products) |>
    separate_rows(top_products, sep = ",\\s*") |>
    separate(top_products, into = c("metabolite", "flux"),
             sep = ":", convert = FALSE) |>
    mutate(metabolite = trimws(metabolite),
           metabolite = shorten_compound(metabolite),
           flux       = suppressWarnings(as.numeric(flux))) |>
    filter(!is.na(flux), flux >= TRACE_FLUX_CUTOFF) |>
    arrange(desc(flux)) |>
    slice_head(n = met_slots)
  n_met_real <- nrow(met_m)
  if (n_met_real < met_slots) {
    met_m <- bind_rows(
      met_m,
      tibble(organism = org_id,
             metabolite = paste0("__pad_", seq_len(met_slots - n_met_real)),
             flux = NA_real_))
  }
  met_m <- met_m |> mutate(metabolite = factor(metabolite, levels = rev(metabolite)))

  p_met <- ggplot(met_m, aes(x = flux, y = metabolite)) +
    geom_col(fill = product_col, colour = "grey20",
             linewidth = 0.4, width = 0.7, na.rm = TRUE) +
    geom_text(data = . %>% filter(!is.na(flux)),
              aes(label = sprintf("%.2f", flux)),
              hjust = -0.15, size = 4.0, colour = "grey20") +  # was 2.7
    scale_x_continuous(expand = expansion(mult = c(0, 0.20))) +
    scale_y_discrete(labels = function(x) ifelse(grepl("^__pad_", x), "", x),
                     drop = FALSE) +
    labs(title = paste0("Predicted secretion (top ", n_met_real, ")"),
         x = expression(paste("Flux (mmol gDW"^-1, " h"^-1, ")")),
         y = NULL) +
    theme_classic(base_size = 14) +  # was 10
    theme(axis.text.y = element_text(size = 12, colour = "grey20"),  # was 8.5
          axis.text.x = element_text(size = 11),
          axis.title.x = element_text(size = 12),
          plot.title  = element_text(size = 14, face = "bold",         # was 10.5
                                     margin = margin(b = 6)))

  H_AUX <- 1.0
  H_SUB <- 1.4 * sub_slots / SLOT_MAX
  H_MET <- 1.4 * met_slots / SLOT_MAX

  # Card-level main title (strain-coloured "SH_XXXX (Family)") is rendered as
  # its own subplot at the top of each card so it can't be clipped by outer
  # composite labels (previous plot_annotation attempts were hidden by the
  # cowplot "A" label). The three subplots below carry the plain subtitles
  # ("Essential factors (n = X)", "Preferred substrates (top X)",
  # "Predicted secretion (top X)").
  p_title <- ggplot() +
    labs(title = display_label) +
    theme_void() +
    theme(plot.title = element_markdown(size   = TITLE_SIZE_PT,
                                        face   = "bold",
                                        colour = title_colour,
                                        hjust  = 0.5,
                                        margin = margin(t = 2, b = 4)),
          plot.margin = margin(t = 2, r = 4, b = 2, l = 4))

  H_TITLE <- 0.12  # ~12% of the vertical budget reserved for the strain title

  (p_title / p_aux / p_sub / p_met +
     plot_layout(heights = c(H_TITLE, H_AUX, H_SUB, H_MET)))
}

organism_cfg <- list(
  list(id = "SH_0542", label = "SH_0542 (Coprobacillaceae)"),
  list(id = "SH_0717", label = "SH_0717 (Ruminococcaceae)"),
  list(id = "SH_0760", label = "SH_0760 (Erysipelotrichaceae)")
)
sub_slots <- max(sapply(organism_cfg, function(o) count_real_substrates(o$id)))
met_slots <- max(sapply(organism_cfg, function(o) count_real_metabolites(o$id)))
cards <- lapply(organism_cfg, function(o) {
  build_card(o$id, o$label, test_palette[[o$id]], sub_slots, met_slots)
})

aux_legend <- cowplot::get_legend(
  ggplot(data.frame(cat = factor(names(aux_pal), levels = names(aux_pal)),
                    x = 1, y = 1),
         aes(x = x, y = y, fill = cat)) +
    geom_col() +
    scale_fill_manual(values = aux_pal, name = NULL, drop = FALSE) +
    theme_void() +
    theme(legend.position = "bottom",
          legend.key.size = unit(0.6, "cm"),           # was 0.4
          legend.text     = element_text(size = 13))    # was 9
)

# Panel A: three cards laid out HORIZONTALLY (side-by-side). Font sizes have
# been increased across the card elements (title, axis text, labels, geom_text)
# to improve legibility at journal size; overall figure width is enlarged
# below to accommodate the larger fonts.
# Each card carries its own strain title as its top subplot (see build_card()),
# so no extra top margin is needed here to protect it from the composite "A"
# label.
panel_A <- ((cards[[1]] | cards[[2]] | cards[[3]]) &
             theme(plot.margin = margin(t = 4, r = 4, b = 4, l = 4))) /
           patchwork::wrap_elements(aux_legend) +
           plot_layout(heights = c(1, 0.05))

# =============================================================
# ------------- PANEL B: SH vs family (I / II / III / IV) -----
# =============================================================
map_df <- read.delim(map_tsv, sep = "\t", stringsAsFactors = FALSE)
org2fam <- unique(rbind(
  data.frame(organism = map_df$sgb_id, family = map_df$family, stringsAsFactors = FALSE),
  data.frame(organism = map_df$mag_id, family = map_df$family, stringsAsFactors = FALSE)
))
org2fam <- org2fam[!duplicated(org2fam$organism), ]

sh_ef   <- read_excel(xlsx_main, sheet = "S10b_Essential_factors",         skip = 3)
sh_sub  <- read_excel(xlsx_main, sheet = "S10c_All_substrates",            skip = 3)
# NB: Section 1 revision moved matched-SGB comparison into S11a, so the
# family sheets shifted one letter: family_ess_factors is now S11c and
# family_substrates_top15 is now S11d (see xlsx TOC row 15).
fam_ef  <- read_excel(xlsx_main, sheet = "S11c_SH_family_ess_factors",     skip = 3)
fam_sub <- read_excel(xlsx_main, sheet = "S11d_SH_family_substrates_top15", skip = 3)
fam_ef  <- fam_ef  %>% left_join(org2fam, by = c("organism" = "organism"))
fam_sub <- fam_sub %>% left_join(org2fam, by = c("organism" = "organism"))

shorten <- function(x) {
  x <- str_replace_all(x, "starch \\(n=27, .*\\)",  "Starch (n=27)")
  x <- str_replace_all(x, "starch \\(n=19, .*\\)",  "Starch (n=19)")
  x <- str_replace_all(x, "^N-Acetyl-beta-D-glucosaminyl-1,6.*", "GlcNAc\u03b21-6 (glycan)")
  x <- str_replace_all(x, "^N-Acetyl-D-glucosamine",             "GlcNAc")
  ifelse(nchar(x) > 34, paste0(substr(x, 1, 32), "\u2026"), x)
}

build_panel_df <- function(sh_df, fam_df, strain_id, family_name,
                           top_sh_show = 15, top_family_only = 5) {
  sh_items <- sh_df %>%
    filter(strain_id == !!strain_id) %>%
    pull(compound) %>% unique()

  fam_this <- fam_df %>% filter(family == !!family_name)
  fam_sgbs <- length(unique(fam_this$organism))
  if (fam_sgbs == 0) stop("No family SGBs found for ", family_name)

  prev <- fam_this %>%
    distinct(organism, compound) %>%
    count(compound) %>%
    mutate(prev_pct = 100 * n / fam_sgbs) %>%
    select(compound, prev_pct, n_present = n)

  sh_df_panel <- data.frame(compound = sh_items, has_sh = TRUE, stringsAsFactors = FALSE) %>%
    left_join(prev, by = "compound") %>%
    mutate(prev_pct  = ifelse(is.na(prev_pct), 0, prev_pct),
           n_present = ifelse(is.na(n_present), 0, n_present)) %>%
    arrange(desc(prev_pct)) %>%
    head(top_sh_show)

  fam_only <- prev %>%
    filter(!(compound %in% sh_items)) %>%
    arrange(desc(prev_pct)) %>%
    head(top_family_only) %>%
    mutate(has_sh = FALSE)

  bind_rows(sh_df_panel, fam_only) %>%
    mutate(compound_lbl = shorten(compound),
           group      = ifelse(has_sh, "In SH strain", "Family-common, absent in SH"),
           n_fam_sgbs = fam_sgbs)
}

d_A <- build_panel_df(sh_ef,  fam_ef,  "SH_0717", "Ruminococcaceae")
d_B <- build_panel_df(sh_sub, fam_sub, "SH_0717", "Ruminococcaceae")
# Section 1 revision: SH_0760 was previously compared against the whole
# Erysipelotrichaceae family (panels III / IV). It now collapses to a
# matched-SGB one-vs-one comparison against MGYG000290613 (singleton),
# reported quantitatively in Supplementary Table S11e rather than as a
# figure panel. SH_0542 similarly collapses to MGYG000290992 (S11e).
# Panels III / IV are therefore no longer built.

grp_pal <- c("In SH strain" = "#2E8B6E", "Family-common, absent in SH" = "#B0B0B0")

make_panel <- function(df, subtitle, palette, x_max = 100, base_size = 15,
                       show_x_lab = TRUE) {
  df <- df %>%
    group_by(group) %>%
    mutate(order_in_grp = row_number()) %>%
    ungroup() %>%
    mutate(
      global_order = ifelse(group == "In SH strain", order_in_grp,
                            max(order_in_grp[group == "In SH strain"]) + order_in_grp),
      compound_f = factor(compound_lbl,
                          levels = rev(unique(compound_lbl[order(global_order)])))
    )
  n_fam <- df$n_fam_sgbs[1]
  ggplot(df, aes(x = prev_pct, y = compound_f, fill = group)) +
    geom_col(width = 0.72) +
    geom_text(data = subset(df, prev_pct >= 15),
              aes(label = sprintf("%.0f%%", prev_pct)),
              hjust = 1.10, size = 4.4, colour = "white", fontface = "bold") +
    geom_text(data = subset(df, prev_pct <  15),
              aes(label = sprintf("%.0f%%", prev_pct)),
              hjust = -0.15, size = 4.4, colour = "grey20") +
    scale_x_continuous(limits = c(0, x_max + 2), breaks = c(0, 25, 50, 75, 100),
                       labels = function(x) paste0(x, "%"), expand = c(0, 0)) +
    scale_y_discrete(drop = FALSE) +
    scale_fill_manual(values = palette, name = NULL, drop = FALSE) +
    labs(subtitle = subtitle,
         x = if (show_x_lab) sprintf("Prevalence in %d family SGBs", n_fam) else NULL,
         y = NULL) +
    theme_minimal(base_size = base_size) +
    theme(
      plot.subtitle       = element_text(face = "bold", size = base_size + 2.5,
                                         margin = margin(b = 6)),
      axis.text.y         = element_text(size = base_size + 0.5, colour = "grey15"),
      axis.text.x         = element_text(size = base_size - 1),
      axis.title.x        = element_text(size = base_size, colour = "grey25",
                                         margin = margin(t = 4)),
      panel.grid.major.y  = element_blank(),
      panel.grid.minor    = element_blank(),
      panel.grid.major.x  = element_line(colour = "grey92"),
      legend.position     = "bottom",
      legend.text         = element_text(size = base_size + 2),
      legend.key.size     = unit(0.9, "cm"),
      plot.margin         = margin(2, 22, 2, 22)  # wider L/R padding so I|II and III|IV columns sit further apart
    )
}

# Two-line subtitles with inline I / II tags — tags sit on the same
# horizontal baseline as the strain-vs-family text.
p_I  <- make_panel(d_A,
                   "I    SH_0717 vs Ruminococcaceae\nPredicted essential factors",
                   grp_pal)
p_II <- make_panel(d_B,
                   "II   SH_0717 vs Ruminococcaceae\nTop preferred substrates",
                   grp_pal)

# Panel B is now a single row of two side-by-side plots (SH_0717 only);
# SH_0760 and SH_0542 matched-SGB comparisons are in Supplementary Table S11e.
panel_B <- (p_I | p_II) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom",
        legend.text     = element_text(size = 17),
        legend.key.size = unit(0.9, "cm"))

# =============================================================
# ------------- Composite A over B ----------------------------
# =============================================================
# Use cowplot::plot_grid for the outer A/B labels so panel_A's card layout
# and panel_B's inner I/II/III/IV tags both survive. label_x pushed hard
# to the left so the outer tags do not overlap with any panel content.
composite <- cowplot::plot_grid(
  panel_A, panel_B,
  labels        = c("A", "B"),
  label_size    = 34,
  label_fontface = "bold",
  ncol          = 1,
  rel_heights   = c(1, 0.7),  # standard 3 cards side-by-side; top padding is
                              # added inside panel_A above (not by enlarging
                              # this ratio, which would shrink the fonts)
  hjust         = 0,
  vjust         = 1,
  label_x       = 0.002,
  label_y       = 0.995
)

# Output dimensions: 26 × 18 (landscape) — matches the font scale set for
# card content; extra top padding inside panel_A protects the strain titles
# from being clipped by the composite "A" label without stretching the whole
# figure taller (which would visually shrink the fonts).
out_png <- file.path(out_dir, "Fig5.png")
ggsave(out_png, composite, width = 26, height = 18, dpi = 400, bg = "white")
cat("Wrote:", out_png, "\n")
