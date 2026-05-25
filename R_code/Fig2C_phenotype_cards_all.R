# =============================================================
# Metabolic phenotype cards (single script, four figures).
#
#   Fig1C_ARG_MAG_metabolic_phenotype.png : MGYG000292883 + MGYG000295553
#   Fig_pair_Pseudomonadota.png           : E. coli K-12 MG1655 vs E. coli PA-3
#   Fig_pair_Bacteroidota.png             : B. thetaiotaomicron VPI-5482 vs KPPR-3
#   Fig_pair_Bacillota.png                : E. faecalis ATCC 19433 vs E. faecalis 68A
#
# Each card stacks three sub-panels:
#   Top    | Essential growth factors (auxotrophies), as coloured chips
#   Middle | Top preferred substrates, horizontal bars (Net_Growth)
#   Bottom | Predicted secreted metabolites, horizontal bars (flux)
#
# Title colour reflects group:
#   MAG Target (red), PCG Reference (gold), Rumen Reference (green).
#
# Behaviours worth flagging up front:
#   1. shorten_compound() rewrites verbose gapseq names into compact ones.
#   2. Predicted secretion fluxes < TRACE_FLUX_CUTOFF mmol gDW^-1 h^-1 are
#      filtered out as numerical trace by-products of pFBA.
#   3. Substrate and metabolite panels use pair-level dynamic slot counts:
#      both cards in a pair are padded to the larger of the two real
#      counts. Figure height scales proportionally so empty slots are not
#      reserved for a pair that does not need them.
# =============================================================

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
  library(cowplot)
})

# ---------- Configuration ----------
xlsx_path <- "C:/Users/CFL/OneDrive - Queen's University Belfast/PhD_onedrive/Manuscript4/M4_Supplement.xlsx"
out_dir   <- "C:/Users/CFL/OneDrive - Queen's University Belfast/PhD_onedrive/Manuscript4/figures"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# Group title colours
mag_red     <- "#B22222"
pcg_gold    <- "#D4A017"
rumen_green <- "#2E7D32"

# Bar colours
substrate_col <- "#5B9BD5"   # blue, inputs
product_col   <- "#C9622B"   # orange, outputs

# Auxotrophy category palette
aux_pal <- c(
  "Amino acid"       = "#4C9F70",
  "Vitamin/cofactor" = "#9B59B6",
  "Metal"            = "#D4A017",
  "Other"            = "#7F8C8D"
)

# Display constants
SLOT_MAX          <- 10       # cap on rows per panel (i.e. top-N cap)
TRACE_FLUX_CUTOFF <- 0.005    # mmol gDW^-1 h^-1; below this, secretion is treated as trace
REF_FIG_HEIGHT    <- 9.0      # figure height (inches) when both panels use SLOT_MAX rows
REF_FIG_WIDTH     <- 12.0     # figure width (inches), kept fixed across pairs

# Per-panel base height ratios (for SLOT_MAX rows)
H_AUX <- 1.0
H_SUB <- 1.4
H_MET <- 1.4
REF_TOTAL_RATIO <- H_AUX + H_SUB + H_MET

# ---------- Compound name shortener ----------
shorten_compound <- function(x) {
  case_when(
    grepl("^N-Acetyl-beta-D-glucosaminyl-1,6", x) ~ "GlcNAcβ1-6",
    grepl("^starch \\(n=27", x)                    ~ "Starch (n=27)",
    grepl("^starch \\(n=19", x)                    ~ "Starch (n=19)",
    grepl("^Inulin",         x)                    ~ "Inulin",
    x == "Xylooligo-b-1-4"                          ~ "Xylooligo-b-1-4",
    x == "Xylan-b-1-4"                              ~ "Xylan-b-1-4",
    x == "5-Methylthio-D-ribose"                    ~ "5-Methylthio-D-ribose",
    x == "N-acetylneuraminate"                      ~ "N-acetylneuraminate",
    TRUE ~ x
  )
}

# ---------- Auxotrophy classifier ----------
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

# ---------- Figure definitions ----------
figures_cfg <- list(
  list(
    out_name = "Fig3C_ARG_MAG_metabolic_phenotype.png",
    cards = list(
      list(id = "MGYG000292883", label = "MGYG000292883", colour = mag_red),
      list(id = "MGYG000295553", label = "MGYG000295553", colour = mag_red)
    )
  ),
  list(
    out_name = "FigS1_pair_Pseudomonadota.png",
    cards = list(
      list(id = "Ecoli_K12_MG1655", label = "E. coli K-12 MG1655", colour = pcg_gold),
      list(id = "Ecoli_PA3",        label = "E. coli PA-3",        colour = rumen_green)
    )
  ),
  list(
    out_name = "FigS1_pair_Bacteroidota.png",
    cards = list(
      list(id = "Btheta_VPI5482", label = "B. thetaiotaomicron VPI-5482", colour = pcg_gold),
      list(id = "Btheta_KPPR3",   label = "B. thetaiotaomicron KPPR-3",   colour = rumen_green)
    )
  ),
  list(
    out_name = "FigS1_pair_Bacillota.png",
    cards = list(
      list(id = "Efaecalis_ATCC19433", label = "E. faecalis ATCC 19433", colour = pcg_gold),
      list(id = "Efaecalis_68A",       label = "E. faecalis 68A",        colour = rumen_green)
    )
  )
)

# ---------- Read data once ----------
aux   <- read_excel(xlsx_path, sheet = "single_all_auxotrophies")
sub   <- read_excel(xlsx_path, sheet = "single_top15_substrates")
recon <- read_excel(xlsx_path, sheet = "single_all_reconstruction")

# ---------- Helpers: count real items per organism ----------
count_real_substrates <- function(org_id) {
  n <- sub %>%
    filter(organism == org_id) %>%
    select(compound, Net_Growth) %>%
    arrange(desc(Net_Growth)) %>%
    nrow()
  min(n, SLOT_MAX)
}

count_real_metabolites <- function(org_id) {
  n <- recon %>%
    filter(organism == org_id) %>%
    select(top_products) %>%
    separate_rows(top_products, sep = ",\\s*") %>%
    separate(top_products, into = c("metabolite", "flux"),
             sep = ":", convert = FALSE) %>%
    mutate(flux = as.numeric(flux)) %>%
    filter(!is.na(flux), flux >= TRACE_FLUX_CUTOFF) %>%
    nrow()
  min(n, SLOT_MAX)
}

# ---------- Card builder ----------
build_card <- function(org_id, display_label, title_colour,
                       sub_slots, met_slots) {
  
  # ----- Top: auxotrophy chips -----
  aux_m <- aux %>%
    filter(organism == org_id) %>%
    select(organism, compound) %>%
    mutate(compound = shorten_compound(compound),
           category = factor(aux_class(compound), levels = names(aux_pal))) %>%
    arrange(category, compound) %>%
    group_by(category) %>%
    mutate(col_idx = row_number()) %>%
    ungroup() %>%
    mutate(category = factor(category, levels = rev(names(aux_pal))))
  
  if (nrow(aux_m) == 0) {
    p_aux <- ggplot() +
      annotate("text", x = 0, y = 0,
               label = paste0(display_label, "  |  no auxotrophies detected"),
               size = 4, fontface = "bold", colour = title_colour) +
      theme_void()
  } else {
    max_col <- max(aux_m$col_idx)
    p_aux <- ggplot(aux_m,
                    aes(x = col_idx, y = category,
                        fill = category, label = compound)) +
      geom_label(colour = "white", size = 2.4, fontface = "bold",
                 label.padding = unit(0.18, "lines"),
                 label.r       = unit(0.10, "lines")) +
      scale_fill_manual(values = aux_pal, drop = FALSE) +
      scale_x_continuous(limits = c(0.4, max_col + 0.6), expand = c(0, 0)) +
      scale_y_discrete(drop = FALSE) +
      labs(title = paste0(display_label,
                          "  |  Essential growth factors (n = ",
                          nrow(aux_m), ")"),
           x = NULL, y = NULL) +
      theme_classic(base_size = 10) +
      theme(
        axis.text.x  = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line.x  = element_blank(),
        axis.text.y  = element_text(size = 9, colour = "grey20"),
        axis.ticks.y = element_blank(),
        axis.line.y  = element_blank(),
        legend.position = "none",
        plot.title   = element_text(size = 11, face = "bold",
                                    colour = title_colour,
                                    margin = margin(b = 6))
      )
  }
  
  # ----- Middle: top preferred substrates (pad to sub_slots) -----
  sub_m <- sub %>%
    filter(organism == org_id) %>%
    select(organism, compound, Net_Growth) %>%
    mutate(compound = shorten_compound(compound)) %>%
    arrange(desc(Net_Growth)) %>%
    slice_head(n = sub_slots)
  
  n_sub_real <- nrow(sub_m)
  if (n_sub_real < sub_slots) {
    sub_m <- bind_rows(
      sub_m,
      tibble(organism = org_id,
             compound = paste0("__pad_", seq_len(sub_slots - n_sub_real)),
             Net_Growth = NA_real_)
    )
  }
  sub_m <- sub_m %>%
    mutate(compound = factor(compound, levels = rev(compound)))
  
  p_sub <- ggplot(sub_m, aes(x = Net_Growth, y = compound)) +
    geom_col(fill = substrate_col, colour = "grey20",
             linewidth = 0.3, width = 0.7, na.rm = TRUE) +
    geom_text(data = . %>% filter(!is.na(Net_Growth)),
              aes(label = sprintf("%.3f", Net_Growth)),
              hjust = -0.15, size = 2.8, colour = "grey20") +
    scale_x_continuous(expand = expansion(mult = c(0, 0.22))) +
    scale_y_discrete(labels = function(x) ifelse(grepl("^__pad_", x), "", x),
                     drop = FALSE) +
    labs(title = paste0("Preferred substrates (top ", n_sub_real, ")"),
         x = expression(paste("Net growth (h"^-1, ")")),
         y = NULL) +
    theme_classic(base_size = 10) +
    theme(
      axis.text.y = element_text(size = 8, colour = "grey20"),
      plot.title  = element_text(size = 10, face = "bold",
                                 margin = margin(b = 4))
    )
  
  # ----- Bottom: predicted secreted metabolites (pad to met_slots) -----
  met_m <- recon %>%
    filter(organism == org_id) %>%
    select(organism, top_products) %>%
    separate_rows(top_products, sep = ",\\s*") %>%
    separate(top_products, into = c("metabolite", "flux"),
             sep = ":", convert = FALSE) %>%
    mutate(metabolite = trimws(metabolite),
           metabolite = shorten_compound(metabolite),
           flux       = as.numeric(flux)) %>%
    filter(!is.na(flux), flux >= TRACE_FLUX_CUTOFF) %>%
    arrange(desc(flux)) %>%
    slice_head(n = met_slots)
  
  n_met_real <- nrow(met_m)
  if (n_met_real < met_slots) {
    met_m <- bind_rows(
      met_m,
      tibble(organism = org_id,
             metabolite = paste0("__pad_", seq_len(met_slots - n_met_real)),
             flux = NA_real_)
    )
  }
  met_m <- met_m %>%
    mutate(metabolite = factor(metabolite, levels = rev(metabolite)))
  
  p_met <- ggplot(met_m, aes(x = flux, y = metabolite)) +
    geom_col(fill = product_col, colour = "grey20",
             linewidth = 0.3, width = 0.7, na.rm = TRUE) +
    geom_text(data = . %>% filter(!is.na(flux)),
              aes(label = sprintf("%.2f", flux)),
              hjust = -0.15, size = 2.8, colour = "grey20") +
    scale_x_continuous(expand = expansion(mult = c(0, 0.20))) +
    scale_y_discrete(labels = function(x) ifelse(grepl("^__pad_", x), "", x),
                     drop = FALSE) +
    labs(title = paste0("Predicted secretion (top ", n_met_real, ")"),
         x = expression(paste("Flux (mmol gDW"^-1, " h"^-1, ")")),
         y = NULL) +
    theme_classic(base_size = 10) +
    theme(
      axis.text.y = element_text(size = 8, colour = "grey20"),
      plot.title  = element_text(size = 10, face = "bold",
                                 margin = margin(b = 4))
    )
  
  # Card-level layout: heights scale with slot counts so bar thickness
  # stays consistent across pairs.
  card_heights <- c(H_AUX,
                    H_SUB * sub_slots / SLOT_MAX,
                    H_MET * met_slots / SLOT_MAX)
  
  p_aux / p_sub / p_met + plot_layout(heights = card_heights)
}

# ---------- Shared aux category legend ----------
aux_legend <- cowplot::get_legend(
  ggplot(data.frame(cat = factor(names(aux_pal), levels = names(aux_pal)),
                    x = 1, y = 1),
         aes(x = x, y = y, fill = cat)) +
    geom_col() +
    scale_fill_manual(values = aux_pal, name = NULL, drop = FALSE) +
    theme_void() +
    theme(legend.position = "bottom",
          legend.key.size = unit(0.4, "cm"),
          legend.text     = element_text(size = 9))
)

# ---------- Loop over figures, build and save ----------
for (cfg in figures_cfg) {
  id_a <- cfg$cards[[1]]$id
  id_b <- cfg$cards[[2]]$id
  
  # Pair-level slot counts (max across the two cards)
  sub_slots <- max(count_real_substrates(id_a),  count_real_substrates(id_b))
  met_slots <- max(count_real_metabolites(id_a), count_real_metabolites(id_b))
  
  # Build cards
  card_a <- build_card(id_a, cfg$cards[[1]]$label, cfg$cards[[1]]$colour,
                       sub_slots, met_slots)
  card_b <- build_card(id_b, cfg$cards[[2]]$label, cfg$cards[[2]]$colour,
                       sub_slots, met_slots)
  
  # Combine cards plus shared legend
  fig <- (card_a | card_b) /
         patchwork::wrap_elements(aux_legend) +
         plot_layout(heights = c(1, 0.04))
  
  # Scale figure height proportionally to content
  content_ratio <- (H_AUX +
                    H_SUB * sub_slots / SLOT_MAX +
                    H_MET * met_slots / SLOT_MAX) / REF_TOTAL_RATIO
  fig_h <- REF_FIG_HEIGHT * content_ratio
  
  out_path <- file.path(out_dir, cfg$out_name)
  ggsave(out_path, fig,
         width = REF_FIG_WIDTH, height = fig_h, units = "in",
         dpi = 1000, bg = "white")
  
  cat(sprintf("Saved: %s   (sub_slots=%d, met_slots=%d, height=%.2f in)\n",
              out_path, sub_slots, met_slots, fig_h))
}

