# =============================================================
# In-vitro phenotype cards
#
# Four test organisms reconstructed via the same pipeline as the
# 42-organism manuscript set. Data are stored in the 'in-vitro'
# sheet as four stacked blocks (reconstruction summary, auxotrophies,
# substrates, model summary), separated by blank rows.
#
# Output: one combined 2x2 figure with four three-panel cards.
#   Fig6_invitro_4cards.png  (forms Panel A of Figure 6; per-strain
#   GEM predictions for the four in-vitro panel models.)
#
# Each card mirrors the template used for Fig 1C and the pair
# figures: aux chips (top) -> top substrates (middle) ->
# predicted secretion (bottom). Shared auxotrophy category legend
# at the figure bottom.
# =============================================================

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
  library(cowplot)
  library(ggtext)        # for italic species names in titles
})

# ---------- Configuration ----------
xlsx_path <- "C:/Users/CFL/OneDrive - Queen's University Belfast/PhD_onedrive/Manuscript4/others/M4_Supplement.xlsx"
out_dir   <- "C:/Users/CFL/OneDrive - Queen's University Belfast/PhD_onedrive/Manuscript4/figures"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# Title colours: one per organism (distinct palette, none colliding with
# manuscript main-set role coding red/gold/green)
test_palette <- c(
  "GCF_000426565.1" = "#4F8FA8",   # muted blue
  "SH_0760"         = "#7A5C9E",   # purple
  "SH_0717"         = "#2E8B6E",   # teal-green
  "SH_0542"         = "#C46A39"    # warm orange-brown
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

# ---------- Compound name shortener (mirror template) ----------
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

# ---------- Read in-vitro sheet, parse 4 blocks --------------------------
# Block 1 (rows 1-4 after header)   = reconstruction summary (incl. top_products)
# Block 2 (rows 8-50 after header)  = auxotrophies
# Block 3 (rows 54-112 after hdr)   = substrates
# Block 4 (rows 116-119 after hdr)  = model summary

recon <- read_excel(xlsx_path, sheet = "in-vitro",
                    skip = 0, n_max = 4)
aux   <- read_excel(xlsx_path, sheet = "in-vitro",
                    skip = 6, n_max = 44) %>%
         filter(!is.na(compound))
sub   <- read_excel(xlsx_path, sheet = "in-vitro",
                    skip = 52, n_max = 60) %>%
         filter(!is.na(compound)) %>%
         mutate(Net_Growth = suppressWarnings(as.numeric(Net_Growth))) %>%
         filter(Net_Growth >= SUB_THRESHOLD)

cat("=== in-vitro data loaded ===\n")
cat(sprintf("  Reconstruction summary: %d organisms\n", nrow(recon)))
cat(sprintf("  Auxotrophies:           %d rows across %d organisms\n",
            nrow(aux), n_distinct(aux$organism)))
cat(sprintf("  Active substrates:      %d rows across %d organisms\n",
            nrow(sub), n_distinct(sub$organism)))

# ---------- Helpers: count real items per organism -----------------------
count_real_substrates <- function(org_id) {
  n <- sub %>% filter(organism == org_id) %>% nrow()
  min(n, SLOT_MAX)
}

count_real_metabolites <- function(org_id) {
  n <- recon %>%
    filter(organism == org_id) %>%
    select(top_products) %>%
    separate_rows(top_products, sep = ",\\s*") %>%
    separate(top_products, into = c("metabolite", "flux"),
             sep = ":", convert = FALSE) %>%
    mutate(flux = suppressWarnings(as.numeric(flux))) %>%
    filter(!is.na(flux), flux >= TRACE_FLUX_CUTOFF) %>%
    nrow()
  min(n, SLOT_MAX)
}

# ---------- Card builder (mirror template) -------------------------------
build_card <- function(org_id, display_label, title_colour,
                       sub_slots, met_slots) {

  # Aux chips
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
        plot.title   = element_markdown(size = 14, face = "bold",
                                        colour = title_colour,
                                        hjust = 0.5,
                                        margin = margin(b = 6))
      )
  }

  # Substrates (pad to sub_slots)
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

  # Secreted metabolites (pad to met_slots)
  met_m <- recon %>%
    filter(organism == org_id) %>%
    select(organism, top_products) %>%
    separate_rows(top_products, sep = ",\\s*") %>%
    separate(top_products, into = c("metabolite", "flux"),
             sep = ":", convert = FALSE) %>%
    mutate(metabolite = trimws(metabolite),
           metabolite = shorten_compound(metabolite),
           flux       = suppressWarnings(as.numeric(flux))) %>%
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

  # Card-level height ratios
  H_AUX <- 1.0
  H_SUB <- 1.4 * sub_slots / SLOT_MAX
  H_MET <- 1.4 * met_slots / SLOT_MAX

  p_aux / p_sub / p_met + plot_layout(heights = c(H_AUX, H_SUB, H_MET))
}

# ---------- Build cards ------------------------------------------------------
organism_cfg <- list(
  list(id = "GCF_000426565.1",
       label = "*Segatella albensis* DSM 11370"),
  list(id = "SH_0760", label = "SH_0760"),
  list(id = "SH_0717", label = "SH_0717"),
  list(id = "SH_0542", label = "SH_0542")
)

# Figure-wide slot sizes (use the max across all 4 to keep proportions even)
sub_slots <- max(sapply(organism_cfg, function(o) count_real_substrates(o$id)))
met_slots <- max(sapply(organism_cfg, function(o) count_real_metabolites(o$id)))

cat(sprintf("\nFigure slot sizes: substrates = %d, metabolites = %d\n",
            sub_slots, met_slots))

cards <- lapply(organism_cfg, function(o) {
  build_card(o$id, o$label, test_palette[[o$id]], sub_slots, met_slots)
})

# Shared aux legend
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

# 2x2 layout: cards 1+2 top row, 3+4 bottom row, legend below
fig <- ((cards[[1]] | cards[[2]]) / (cards[[3]] | cards[[4]])) /
       patchwork::wrap_elements(aux_legend) +
       plot_layout(heights = c(1, 1, 0.04))

out_path <- file.path(out_dir, "Fig6_invitro_4cards.png")
ggsave(out_path, fig,
       width = 18, height = 16, units = "in",
       dpi = 1000, bg = "white")

cat(sprintf("\nSaved: %s\n", out_path))
cat(sprintf("Card slot sizes: substrates = %d, metabolites = %d\n",
            sub_slots, met_slots))
