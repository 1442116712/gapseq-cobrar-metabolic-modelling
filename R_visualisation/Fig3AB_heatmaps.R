# =============================================================
# Figure 3 - panels A and B (heatmaps)
#
# Panel A: Auxotrophy heatmap. Cells are coloured by category.
#
# Panel B: Substrate preference heatmap.
#   - Reads single_all_substrates and pandraft_all_substrates and
#     applies the active-substrate threshold Net_Growth >= 0.01 h^-1
#     (gapseq default). Below-threshold rows in the supplementary
#     sheets are exchange-flux artefacts (CO2, urea, methane, fatty
#     acids etc.) and are excluded from the heatmap.
#   - 117 surviving compounds are mapped onto 8 manual categories
#     (see Supplementary Table SX), shown as a coloured strip above
#     the heatmap. Empty columns separate adjacent categories.
#   - Within each category, compounds are ordered by the number
#     of organisms with a hit (most to least).
#   - Cells use a pseudo_log gradient (sigma = 0.1) on Net_Growth
#     so small values stay visible alongside large ones.
#
# Y-axis order (top to bottom):
#   References (Gammaproteobacteria, Bacilli, Bacteroidia)
#   Bacteroidia pan models
#   Clostridia pan models
#   ARG-MAGs (Verrucomicrobiota)
#
# Output: Fig3_AB_heatmaps.png
# =============================================================

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(patchwork)
  library(scales)
})

# ---------- Configuration ----------
xlsx_path <- "C:/Users/CFL/OneDrive - Queen's University Belfast/PhD_onedrive/Manuscript4/M4_Supplement.xlsx"
out_dir   <- "C:/Users/CFL/OneDrive - Queen's University Belfast/PhD_onedrive/Manuscript4/figures"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

LOG_SIGMA     <- 0.1    # pseudo-log compression for Panel B fill
SUB_THRESHOLD <- 0.01   # active-substrate threshold (gapseq default, h^-1)

# ---------- Palettes ----------
aux_pal <- c(
  "Amino acid"       = "#4C9F70",
  "Vitamin/cofactor" = "#9B59B6",
  "Metal"            = "#D4A017",
  "Other"            = "#7F8C8D"
)

cat_pal <- c(
  "Amino acids"               = "#88BBA1",
  "Mono- & disaccharides"     = "#E8B080",
  "Amino sugars"              = "#D88FA8",
  "Oligo- & polysaccharides"  = "#C68272",
  "Sugar alcohols & alcohols" = "#B59ECC",
  "Organic & sugar acids"     = "#7FA5BE",
  "Nucleosides & bases"       = "#80BAA5",
  "Other"                     = "#B5BABC"
)
cat_levels <- names(cat_pal)

# ---------- Compound -> category dictionary ----------
# Sourced from Supplementary Table SX. Mapping derived from KEGG BRITE
# br08001 + ModelSEED Biochemistry, with two adaptations: (i) carbohydrates
# split by degree of polymerisation; (ii) cofactor / aldehyde / bile acid
# / mucin-glycan compounds collected into "Other".
compound_category <- c(
  # Amino acids (19): proteinogenic + Ornithine + GABA + D-isomers
  cpd00023 = "Amino acids", cpd00033 = "Amino acids",
  cpd00035 = "Amino acids", cpd00039 = "Amino acids",
  cpd00041 = "Amino acids", cpd00051 = "Amino acids",
  cpd00053 = "Amino acids", cpd00054 = "Amino acids",
  cpd00060 = "Amino acids", cpd00064 = "Amino acids",
  cpd00065 = "Amino acids", cpd00084 = "Amino acids",
  cpd00117 = "Amino acids", cpd00129 = "Amino acids",
  cpd00132 = "Amino acids", cpd00156 = "Amino acids",
  cpd00161 = "Amino acids", cpd00281 = "Amino acids",
  cpd00550 = "Amino acids",
  # Mono- & disaccharides (17): hexoses, pentoses, deoxy-sugars + DP=2
  cpd00027 = "Mono- & disaccharides", cpd00076 = "Mono- & disaccharides",
  cpd00082 = "Mono- & disaccharides", cpd00105 = "Mono- & disaccharides",
  cpd00108 = "Mono- & disaccharides", cpd00138 = "Mono- & disaccharides",
  cpd00154 = "Mono- & disaccharides", cpd00158 = "Mono- & disaccharides",
  cpd00179 = "Mono- & disaccharides", cpd00185 = "Mono- & disaccharides",
  cpd00208 = "Mono- & disaccharides", cpd00224 = "Mono- & disaccharides",
  cpd00396 = "Mono- & disaccharides", cpd00751 = "Mono- & disaccharides",
  cpd00794 = "Mono- & disaccharides", cpd01354 = "Mono- & disaccharides",
  cpd03198 = "Mono- & disaccharides",
  # Amino sugars (4)
  cpd00122 = "Amino sugars", cpd00232 = "Amino sugars",
  cpd00276 = "Amino sugars", cpd00832 = "Amino sugars",
  # Oligo- & polysaccharides (17): DP >= 3
  cpd00155 = "Oligo- & polysaccharides", cpd00382 = "Oligo- & polysaccharides",
  cpd01133 = "Oligo- & polysaccharides", cpd01262 = "Oligo- & polysaccharides",
  cpd01399 = "Oligo- & polysaccharides", cpd11976 = "Oligo- & polysaccharides",
  cpd15494 = "Oligo- & polysaccharides", cpd28763 = "Oligo- & polysaccharides",
  cpd90003 = "Oligo- & polysaccharides", cpd90004 = "Oligo- & polysaccharides",
  cpd90005 = "Oligo- & polysaccharides", cpd90006 = "Oligo- & polysaccharides",
  cpd90007 = "Oligo- & polysaccharides", cpd90008 = "Oligo- & polysaccharides",
  cpd90020 = "Oligo- & polysaccharides", cpd90021 = "Oligo- & polysaccharides",
  cpd90022 = "Oligo- & polysaccharides",
  # Sugar alcohols & alcohols (11): polyols + simple alcohols
  cpd00080 = "Sugar alcohols & alcohols", cpd00100 = "Sugar alcohols & alcohols",
  cpd00116 = "Sugar alcohols & alcohols", cpd00306 = "Sugar alcohols & alcohols",
  cpd00314 = "Sugar alcohols & alcohols", cpd00363 = "Sugar alcohols & alcohols",
  cpd00366 = "Sugar alcohols & alcohols", cpd00453 = "Sugar alcohols & alcohols",
  cpd00588 = "Sugar alcohols & alcohols", cpd01861 = "Sugar alcohols & alcohols",
  cpd03662 = "Sugar alcohols & alcohols",
  # Organic & sugar acids (20): TCA intermediates + SCFAs + sugar acids
  cpd00020 = "Organic & sugar acids", cpd00024 = "Organic & sugar acids",
  cpd00029 = "Organic & sugar acids", cpd00036 = "Organic & sugar acids",
  cpd00059 = "Organic & sugar acids", cpd00106 = "Organic & sugar acids",
  cpd00130 = "Organic & sugar acids", cpd00137 = "Organic & sugar acids",
  cpd00139 = "Organic & sugar acids", cpd00141 = "Organic & sugar acids",
  cpd00159 = "Organic & sugar acids", cpd00176 = "Organic & sugar acids",
  cpd00180 = "Organic & sugar acids", cpd00221 = "Organic & sugar acids",
  cpd00222 = "Organic & sugar acids", cpd00280 = "Organic & sugar acids",
  cpd00573 = "Organic & sugar acids", cpd00609 = "Organic & sugar acids",
  cpd00653 = "Organic & sugar acids", cpd22614 = "Organic & sugar acids",
  # Nucleosides & bases (11)
  cpd00092 = "Nucleosides & bases", cpd00182 = "Nucleosides & bases",
  cpd00184 = "Nucleosides & bases", cpd00246 = "Nucleosides & bases",
  cpd00249 = "Nucleosides & bases", cpd00307 = "Nucleosides & bases",
  cpd00309 = "Nucleosides & bases", cpd00311 = "Nucleosides & bases",
  cpd00355 = "Nucleosides & bases", cpd00367 = "Nucleosides & bases",
  cpd01217 = "Nucleosides & bases",
  # Other (19): cofactor artefacts, ROS, gases, aldehydes, sulfonates,
  # aromatic intermediates, bile acids, mucin glycans
  cpd00003 = "Other", cpd00006 = "Other", cpd00007 = "Other",
  cpd00010 = "Other", cpd00012 = "Other", cpd11606 = "Other",
  cpd00025 = "Other", cpd00055 = "Other",
  cpd00071 = "Other", cpd00229 = "Other", cpd00371 = "Other",
  cpd00448 = "Other", cpd00204 = "Other", cpd00210 = "Other",
  cpd00216 = "Other", cpd01318 = "Other", cpd03247 = "Other",
  cpd02992 = "Other", cpd11842 = "Other"
)

# ---------- Organism lists and metadata ----------
single_targets <- c("Ecoli_K12_MG1655", "Ecoli_PA3",
                    "Btheta_VPI5482",   "Btheta_KPPR3",
                    "Efaecalis_ATCC19433", "Efaecalis_68A",
                    "MGYG000292883", "MGYG000295553")

pan_targets <- c("MGYG000290784_pan", "MGYG000291361_pan",
                 "MGYG000291777_pan", "MGYG000292637_pan",
                 "MGYG000293427_pan", "MGYG000294127_pan",
                 "MGYG000295164_pan", "MGYG000295308_pan",
                 "MGYG000295316_pan")

all_meta <- tribble(
  ~organism,             ~display,             ~class,
  "Ecoli_K12_MG1655",    "E. coli K-12",        "Gammaproteobacteria",
  "Ecoli_PA3",           "E. coli PA-3",        "Gammaproteobacteria",
  "Btheta_VPI5482",      "B. theta VPI-5482",   "Bacteroidia",
  "Btheta_KPPR3",        "B. theta KPPR-3",     "Bacteroidia",
  "Efaecalis_ATCC19433", "E. faecalis ATCC",    "Bacilli",
  "Efaecalis_68A",       "E. faecalis 68A",     "Bacilli",
  "MGYG000292883",       "MGYG000292883",       "Lentisphaeria",
  "MGYG000295553",       "MGYG000295553",       "Kiritimatiellia",
  "MGYG000290784_pan",   "MGYG000290784",       "Bacteroidia",
  "MGYG000291361_pan",   "MGYG000291361",       "Clostridia",
  "MGYG000291777_pan",   "MGYG000291777",       "Clostridia",
  "MGYG000292637_pan",   "MGYG000292637",       "Bacteroidia",
  "MGYG000293427_pan",   "MGYG000293427",       "Clostridia",
  "MGYG000294127_pan",   "MGYG000294127",       "Clostridia",
  "MGYG000295164_pan",   "MGYG000295164",       "Clostridia",
  "MGYG000295308_pan",   "MGYG000295308",       "Clostridia",
  "MGYG000295316_pan",   "MGYG000295316",       "Clostridia"
)

y_axis_order <- c(
  "Ecoli_K12_MG1655", "Ecoli_PA3",
  "Efaecalis_ATCC19433", "Efaecalis_68A",
  "Btheta_VPI5482", "Btheta_KPPR3",
  "MGYG000290784_pan", "MGYG000292637_pan",
  "MGYG000291361_pan", "MGYG000291777_pan",
  "MGYG000293427_pan", "MGYG000294127_pan",
  "MGYG000295164_pan", "MGYG000295308_pan", "MGYG000295316_pan",
  "MGYG000292883", "MGYG000295553"
)

# ---------- Helpers ----------
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

shorten_compound <- function(x) {
  case_when(
    grepl("^N-Acetyl-beta-D-glucosaminyl-1,6", x) ~ "GlcNAc\u03b21-6",
    grepl("^N-Acetyl-beta-D-glucosaminyl-1,3", x) ~ "GlcNAc\u03b21-3",
    grepl("^starch \\(n=27", x)         ~ "Starch (n=27)",
    grepl("^starch \\(n=19", x)         ~ "Starch (n=19)",
    grepl("^Inulin", x)                  ~ "Inulin",
    x == "5-Methylthio-D-ribose"         ~ "5-MTR",
    x == "N-acetylneuraminate"           ~ "Neu5Ac",
    x == "N-Acetyl-D-glucosamine"        ~ "GlcNAc",
    x == "N-Acetyl-D-chondrosamine"      ~ "GalNAc",
    x == "Glycerol-3-phosphate"          ~ "Glycerol-3-P",
    x == "(R)-1,2-Propanediol"           ~ "(R)-1,2-PDO",
    x == "Nicotinamide ribonucleotide"   ~ "NMN",
    x == "2-keto-3-deoxygluconate"       ~ "KDG",
    x == "Glycochenodeoxycholate"        ~ "GCDCA",
    TRUE ~ x
  )
}

build_y_labels <- function(orgs) {
  m <- all_meta[match(orgs, all_meta$organism), ]
  setNames(sprintf("%s (%s)", m$display, m$class), orgs)
}
y_labels_named <- build_y_labels(y_axis_order)

# =============================================================
# Panel A - Auxotrophy heatmap
# =============================================================
pan_aux    <- read_excel(xlsx_path, sheet = "pandraft_all_auxotrophies")
single_aux <- read_excel(xlsx_path, sheet = "single_all_auxotrophies")

aux_all <- bind_rows(
  pan_aux    %>% filter(organism %in% pan_targets) %>%
    select(organism, compound),
  single_aux %>% filter(organism %in% single_targets) %>%
    select(organism, compound)
) %>%
  mutate(category = factor(aux_class(compound), levels = names(aux_pal)))

compound_order_A <- aux_all %>%
  count(compound, category, name = "n_orgs") %>%
  arrange(category, desc(n_orgs), compound) %>%
  pull(compound)

aux_long_A <- aux_all %>%
  mutate(organism = factor(organism, levels = rev(y_axis_order)),
         compound = factor(compound, levels = compound_order_A))

p_a <- ggplot(aux_long_A, aes(x = compound, y = organism, fill = category)) +
  geom_tile(colour = "grey85", linewidth = 0.4) +
  scale_fill_manual(values = aux_pal, name = "Category", drop = FALSE) +
  scale_y_discrete(labels = y_labels_named, expand = c(0, 0)) +
  scale_x_discrete(expand = c(0, 0)) +
  labs(title = "A. Auxotrophies across all metabolic models",
       x = NULL, y = NULL) +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid          = element_blank(),
    panel.background    = element_rect(fill = "grey97", colour = NA),
    plot.title          = element_text(face = "bold", size = 14, hjust = 0.5,
                                       margin = margin(b = 8)),
    plot.title.position = "plot",
    axis.text.x         = element_text(angle = 45, hjust = 1, size = 9, face = "bold"),
    axis.text.y         = element_text(face = "bold", colour = "grey15", size = 10),
    axis.ticks.y        = element_blank(),
    legend.title        = element_text(face = "bold", size = 12),
    legend.text         = element_text(face = "bold", size = 11),
    legend.key.size     = unit(0.7, "cm"),
    legend.position     = "right"
  )

# =============================================================
# Panel B - Substrate heatmap with category strip
# =============================================================
pan_sub    <- read_excel(xlsx_path, sheet = "pandraft_all_substrates")
single_sub <- read_excel(xlsx_path, sheet = "single_all_substrates")

# Apply active-substrate threshold (gapseq default 0.01 h^-1).
# Below-threshold rows in the supplementary sheets are exchange-flux
# artefacts (CO2, urea, methane, fatty acids, etc.) and are excluded.
sub_all <- bind_rows(
  pan_sub    %>% filter(organism %in% pan_targets,
                        Net_Growth >= SUB_THRESHOLD) %>%
    select(organism, cpd_id, compound, Net_Growth),
  single_sub %>% filter(organism %in% single_targets,
                        Net_Growth >= SUB_THRESHOLD) %>%
    select(organism, cpd_id, compound, Net_Growth)
) %>%
  mutate(category = factor(compound_category[cpd_id], levels = cat_levels))

cat(sprintf("Active substrates after threshold (>= %.2f h^-1): %d rows, %d compounds\n",
            SUB_THRESHOLD, nrow(sub_all), n_distinct(sub_all$cpd_id)))

# Sanity check: any compound not assigned a category?
unmapped <- sub_all %>% filter(is.na(category)) %>% distinct(cpd_id, compound)
if (nrow(unmapped) > 0) {
  warning(sprintf("Unmapped compounds (no category): %s",
                  paste(unmapped$cpd_id, collapse = ", ")))
  cat("\n=== UNMAPPED compounds (need to be added to compound_category) ===\n")
  print(unmapped)
} else {
  cat("All active substrates are categorised.\n")
}

# Within-category ordering by organism count
compound_hits <- sub_all %>%
  group_by(cpd_id, compound, category) %>%
  summarise(n_hits = n_distinct(organism), .groups = "drop") %>%
  arrange(category, desc(n_hits), compound)

compound_order_real <- compound_hits$compound

# Insert one blank placeholder column between adjacent categories
compound_with_gaps <- character()
prev_cat <- NA_character_
gap_id   <- 0
for (i in seq_len(nrow(compound_hits))) {
  cur_cat <- as.character(compound_hits$category[i])
  if (!is.na(prev_cat) && cur_cat != prev_cat) {
    gap_id <- gap_id + 1
    compound_with_gaps <- c(compound_with_gaps, paste0("__gap_", gap_id))
  }
  compound_with_gaps <- c(compound_with_gaps, compound_hits$compound[i])
  prev_cat <- cur_cat
}

x_label_fn <- function(x) ifelse(grepl("^__gap_", x), "", shorten_compound(x))

# ---------- Panel B strip ----------
strip_df <- compound_hits %>%
  mutate(compound = factor(compound, levels = compound_with_gaps),
         y_pos    = "")

p_strip_B <- ggplot(strip_df, aes(x = compound, y = y_pos, fill = category)) +
  geom_tile(colour = "white", linewidth = 0.3) +
  scale_fill_manual(values = cat_pal, name = "Category", drop = FALSE) +
  scale_x_discrete(drop = FALSE, expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  labs(title = "B. Predicted substrate preference across all metabolic models",
       x = NULL, y = NULL) +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid          = element_blank(),
    panel.background    = element_blank(),
    plot.title          = element_text(face = "bold", size = 14, hjust = 0.5,
                                       margin = margin(b = 6)),
    plot.title.position = "plot",
    axis.text           = element_blank(),
    axis.ticks          = element_blank(),
    legend.title        = element_text(face = "bold", size = 12),
    legend.text         = element_text(face = "bold", size = 10),
    legend.key.size     = unit(0.55, "cm"),
    legend.position     = "right",
    plot.margin         = margin(0, 0, 0, 0)
  )

# ---------- Panel B main heatmap ----------
heat_df_B <- expand.grid(organism = y_axis_order,
                         compound = compound_order_real,
                         stringsAsFactors = FALSE) %>%
  left_join(sub_all %>% select(organism, compound, Net_Growth),
            by = c("organism", "compound")) %>%
  mutate(Net_Growth = ifelse(is.na(Net_Growth), 0, Net_Growth),
         organism   = factor(organism, levels = rev(y_axis_order)),
         compound   = factor(compound, levels = compound_with_gaps))

max_growth  <- max(heat_df_B$Net_Growth)
fill_breaks <- c(0, 0.1, 0.3, 1, 3, 6)
fill_breaks <- fill_breaks[fill_breaks <= max_growth * 1.05]

p_heat_B <- ggplot(heat_df_B, aes(x = compound, y = organism, fill = Net_Growth)) +
  geom_tile(colour = "grey92", linewidth = 0.2) +
  scale_fill_gradient(
    low    = "grey97",
    high   = "#1A5490",
    name   = expression(bold(paste("Net growth (h"^-1, ")"))),
    trans  = pseudo_log_trans(sigma = LOG_SIGMA, base = 10),
    breaks = fill_breaks,
    labels = function(x) ifelse(x == 0, "0", sprintf("%g", x)),
    limits = c(0, max_growth)
  ) +
  scale_x_discrete(drop = FALSE, expand = c(0, 0), labels = x_label_fn) +
  scale_y_discrete(labels = y_labels_named, expand = c(0, 0)) +
  labs(x = NULL, y = NULL) +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid          = element_blank(),
    panel.background    = element_rect(fill = "grey99", colour = NA),
    axis.text.x         = element_text(angle = 90, hjust = 1, vjust = 0.5,
                                       size = 6.5, face = "bold"),
    axis.text.y         = element_text(face = "bold", colour = "grey15",
                                       size = 10),
    axis.ticks          = element_blank(),
    legend.title        = element_text(face = "bold", size = 12),
    legend.text         = element_text(face = "bold", size = 11),
    legend.key.height   = unit(1.4, "cm"),
    legend.key.width    = unit(0.5, "cm"),
    legend.position     = "right",
    plot.margin         = margin(0, 0, 0, 0)
  )

panel_B <- p_strip_B / p_heat_B + plot_layout(heights = c(0.05, 1))

# =============================================================
# Final figure
# =============================================================
fig <- p_a / panel_B + plot_layout(heights = c(0.7, 1.4))

out_path <- file.path(out_dir, "Fig3_AB_heatmaps.png")
ggsave(out_path, fig,
       width = 16, height = 14, units = "in",
       dpi = 1000, bg = "white")

cat(sprintf("\nSaved: %s\n", out_path))
cat(sprintf("Panel B: %d compounds, %d categories\n",
            length(compound_order_real), length(unique(compound_hits$category))))
cat(sprintf("Pseudo-log: sigma = %g, base = 10\n", LOG_SIGMA))
cat(sprintf("Max Net_Growth = %.3f\n", max_growth))

# ---------- Diagnostics ----------
cat("\n=== Compound counts per category (post-threshold) ===\n")
print(compound_hits %>% count(category, name = "n_compounds"))

cat("\n=== Y-axis order ===\n")
m_summary <- all_meta[match(y_axis_order, all_meta$organism), ]
for (i in seq_along(y_axis_order)) {
  cat(sprintf("  %2d. %-30s  %s\n",
              i, y_axis_order[i], m_summary$class[i]))
}
