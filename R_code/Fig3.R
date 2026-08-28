# =============================================================
# Figure 3 - panels A and B (heatmaps)
#
# Sample set: 42 organisms
#   6 reference single models (3 paired clinical-rumen)
#   34 SGB pan-Draft models
#   2 ARG-MAG single models (placed at the bottom for clarity)
#
# Compound classification: 126 active compounds across 11 manual
# categories (KEGG BRITE br08001 + ModelSEED Biochemistry sourcing,
# see Supplementary Table S1). Two adaptations relative to KEGG:
#   - Carbohydrates split by degree of polymerisation
#   - Aldehydes, Inorganic gases and Host-derived substrates pulled
#     out of the previous "Other" pool.
# Cofactor-like exchange artefacts (NAD, NADP, CoA, PPi, H2O2,
# chorismate, three menaquinone variants) remain as "Other".
# Long-chain fatty acids (C6-C18) were not predicted as active sole
# carbon sources in any of the 42 GEMs (consistent with anaerobic
# rumen biology), so the Fatty acids category is omitted.
#
# Y-axis order (top to bottom), grouped by Class with paired
# references adjacent to same-Class pan models, then ARG-MAGs:
#   Bacteroidia (2 ref + 11 pan)
#   Clostridia  (16 pan)
#   Bacilli     (2 ref +  3 pan)
#   Negativicutes (1 pan)
#   Coriobacteriia (1 pan)
#   Gammaproteobacteria (2 ref)
#   Methanobacteria (2 archaea pan)
#   ARG-MAG (Lentisphaeria, Kiritimatiellia)
#
# Panel A: Auxotrophy heatmap, cells coloured by category.
# Panel B: Substrate preference heatmap.
#   - substrate threshold Net_Growth >= 0.01 h^-1 (gapseq default)
#   - 6 fba-only SGBs have Total_Flux NA but Net_Growth comparable
#   - within each category, compounds ordered by hit count
#   - blank columns separate categories
#   - pseudo_log fill (sigma = 0.1) preserves small values
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
# Canonical Supplementary Tables xlsx (top-level). Sheets have 3 header
# rows before the column header (skip = 3). S3a / S3b combine Pan-Draft
# and Single models; split by model_type below and re-append "_pan" suffix
# on Pan-Draft ids so pan_targets / y_axis_order / all_meta continue to
# match (they use the "_pan" convention throughout).
xlsx_path <- "C:/Users/CFL/OneDrive - Queen's University Belfast/PhD_onedrive/Manuscript4_gapseq/M4_Supplementary_Tables.xlsx"
out_dir   <- "C:/Users/CFL/OneDrive - Queen's University Belfast/PhD_onedrive/Manuscript4_gapseq/figures"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

LOG_SIGMA     <- 0.1
SUB_THRESHOLD <- 0.01

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
  "Aldehydes"                 = "#D6B05E",
  "Nucleosides & bases"       = "#4F8F9C",
  "Host-derived substrates"   = "#A05A66",
  "Inorganic gases"           = "#6BA292",
  "Other"                     = "#B5BABC"
)
cat_levels <- names(cat_pal)

# ---------- Compound -> category dictionary (131 entries) ----------
compound_category <- c(
  # Amino acids (19)
  cpd00023="Amino acids", cpd00033="Amino acids", cpd00035="Amino acids",
  cpd00039="Amino acids", cpd00041="Amino acids", cpd00051="Amino acids",
  cpd00053="Amino acids", cpd00054="Amino acids", cpd00060="Amino acids",
  cpd00064="Amino acids", cpd00065="Amino acids", cpd00084="Amino acids",
  cpd00117="Amino acids", cpd00129="Amino acids", cpd00132="Amino acids",
  cpd00156="Amino acids", cpd00161="Amino acids", cpd00281="Amino acids",
  cpd00550="Amino acids",
  # Mono- & disaccharides (17)
  cpd00027="Mono- & disaccharides", cpd00076="Mono- & disaccharides",
  cpd00082="Mono- & disaccharides", cpd00105="Mono- & disaccharides",
  cpd00108="Mono- & disaccharides", cpd00138="Mono- & disaccharides",
  cpd00154="Mono- & disaccharides", cpd00158="Mono- & disaccharides",
  cpd00179="Mono- & disaccharides", cpd00185="Mono- & disaccharides",
  cpd00208="Mono- & disaccharides", cpd00224="Mono- & disaccharides",
  cpd00396="Mono- & disaccharides", cpd00751="Mono- & disaccharides",
  cpd00794="Mono- & disaccharides", cpd01354="Mono- & disaccharides",
  cpd03198="Mono- & disaccharides",
  # Amino sugars (5; +Lacto-N-biose)
  cpd00122="Amino sugars", cpd00232="Amino sugars",
  cpd00276="Amino sugars", cpd00832="Amino sugars",
  cpd03808="Amino sugars",                   # Lacto-N-biose
  # Oligo- & polysaccharides (17)
  cpd00155="Oligo- & polysaccharides", cpd00382="Oligo- & polysaccharides",
  cpd01133="Oligo- & polysaccharides", cpd01262="Oligo- & polysaccharides",
  cpd01399="Oligo- & polysaccharides", cpd11976="Oligo- & polysaccharides",
  cpd15494="Oligo- & polysaccharides", cpd28763="Oligo- & polysaccharides",
  cpd90003="Oligo- & polysaccharides", cpd90004="Oligo- & polysaccharides",
  cpd90005="Oligo- & polysaccharides", cpd90006="Oligo- & polysaccharides",
  cpd90007="Oligo- & polysaccharides", cpd90008="Oligo- & polysaccharides",
  cpd90020="Oligo- & polysaccharides", cpd90021="Oligo- & polysaccharides",
  cpd90022="Oligo- & polysaccharides",
  # Sugar alcohols & alcohols (11)
  cpd00080="Sugar alcohols & alcohols", cpd00100="Sugar alcohols & alcohols",
  cpd00116="Sugar alcohols & alcohols", cpd00306="Sugar alcohols & alcohols",
  cpd00314="Sugar alcohols & alcohols", cpd00363="Sugar alcohols & alcohols",
  cpd00366="Sugar alcohols & alcohols", cpd00453="Sugar alcohols & alcohols",
  cpd00588="Sugar alcohols & alcohols", cpd01861="Sugar alcohols & alcohols",
  cpd03662="Sugar alcohols & alcohols",
  # Organic & sugar acids (21; +Formate)
  cpd00020="Organic & sugar acids", cpd00024="Organic & sugar acids",
  cpd00029="Organic & sugar acids", cpd00036="Organic & sugar acids",
  cpd00047="Organic & sugar acids",          # Formate
  cpd00059="Organic & sugar acids", cpd00106="Organic & sugar acids",
  cpd00130="Organic & sugar acids", cpd00137="Organic & sugar acids",
  cpd00139="Organic & sugar acids", cpd00141="Organic & sugar acids",
  cpd00159="Organic & sugar acids", cpd00176="Organic & sugar acids",
  cpd00180="Organic & sugar acids", cpd00221="Organic & sugar acids",
  cpd00222="Organic & sugar acids", cpd00280="Organic & sugar acids",
  cpd00573="Organic & sugar acids", cpd00609="Organic & sugar acids",
  cpd00653="Organic & sugar acids", cpd22614="Organic & sugar acids",
  # Aldehydes (5)
  cpd00055="Aldehydes", cpd00071="Aldehydes",
  cpd00229="Aldehydes", cpd00371="Aldehydes",
  cpd00448="Aldehydes",
  # Nucleosides & bases (11)
  cpd00092="Nucleosides & bases", cpd00182="Nucleosides & bases",
  cpd00184="Nucleosides & bases", cpd00246="Nucleosides & bases",
  cpd00249="Nucleosides & bases", cpd00307="Nucleosides & bases",
  cpd00309="Nucleosides & bases", cpd00311="Nucleosides & bases",
  cpd00355="Nucleosides & bases", cpd00367="Nucleosides & bases",
  cpd01217="Nucleosides & bases",
  # Host-derived substrates (5)
  cpd00210="Host-derived substrates", cpd01318="Host-derived substrates",
  cpd03247="Host-derived substrates", cpd02992="Host-derived substrates",
  cpd11842="Host-derived substrates",
  # Inorganic gases (6)
  cpd00007="Inorganic gases",                # O2
  cpd00204="Inorganic gases",                # CO
  cpd00011="Inorganic gases",                # CO2
  cpd11640="Inorganic gases",                # H2
  cpd00239="Inorganic gases",                # H2S
  cpd00324="Inorganic gases",                # Methanethiol
  # Other (9): cofactor / inorganic exchange artefacts
  cpd00003="Other",                          # NAD
  cpd00006="Other",                          # NADP
  cpd00010="Other",                          # CoA
  cpd00012="Other",                          # PPi
  cpd00025="Other",                          # H2O2 (ROS, not a gas)
  cpd00216="Other",                          # Chorismate
  cpd11606="Other",                          # Menaquinone-7
  cpd11451="Other",                          # Menaquinol-7
  cpd17026="Other"                           # Methylmenaquinol-7
)

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
    x == "Glycocholate"                  ~ "GCA",
    x == "Lacto-N-biose"                 ~ "LNB",
    x == "Carbon monoxide"               ~ "CO",
    x == "Methanethiol"                  ~ "CH3SH",
    x == "Hydrogen"                      ~ "H2",
    x == "Menaquinone 7"                 ~ "MK7",
    x == "Menaquinol 7"                  ~ "MK7H2",
    x == "methymenaquinol 7"             ~ "MMK7",
    x == "Pyrophosphate"                 ~ "PPi",
    x == "Hexanoate"                     ~ "C6",
    x == "octanoate"                     ~ "C8",
    x == "Decanoate"                     ~ "C10",
    x == "Dodecanoic acid"               ~ "C12",
    x == "Octadecanoate"                 ~ "C18",
    TRUE ~ x
  )
}

# ---------- Organism lists ----------
single_targets <- c("Ecoli_K12_MG1655", "Ecoli_PA3",
                    "Btheta_VPI5482",   "Btheta_KPPR3",
                    "Efaecalis_ATCC19433", "Efaecalis_68A",
                    "MGYG000292883", "MGYG000295553")

# 34 SGB pan models
pan_targets <- c(
  "MGYG000290784_pan","MGYG000290832_pan","MGYG000290859_pan","MGYG000291338_pan",
  "MGYG000291361_pan","MGYG000291422_pan","MGYG000291777_pan","MGYG000291885_pan",
  "MGYG000291990_pan","MGYG000291991_pan","MGYG000292023_pan","MGYG000292637_pan",
  "MGYG000292896_pan","MGYG000293201_pan","MGYG000293286_pan","MGYG000293427_pan",
  "MGYG000293532_pan","MGYG000293741_pan","MGYG000293943_pan","MGYG000294070_pan",
  "MGYG000294127_pan","MGYG000294153_pan","MGYG000294175_pan","MGYG000294280_pan",
  "MGYG000294398_pan","MGYG000294423_pan","MGYG000294647_pan","MGYG000294799_pan",
  "MGYG000295164_pan","MGYG000295175_pan","MGYG000295242_pan","MGYG000295308_pan",
  "MGYG000295316_pan","MGYG000295318_pan"
)

# ---------- Metadata: 42 organisms with class + family ----------
all_meta <- tribble(
  ~organism,             ~display,            ~class,                 ~role,
  # 6 references
  "Btheta_VPI5482",      "B. theta VPI-5482", "Bacteroidia",          "ref",
  "Btheta_KPPR3",        "B. theta KPPR-3",   "Bacteroidia",          "ref",
  "Efaecalis_ATCC19433", "E. faecalis ATCC",  "Bacilli",              "ref",
  "Efaecalis_68A",       "E. faecalis 68A",   "Bacilli",              "ref",
  "Ecoli_K12_MG1655",    "E. coli K-12",      "Gammaproteobacteria",  "ref",
  "Ecoli_PA3",           "E. coli PA-3",      "Gammaproteobacteria",  "ref",
  # 11 Bacteroidia pan
  "MGYG000290784_pan",   "MGYG000290784",     "Bacteroidia",          "pan",
  "MGYG000292023_pan",   "MGYG000292023",     "Bacteroidia",          "pan",
  "MGYG000293741_pan",   "MGYG000293741",     "Bacteroidia",          "pan",
  "MGYG000294153_pan",   "MGYG000294153",     "Bacteroidia",          "pan",
  "MGYG000295242_pan",   "MGYG000295242",     "Bacteroidia",          "pan",
  "MGYG000294070_pan",   "MGYG000294070",     "Bacteroidia",          "pan",
  "MGYG000291991_pan",   "MGYG000291991",     "Bacteroidia",          "pan",
  "MGYG000293943_pan",   "MGYG000293943",     "Bacteroidia",          "pan",
  "MGYG000292637_pan",   "MGYG000292637",     "Bacteroidia",          "pan",
  "MGYG000293201_pan",   "MGYG000293201",     "Bacteroidia",          "pan",
  "MGYG000294423_pan",   "MGYG000294423",     "Bacteroidia",          "pan",
  # 16 Clostridia pan (Acutalibacteraceae 7, Lachnospiraceae 6, Oscillospiraceae 3)
  "MGYG000290832_pan",   "MGYG000290832",     "Clostridia",           "pan",
  "MGYG000290859_pan",   "MGYG000290859",     "Clostridia",           "pan",
  "MGYG000291777_pan",   "MGYG000291777",     "Clostridia",           "pan",
  "MGYG000294175_pan",   "MGYG000294175",     "Clostridia",           "pan",
  "MGYG000294398_pan",   "MGYG000294398",     "Clostridia",           "pan",
  "MGYG000295164_pan",   "MGYG000295164",     "Clostridia",           "pan",
  "MGYG000295316_pan",   "MGYG000295316",     "Clostridia",           "pan",
  "MGYG000291422_pan",   "MGYG000291422",     "Clostridia",           "pan",
  "MGYG000293286_pan",   "MGYG000293286",     "Clostridia",           "pan",
  "MGYG000293427_pan",   "MGYG000293427",     "Clostridia",           "pan",
  "MGYG000294127_pan",   "MGYG000294127",     "Clostridia",           "pan",
  "MGYG000294799_pan",   "MGYG000294799",     "Clostridia",           "pan",
  "MGYG000295308_pan",   "MGYG000295308",     "Clostridia",           "pan",
  "MGYG000291338_pan",   "MGYG000291338",     "Clostridia",           "pan",
  "MGYG000291361_pan",   "MGYG000291361",     "Clostridia",           "pan",
  "MGYG000295318_pan",   "MGYG000295318",     "Clostridia",           "pan",
  # 3 Bacilli pan
  "MGYG000291885_pan",   "MGYG000291885",     "Bacilli",              "pan",
  "MGYG000294280_pan",   "MGYG000294280",     "Bacilli",              "pan",
  "MGYG000295175_pan",   "MGYG000295175",     "Bacilli",              "pan",
  # 1 Negativicutes pan
  "MGYG000294647_pan",   "MGYG000294647",     "Negativicutes",        "pan",
  # 1 Coriobacteriia pan
  "MGYG000292896_pan",   "MGYG000292896",     "Coriobacteriia",       "pan",
  # 2 Methanobacteria archaea pan
  "MGYG000291990_pan",   "MGYG000291990",     "Methanobacteria",      "pan",
  "MGYG000293532_pan",   "MGYG000293532",     "Methanobacteria",      "pan",
  # 2 ARG-MAGs (bottom, separated)
  "MGYG000292883",       "MGYG000292883",     "Lentisphaeria",        "argmag",
  "MGYG000295553",       "MGYG000295553",     "Kiritimatiellia",      "argmag"
)

# Y-axis order = the order of rows in all_meta above (top -> bottom)
y_axis_order <- all_meta$organism

build_y_labels <- function(orgs) {
  m <- all_meta[match(orgs, all_meta$organism), ]
  setNames(sprintf("%s (%s)", m$display, m$class), orgs)
}
y_labels_named <- build_y_labels(y_axis_order)

# =============================================================
# Panel A - Auxotrophy heatmap
# =============================================================
# S3a combines Pan-Draft + Single auxotrophies. Read once, split by
# model_type; re-append "_pan" to Pan-Draft ids so pan_targets matches.
aux_raw <- read_excel(xlsx_path, sheet = "S3a_Essential_factors", skip = 3)
pan_aux <- aux_raw %>%
  filter(model_type == "Pan-Draft") %>%
  mutate(organism = paste0(organism, "_pan")) %>%
  select(organism, compound)
single_aux <- aux_raw %>%
  filter(model_type == "Single-genome") %>%
  select(organism, compound)

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
  labs(title = "Predicted essential-factor profiles across all metabolic models",
       x = NULL, y = NULL) +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid          = element_blank(),
    panel.background    = element_rect(fill = "grey97", colour = NA),
    plot.title          = element_text(face = "bold", size = 20, hjust = 0.5,
                                       margin = margin(b = 8)),
    plot.title.position = "plot",
    axis.text.x         = element_text(angle = 45, hjust = 1, size = 12, face = "bold"),
    axis.text.y         = element_text(face = "bold", colour = "grey15", size = 12),
    axis.ticks.y        = element_blank(),
    legend.title        = element_text(face = "bold", size = 15),
    legend.text         = element_text(face = "bold", size = 13),
    legend.key.size     = unit(0.9, "cm"),
    legend.position     = "right"
  )

# =============================================================
# Panel B - Substrate heatmap with category strip
# =============================================================
# S3b combines Pan-Draft + Single active substrates. Read once, split by
# model_type. Column `net_growth` is renamed to `Net_Growth` to keep the
# rest of this script (case-sensitive) unchanged.
sub_raw <- read_excel(xlsx_path, sheet = "S3b_Substrates_all", skip = 3) %>%
  rename(Net_Growth  = net_growth,
         Total_Flux  = total_flux,
         Growth_Rate = growth_rate,
         Reaction_ID = reaction_id)
pan_sub <- sub_raw %>%
  filter(model_type == "Pan-Draft") %>%
  mutate(organism = paste0(organism, "_pan"))
single_sub <- sub_raw %>%
  filter(model_type == "Single-genome")

sub_all <- bind_rows(
  pan_sub    %>% filter(organism %in% pan_targets,
                        Net_Growth >= SUB_THRESHOLD) %>%
    select(organism, cpd_id, compound, Net_Growth),
  single_sub %>% filter(organism %in% single_targets,
                        Net_Growth >= SUB_THRESHOLD) %>%
    select(organism, cpd_id, compound, Net_Growth)
) %>%
  mutate(category = factor(compound_category[cpd_id], levels = cat_levels))

cat(sprintf("Active substrates (>= %.2f h^-1): %d rows, %d compounds\n",
            SUB_THRESHOLD, nrow(sub_all), n_distinct(sub_all$cpd_id)))

unmapped <- sub_all %>% filter(is.na(category)) %>% distinct(cpd_id, compound)
if (nrow(unmapped) > 0) {
  warning(sprintf("Unmapped compounds: %s", paste(unmapped$cpd_id, collapse = ", ")))
  print(unmapped)
} else {
  cat("All active substrates are categorised.\n")
}

compound_hits <- sub_all %>%
  group_by(cpd_id, compound, category) %>%
  summarise(n_hits = n_distinct(organism), .groups = "drop") %>%
  arrange(category, desc(n_hits), compound)

compound_order_real <- compound_hits$compound

# Insert blank columns between categories
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
  labs(title = "Predicted substrate preference across all metabolic models",
       x = NULL, y = NULL) +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid          = element_blank(),
    panel.background    = element_blank(),
    plot.title          = element_text(face = "bold", size = 20, hjust = 0.5,
                                       margin = margin(b = 6)),
    plot.title.position = "plot",
    axis.text           = element_blank(),
    axis.ticks          = element_blank(),
    legend.title        = element_text(face = "bold", size = 15),
    legend.text         = element_text(face = "bold", size = 13),
    legend.key.size     = unit(0.9, "cm"),
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
  geom_tile(colour = "grey92", linewidth = 0.15) +
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
                                       size = 9, face = "bold"),
    axis.text.y         = element_text(face = "bold", colour = "grey15",
                                       size = 12),
    axis.ticks          = element_blank(),
    legend.title        = element_text(face = "bold", size = 15),
    legend.text         = element_text(face = "bold", size = 13),
    legend.key.height   = unit(1.6, "cm"),
    legend.key.width    = unit(0.6, "cm"),
    legend.position     = "right",
    plot.margin         = margin(0, 0, 0, 0)
  )

panel_B <- p_strip_B / p_heat_B +
  plot_layout(heights = c(0.04, 1), guides = "collect") &
  theme(legend.position      = "right",
        legend.justification = "center",
        legend.box           = "vertical",
        legend.box.just      = "left")

# =============================================================
# Final figure
# =============================================================

fig <- p_a / panel_B +
  plot_layout(heights = c(0.9, 1.0)) +
  plot_annotation(tag_levels = list(c("A", "B"))) &
  theme(plot.tag          = element_text(face = "bold", size = 40),
        plot.tag.position = c(0.005, 0.985))
out_path <- file.path(out_dir, "Fig3.png")
ggsave(out_path, fig,
       width = 22, height = 24, units = "in",
       dpi = 1000, bg = "white")

cat(sprintf("\nSaved: %s\n", out_path))
cat(sprintf("Y-axis: %d organisms\n", length(y_axis_order)))
cat(sprintf("Panel B: %d compounds across %d categories\n",
            length(compound_order_real), n_distinct(compound_hits$category)))

# ---------- Diagnostics ----------
cat("\n=== Compound counts per category (post-threshold) ===\n")
print(compound_hits %>% count(category, name = "n_compounds"))

cat("\n=== Y-axis order ===\n")
m_summary <- all_meta[match(y_axis_order, all_meta$organism), ]
for (i in seq_along(y_axis_order)) {
  role_tag <- c(ref="REF", pan="PAN", argmag="ARG-MAG")[m_summary$role[i]]
  cat(sprintf("  %2d. %-8s %-25s  %s\n",
              i, role_tag, m_summary$display[i], m_summary$class[i]))
}
