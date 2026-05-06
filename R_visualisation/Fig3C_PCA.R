# =============================================================
# Figure 3 - Panel C: weighted 3D PCA on metabolic phenotype
#
# Inputs (17 x 10 matrix):
#   3 auxotrophy categories (Amino acid, Vitamin/cofactor, Metal)
#   7 substrate categories (Amino acids, Mono- & disaccharides,
#     Amino sugars, Oligo- & polysaccharides, Sugar alcohols &
#     alcohols, Organic & sugar acids, Nucleosides & bases)
#   "Other" category dropped from both blocks.
#
# Sample weights: sqrt(orf_coverage), sum-normalised to n=17.
#   Pan models -> orf_cov_max from pandraft_all_stage1_summary
#   Single models -> orf_coverage_pct from single_all_reconstruction
#
# Feature weights: block-balancing
#   aux 3 columns weighted 1/sqrt(3)
#   sub 7 columns weighted 1/sqrt(7)
#
# PCA: FactoMineR::PCA with scale.unit = TRUE
#
# Outputs:
#   Fig3C_PCA_3D.png         - main 3D scatter (plot3D::scatter3D)
#   Fig3C_PCA_diagnostic.png - scree + loadings biplots (ggplot)
#   Fig3C_PCA_scores.csv     - PC1-3 scores per organism
# =============================================================

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(FactoMineR)
  library(plot3D)
  library(ggplot2)
  library(patchwork)
})

# ---------- Configuration ----------
xlsx_path <- "C:/Users/CFL/OneDrive - Queen's University Belfast/PhD_onedrive/Manuscript4/M4_Supplement.xlsx"
out_dir   <- "C:/Users/CFL/OneDrive - Queen's University Belfast/PhD_onedrive/Manuscript4/figures"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# Palettes
group_pal <- c(
  "PCG Reference"   = "#D4A017",
  "Rumen Reference" = "#2E7D32",
  "Pan model"       = "#2E86AB",
  "ARG-MAG"         = "#B22222"
)

phylum_pch <- c(
  "Pseudomonadota"    = 16,  # filled circle
  "Bacillota"         = 17,  # filled triangle
  "Bacteroidota"      = 15,  # filled square
  "Verrucomicrobiota" = 18   # filled diamond
)

# ---------- Compound -> category dictionary ----------
compound_category <- c(
  # Amino acids
  cpd00023 = "Amino acids", cpd00033 = "Amino acids",
  cpd00035 = "Amino acids", cpd00041 = "Amino acids",
  cpd00053 = "Amino acids", cpd00054 = "Amino acids",
  cpd00060 = "Amino acids", cpd00064 = "Amino acids",
  cpd00065 = "Amino acids", cpd00084 = "Amino acids",
  cpd00117 = "Amino acids", cpd00132 = "Amino acids",
  cpd00161 = "Amino acids", cpd00550 = "Amino acids",
  # Mono- & disaccharides
  cpd00027 = "Mono- & disaccharides", cpd00076 = "Mono- & disaccharides",
  cpd00082 = "Mono- & disaccharides", cpd00105 = "Mono- & disaccharides",
  cpd00108 = "Mono- & disaccharides", cpd00138 = "Mono- & disaccharides",
  cpd00154 = "Mono- & disaccharides", cpd00158 = "Mono- & disaccharides",
  cpd00179 = "Mono- & disaccharides", cpd00185 = "Mono- & disaccharides",
  cpd00208 = "Mono- & disaccharides", cpd00224 = "Mono- & disaccharides",
  cpd00396 = "Mono- & disaccharides", cpd00751 = "Mono- & disaccharides",
  cpd00794 = "Mono- & disaccharides", cpd01354 = "Mono- & disaccharides",
  cpd03198 = "Mono- & disaccharides",
  # Amino sugars
  cpd00122 = "Amino sugars", cpd00232 = "Amino sugars",
  cpd00276 = "Amino sugars", cpd00832 = "Amino sugars",
  # Oligo- & polysaccharides
  cpd00155 = "Oligo- & polysaccharides", cpd00382 = "Oligo- & polysaccharides",
  cpd01133 = "Oligo- & polysaccharides", cpd01262 = "Oligo- & polysaccharides",
  cpd01399 = "Oligo- & polysaccharides", cpd11976 = "Oligo- & polysaccharides",
  cpd15494 = "Oligo- & polysaccharides", cpd28763 = "Oligo- & polysaccharides",
  cpd90003 = "Oligo- & polysaccharides", cpd90004 = "Oligo- & polysaccharides",
  cpd90005 = "Oligo- & polysaccharides", cpd90006 = "Oligo- & polysaccharides",
  cpd90007 = "Oligo- & polysaccharides", cpd90008 = "Oligo- & polysaccharides",
  cpd90020 = "Oligo- & polysaccharides", cpd90021 = "Oligo- & polysaccharides",
  cpd90022 = "Oligo- & polysaccharides",
  # Sugar alcohols & alcohols
  cpd00080 = "Sugar alcohols & alcohols", cpd00100 = "Sugar alcohols & alcohols",
  cpd00306 = "Sugar alcohols & alcohols", cpd00314 = "Sugar alcohols & alcohols",
  cpd00363 = "Sugar alcohols & alcohols", cpd00366 = "Sugar alcohols & alcohols",
  cpd00453 = "Sugar alcohols & alcohols", cpd00588 = "Sugar alcohols & alcohols",
  cpd01861 = "Sugar alcohols & alcohols", cpd03662 = "Sugar alcohols & alcohols",
  # Organic & sugar acids
  cpd00020 = "Organic & sugar acids", cpd00024 = "Organic & sugar acids",
  cpd00059 = "Organic & sugar acids", cpd00106 = "Organic & sugar acids",
  cpd00130 = "Organic & sugar acids", cpd00137 = "Organic & sugar acids",
  cpd00139 = "Organic & sugar acids", cpd00159 = "Organic & sugar acids",
  cpd00176 = "Organic & sugar acids", cpd00221 = "Organic & sugar acids",
  cpd00222 = "Organic & sugar acids", cpd00280 = "Organic & sugar acids",
  cpd00573 = "Organic & sugar acids", cpd00609 = "Organic & sugar acids",
  cpd00653 = "Organic & sugar acids", cpd22614 = "Organic & sugar acids",
  # Nucleosides & bases
  cpd00092 = "Nucleosides & bases", cpd00182 = "Nucleosides & bases",
  cpd00184 = "Nucleosides & bases", cpd00246 = "Nucleosides & bases",
  cpd00249 = "Nucleosides & bases", cpd00307 = "Nucleosides & bases",
  cpd00311 = "Nucleosides & bases", cpd00355 = "Nucleosides & bases",
  cpd00367 = "Nucleosides & bases", cpd01217 = "Nucleosides & bases",
  # Other (will be filtered out)
  cpd00003 = "Other", cpd00006 = "Other", cpd00007 = "Other",
  cpd00012 = "Other", cpd11606 = "Other",
  cpd00071 = "Other", cpd00371 = "Other", cpd00448 = "Other",
  cpd01318 = "Other", cpd03247 = "Other",
  cpd02992 = "Other", cpd11842 = "Other"
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

# ---------- Organism metadata ----------
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
  ~organism,             ~display,             ~group,            ~phylum,
  "Ecoli_K12_MG1655",    "E. coli K-12",        "PCG Reference",   "Pseudomonadota",
  "Ecoli_PA3",           "E. coli PA-3",        "Rumen Reference", "Pseudomonadota",
  "Btheta_VPI5482",      "B. theta VPI-5482",   "PCG Reference",   "Bacteroidota",
  "Btheta_KPPR3",        "B. theta KPPR-3",     "Rumen Reference", "Bacteroidota",
  "Efaecalis_ATCC19433", "E. faecalis ATCC",    "PCG Reference",   "Bacillota",
  "Efaecalis_68A",       "E. faecalis 68A",     "Rumen Reference", "Bacillota",
  "MGYG000292883",       "MGYG000292883",       "ARG-MAG",         "Verrucomicrobiota",
  "MGYG000295553",       "MGYG000295553",       "ARG-MAG",         "Verrucomicrobiota",
  "MGYG000290784_pan",   "MGYG000290784",       "Pan model",       "Bacteroidota",
  "MGYG000291361_pan",   "MGYG000291361",       "Pan model",       "Bacillota",
  "MGYG000291777_pan",   "MGYG000291777",       "Pan model",       "Bacillota",
  "MGYG000292637_pan",   "MGYG000292637",       "Pan model",       "Bacteroidota",
  "MGYG000293427_pan",   "MGYG000293427",       "Pan model",       "Bacillota",
  "MGYG000294127_pan",   "MGYG000294127",       "Pan model",       "Bacillota",
  "MGYG000295164_pan",   "MGYG000295164",       "Pan model",       "Bacillota",
  "MGYG000295308_pan",   "MGYG000295308",       "Pan model",       "Bacillota",
  "MGYG000295316_pan",   "MGYG000295316",       "Pan model",       "Bacillota"
)

# ---------- Read data ----------
pan_aux      <- read_excel(xlsx_path, sheet = "pandraft_all_auxotrophies")
single_aux   <- read_excel(xlsx_path, sheet = "single_all_auxotrophies")
pan_sub      <- read_excel(xlsx_path, sheet = "pandraft_all_substrates")
single_sub   <- read_excel(xlsx_path, sheet = "single_all_substrates")
single_recon <- read_excel(xlsx_path, sheet = "single_all_reconstruction")
pan_stage1   <- read_excel(xlsx_path, sheet = "pandraft_all_stage1_summary")

# ---------- Build feature matrix (17 x 10) ----------
# Block 1: 3 auxotrophy categories (Other dropped)
aux_cat_levels <- c("Amino acid", "Vitamin/cofactor", "Metal")
aux_features <- bind_rows(
  pan_aux    %>% filter(organism %in% pan_targets) %>%
    select(organism, compound),
  single_aux %>% filter(organism %in% single_targets) %>%
    select(organism, compound)
) %>%
  mutate(category = aux_class(compound)) %>%
  filter(category %in% aux_cat_levels) %>%
  count(organism, category, name = "n") %>%
  pivot_wider(names_from = category, values_from = n, values_fill = 0)

# Ensure all 3 aux categories exist as columns
for (cat in aux_cat_levels) {
  if (!cat %in% names(aux_features)) aux_features[[cat]] <- 0L
}
aux_features <- aux_features %>% select(organism, all_of(aux_cat_levels))

# Block 2: 7 substrate categories (Other dropped)
sub_cat_levels <- c("Amino acids", "Mono- & disaccharides", "Amino sugars",
                    "Oligo- & polysaccharides", "Sugar alcohols & alcohols",
                    "Organic & sugar acids", "Nucleosides & bases")
sub_features <- bind_rows(
  pan_sub    %>% filter(organism %in% pan_targets) %>%
    select(organism, cpd_id),
  single_sub %>% filter(organism %in% single_targets) %>%
    select(organism, cpd_id)
) %>%
  mutate(category = compound_category[cpd_id]) %>%
  filter(category %in% sub_cat_levels) %>%
  count(organism, category, name = "n") %>%
  pivot_wider(names_from = category, values_from = n, values_fill = 0)

for (cat in sub_cat_levels) {
  if (!cat %in% names(sub_features)) sub_features[[cat]] <- 0L
}
sub_features <- sub_features %>% select(organism, all_of(sub_cat_levels))

# Combine into 17 x 10 matrix in metadata order
feature_mat_df <- all_meta %>%
  select(organism) %>%
  left_join(aux_features, by = "organism") %>%
  left_join(sub_features, by = "organism") %>%
  mutate(across(-organism, ~ifelse(is.na(.), 0, .)))

# Compose short feature names for the matrix
short_aux <- c("Amino acid" = "Aux:AA",
               "Vitamin/cofactor" = "Aux:Vit",
               "Metal" = "Aux:Metal")
short_sub <- c("Amino acids" = "Sub:AA",
               "Mono- & disaccharides" = "Sub:MonoDi",
               "Amino sugars" = "Sub:AminoSug",
               "Oligo- & polysaccharides" = "Sub:Poly",
               "Sugar alcohols & alcohols" = "Sub:Alcohol",
               "Organic & sugar acids" = "Sub:Acid",
               "Nucleosides & bases" = "Sub:Nucleo")

X <- feature_mat_df %>% select(-organism) %>% as.matrix()
rownames(X) <- feature_mat_df$organism
colnames(X) <- c(unname(short_aux[aux_cat_levels]),
                 unname(short_sub[sub_cat_levels]))

cat("=== Feature matrix (17 x 10) ===\n")
print(X)

# ---------- Sample weights from ORF coverage ----------
single_orf <- single_recon %>%
  filter(organism %in% single_targets) %>%
  select(organism, orf_cov = orf_coverage_pct)

pan_orf <- pan_stage1 %>%
  rename(organism = pan_id) %>%
  filter(organism %in% pan_targets) %>%
  select(organism, orf_cov = orf_cov_max)

orf_data <- bind_rows(single_orf, pan_orf)
orf_data$weight_raw  <- sqrt(orf_data$orf_cov)
n_samples            <- nrow(orf_data)
orf_data$weight      <- orf_data$weight_raw * n_samples / sum(orf_data$weight_raw)

sample_weights <- orf_data$weight[match(rownames(X), orf_data$organism)]

cat("\n=== Sample weights (sqrt-compressed, sum=n) ===\n")
print(data.frame(organism = rownames(X),
                 orf_cov  = round(orf_data$orf_cov[match(rownames(X), orf_data$organism)], 2),
                 weight   = round(sample_weights, 3)))

# ---------- Feature weights (block balancing) ----------
feature_weights <- c(rep(1/sqrt(length(aux_cat_levels)), length(aux_cat_levels)),
                     rep(1/sqrt(length(sub_cat_levels)), length(sub_cat_levels)))
names(feature_weights) <- colnames(X)
cat("\n=== Feature weights (block-balanced) ===\n")
print(round(feature_weights, 3))

# ---------- Run weighted PCA ----------
pca <- PCA(X,
           row.w        = sample_weights,
           col.w        = feature_weights,
           scale.unit   = TRUE,
           graph        = FALSE,
           ncp          = 5)

# Extract scores
scores <- as.data.frame(pca$ind$coord[, 1:3])
colnames(scores) <- c("PC1", "PC2", "PC3")
scores$organism <- rownames(X)
scores <- scores %>%
  left_join(all_meta, by = "organism") %>%
  mutate(group  = factor(group,  levels = names(group_pal)),
         phylum = factor(phylum, levels = names(phylum_pch)),
         color  = group_pal[as.character(group)],
         pch    = phylum_pch[as.character(phylum)])

# Variance explained
var_exp <- pca$eig[, 2]

cat(sprintf("\n=== Variance explained ===\n"))
cat(sprintf("PC1 = %.1f%%\nPC2 = %.1f%%\nPC3 = %.1f%%\nTotal (PC1-3) = %.1f%%\n",
            var_exp[1], var_exp[2], var_exp[3], sum(var_exp[1:3])))

# Save scores
write.csv(scores %>% select(organism, display, group, phylum, PC1, PC2, PC3),
          file.path(out_dir, "Fig3C_PCA_scores.csv"),
          row.names = FALSE)

# =============================================================
# Plot 1: 3D scatter via plot3D
# =============================================================

# Axis limits with 12% padding (so droplines reach a clean floor)
ax_pad <- function(v, frac = 0.12) {
  r <- range(v)
  pad <- frac * diff(r)
  c(r[1] - pad, r[2] + pad)
}
xlim_full <- ax_pad(scores$PC1)
ylim_full <- ax_pad(scores$PC2)
zlim_full <- ax_pad(scores$PC3)

out_3d <- file.path(out_dir, "Fig3C_PCA_3D.png")
png(out_3d, width = 11, height = 8, units = "in", res = 1000, bg = "white")

# Layout: scatter on left (~80%), legends stacked on right (~20%)
layout(matrix(c(1, 2), nrow = 1), widths = c(4.2, 1))

# ---- Region 1: 3D scatter ----
par(mar = c(2, 2, 3, 1), xpd = NA)

scatter3D(
  x         = scores$PC1, y = scores$PC2, z = scores$PC3,
  pch       = scores$pch,
  col       = scores$color,
  colvar    = NULL,
  cex       = 2.6,
  bty       = "g",
  ticktype  = "detailed",
  theta     = 35,
  phi       = 20,
  xlim      = xlim_full,
  ylim      = ylim_full,
  zlim      = zlim_full,
  xlab      = sprintf("PC1 (%.1f%%)", var_exp[1]),
  ylab      = sprintf("PC2 (%.1f%%)", var_exp[2]),
  zlab      = sprintf("PC3 (%.1f%%)", var_exp[3]),
  main      = "Weighted PCA on metabolic phenotype",
  cex.lab   = 1.15,
  cex.axis  = 0.85,
  cex.main  = 1.4,
  font.lab  = 2,
  font.main = 2
)

# Capture projection matrix BEFORE adding more elements
plist <- getplist()
pmat  <- plist$mat

# Droplines: from each point straight down to floor (PC3 = zlim_full[1])
n_pts   <- nrow(scores)
z_floor <- zlim_full[1]
segments3D(
  x0  = scores$PC1, y0 = scores$PC2, z0 = scores$PC3,
  x1  = scores$PC1, y1 = scores$PC2, z1 = rep(z_floor, n_pts),
  add = TRUE,
  col = "grey55",
  lty = 2,
  lwd = 0.6
)

# Project 3D points to 2D screen coordinates
screen <- trans3D(scores$PC1, scores$PC2, scores$PC3, pmat)
sx <- screen$x
sy <- screen$y

# Label positions: small initial offset to the right
lx <- sx + 0.025
ly <- sy

# Simple force-based repel in screen space (vertical pushes only)
min_dist <- 0.034
for (iter in seq_len(120)) {
  moved <- FALSE
  for (i in seq_along(lx)) {
    for (j in seq_along(lx)) {
      if (i >= j) next
      dx <- lx[i] - lx[j]
      dy <- ly[i] - ly[j]
      d  <- sqrt(dx * dx + dy * dy)
      if (d > 0 && d < min_dist) {
        push <- (min_dist - d) / 2
        if (dy >= 0) {
          ly[i] <- ly[i] + push
          ly[j] <- ly[j] - push
        } else {
          ly[i] <- ly[i] - push
          ly[j] <- ly[j] + push
        }
        moved <- TRUE
      }
    }
  }
  if (!moved) break
}

# Faint connector when label has moved away from its point
for (i in seq_along(lx)) {
  if (abs(ly[i] - sy[i]) > 0.005) {
    segments(sx[i], sy[i], lx[i] - 0.004, ly[i],
             col = "grey60", lwd = 0.35)
  }
}

# Draw labels in screen space (always on top, immune to 3D occlusion)
text(lx, ly,
     labels = scores$display,
     cex    = 0.55,
     font   = 2,
     col    = "grey15",
     adj    = c(0, 0.5))

# ---- Region 2: legends, stacked vertically ----
par(mar = c(0, 0, 0, 0))
plot.new()

# Two legends placed at explicit y positions (plot.new() coords run 0-1).
# Tweak LEGEND_Y_TOP and LEGEND_Y_BOT to move legends closer / further apart.
LEGEND_Y_TOP <- 0.62   # y-centre of Group legend
LEGEND_Y_BOT <- 0.38   # y-centre of Phylum legend

legend(x = 0.3, y = LEGEND_Y_TOP, xjust = 1, yjust = 0.5,
       title      = "Group",
       legend     = names(group_pal),
       fill       = unname(group_pal),
       border     = "white",
       cex        = 0.95,
       bty        = "n",
       text.font  = 2,
       title.font = 2,
       title.adj  = 0,
       y.intersp  = 1.3)

legend(x = 0.3, y = LEGEND_Y_BOT, xjust = 1, yjust = 0.5,
       title      = "Phylum",
       legend     = names(phylum_pch),
       pch        = unname(phylum_pch),
       cex        = 0.95,
       pt.cex     = 1.8,
       bty        = "n",
       text.font  = 2,
       title.font = 2,
       title.adj  = 0,
       y.intersp  = 1.3)

dev.off()
cat(sprintf("\nSaved: %s\n", out_3d))

# =============================================================
# Plot 2: scree + loadings biplots (ggplot)
# =============================================================
# Scree
scree_df <- data.frame(
  PC       = factor(paste0("PC", seq_len(nrow(pca$eig))),
                    levels = paste0("PC", seq_len(nrow(pca$eig)))),
  variance = pca$eig[, 2]
) %>% head(8)

p_scree <- ggplot(scree_df, aes(x = PC, y = variance)) +
  geom_col(fill = "#5D8AA8", colour = "grey20", width = 0.7) +
  geom_text(aes(label = sprintf("%.1f%%", variance)),
            vjust = -0.3, size = 3.4, fontface = "bold") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.18))) +
  labs(title = "Variance explained per PC", x = NULL, y = "% Variance") +
  theme_classic(base_size = 11) +
  theme(
    plot.title  = element_text(face = "bold", size = 12, hjust = 0,
                               margin = margin(b = 6)),
    axis.text   = element_text(face = "bold", colour = "grey15"),
    axis.title  = element_text(face = "bold")
  )

# Loadings (correlations of variables with PCs)
loadings_df <- as.data.frame(pca$var$coord[, 1:3])
colnames(loadings_df) <- c("PC1", "PC2", "PC3")
loadings_df$feature <- rownames(loadings_df)
loadings_df$block   <- ifelse(grepl("^Aux:", loadings_df$feature),
                              "Auxotrophy", "Substrate")

block_pal <- c(Auxotrophy = "#9B59B6", Substrate = "#2E7D32")

make_biplot <- function(df, pc_x, pc_y, title) {
  ggplot(df, aes(x = .data[[pc_x]], y = .data[[pc_y]], colour = block)) +
    geom_segment(aes(x = 0, y = 0, xend = .data[[pc_x]], yend = .data[[pc_y]]),
                 arrow = arrow(length = unit(0.18, "cm")),
                 linewidth = 0.7) +
    geom_text(aes(label = feature),
              size = 3.1, fontface = "bold",
              hjust = -0.08, vjust = -0.08,
              show.legend = FALSE) +
    geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey80") +
    geom_vline(xintercept = 0, linewidth = 0.3, colour = "grey80") +
    scale_colour_manual(values = block_pal, name = "Block") +
    expand_limits(x = c(-1.25, 1.25), y = c(-1.25, 1.25)) +
    coord_equal() +
    labs(title = title, x = pc_x, y = pc_y) +
    theme_classic(base_size = 11) +
    theme(
      plot.title    = element_text(face = "bold", size = 12,
                                   margin = margin(b = 6)),
      axis.text     = element_text(face = "bold", colour = "grey15"),
      axis.title    = element_text(face = "bold"),
      legend.title  = element_text(face = "bold"),
      legend.text   = element_text(face = "bold"),
      legend.position = "bottom"
    )
}

p_load12 <- make_biplot(loadings_df, "PC1", "PC2", "Loadings: PC1 vs PC2")
p_load13 <- make_biplot(loadings_df, "PC1", "PC3", "Loadings: PC1 vs PC3")

diag_fig <- p_scree / (p_load12 | p_load13) +
  plot_layout(heights = c(0.8, 1.6), guides = "collect") &
  theme(legend.position = "bottom")

out_diag <- file.path(out_dir, "Fig3C_PCA_diagnostic.png")
ggsave(out_diag, diag_fig,
       width = 11, height = 9.5, units = "in", dpi = 1000, bg = "white")

cat(sprintf("Saved: %s\n", out_diag))
cat(sprintf("Saved: %s\n", file.path(out_dir, "Fig3C_PCA_scores.csv")))
