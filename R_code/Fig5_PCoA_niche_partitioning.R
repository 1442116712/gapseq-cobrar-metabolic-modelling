# =============================================================================
# Fig4_PCoA_niche_partitioning.R
# Manuscript 4 - Niche partitioning by PCoA + PERMANOVA
# 42 organisms in a single ordination (34 pan-Draft + 6 ref + 2 ARG-MAG).
# Family clusters overlaid:
#   n >= 3 families : stat_ellipse (norm, level = 0.68)
#   n = 2 families  : geom_segment between the two members
#   n = 1 families  : no overlay (single point only)
# Family names labelled at the cluster centroid in italics.
# PERMANOVA + beta-dispersion at phylum level.
# =============================================================================

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(vegan)        # vegdist, adonis2, betadisper
  library(ape)          # pcoa
  library(ggplot2)
  library(patchwork)
  library(ggrepel)
})

# ---- Inputs ------------------------------------------------------------------
xl      <- "C:/Users/CFL/OneDrive - Queen's University Belfast/PhD_onedrive/Manuscript4/M4_Supplement.xlsx"
out_dir <- "C:/Users/CFL/OneDrive - Queen's University Belfast/PhD_onedrive/Manuscript4/figures"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

set.seed(7)
SUB_THRESHOLD <- 0.01

# ---- Read data ---------------------------------------------------------------
sub_pan <- read_excel(xl, sheet = "pandraft_all_substrates")
aux_pan <- read_excel(xl, sheet = "pandraft_all_auxotrophies")
sub_sin <- read_excel(xl, sheet = "single_all_substrates")
aux_sin <- read_excel(xl, sheet = "single_all_auxotrophies")
tg      <- read_excel(xl, sheet = "target-SGBs")

keep_single <- c("Ecoli_K12_MG1655","Ecoli_PA3",
                 "Btheta_VPI5482","Btheta_KPPR3",
                 "Efaecalis_ATCC19433","Efaecalis_68A",
                 "MGYG000292883","MGYG000295553")
sub_sin <- sub_sin %>% filter(organism %in% keep_single)
aux_sin <- aux_sin %>% filter(organism %in% keep_single)

sub_pan <- sub_pan %>%
  mutate(Net_Growth = suppressWarnings(as.numeric(Net_Growth))) %>%
  filter(!is.na(compound))
sub_sin <- sub_sin %>%
  mutate(Net_Growth = suppressWarnings(as.numeric(Net_Growth))) %>%
  filter(!is.na(compound))
aux_pan <- aux_pan %>% filter(!is.na(compound))
aux_sin <- aux_sin %>% filter(!is.na(compound))

floor_thresh <- function(x) ifelse(is.na(x) | x < SUB_THRESHOLD, 0, x)
sub_pan <- sub_pan %>% mutate(Net_Growth = floor_thresh(Net_Growth))
sub_sin <- sub_sin %>% mutate(Net_Growth = floor_thresh(Net_Growth))

# ---- Metadata (with Family) -------------------------------------------------
sp_tax <- tg %>%
  group_by(Species_rep) %>%
  summarise(Domain = first(Domain),
            Phylum = first(Phylum),
            Class  = first(Class),
            Family = first(Family), .groups = "drop") %>%
  rename(sgb_id = Species_rep)

pan_orgs <- unique(sub_pan$organism)
pan_meta <- tibble(organism = pan_orgs) %>%
  mutate(sgb_id = sub("_pan$", "", organism)) %>%
  left_join(sp_tax, by = "sgb_id") %>%
  mutate(role = "pan")

single_meta_lookup <- tribble(
  ~organism,             ~Domain,    ~Phylum,             ~Class,                ~Family,             ~role,
  "Ecoli_K12_MG1655",    "Bacteria", "Pseudomonadota",    "Gammaproteobacteria", "Enterobacteriaceae", "ref_gapseq",
  "Ecoli_PA3",           "Bacteria", "Pseudomonadota",    "Gammaproteobacteria", "Enterobacteriaceae", "ref_rumen",
  "Btheta_VPI5482",      "Bacteria", "Bacteroidota",      "Bacteroidia",         "Bacteroidaceae",     "ref_gapseq",
  "Btheta_KPPR3",        "Bacteria", "Bacteroidota",      "Bacteroidia",         "Bacteroidaceae",     "ref_rumen",
  "Efaecalis_ATCC19433", "Bacteria", "Firmicutes",        "Bacilli",             "Enterococcaceae",    "ref_gapseq",
  "Efaecalis_68A",       "Bacteria", "Firmicutes",        "Bacilli",             "Enterococcaceae",    "ref_rumen",
  "MGYG000292883",       "Bacteria", "Lentisphaerota",    "Lentisphaeria",       "Lentisphaeraceae",   "arg_mag",
  "MGYG000295553",       "Bacteria", "Verrucomicrobiota", "Kiritimatiellia",     "Kiritimatiellaceae", "arg_mag")

all_meta <- bind_rows(
  pan_meta %>% select(organism, Domain, Phylum, Class, Family, role),
  single_meta_lookup
)

# ---- Build matrices ---------------------------------------------------------
build_mat <- function(df, value_col, organisms, fill = 0, binary = FALSE) {
  m <- df %>%
    select(organism, compound, !!value_col) %>%
    pivot_wider(names_from  = compound,
                values_from = !!value_col,
                values_fn   = if (binary) function(x) 1L else mean,
                values_fill = fill) %>%
    column_to_rownames("organism")
  m <- m[match(organisms, rownames(m)), , drop = FALSE]
  m[is.na(m)] <- fill
  as.matrix(m)
}

all42 <- c(pan_orgs, keep_single)

sub_combined <- bind_rows(
  sub_pan %>% select(organism, compound, Net_Growth),
  sub_sin %>% select(organism, compound, Net_Growth))
sub_mat_42   <- build_mat(sub_combined, "Net_Growth", all42)

aux_combined <- bind_rows(
  aux_pan %>% select(organism, compound),
  aux_sin %>% select(organism, compound)) %>% mutate(present = 1L)
aux_mat_42   <- build_mat(aux_combined, "present", all42, fill = 0L, binary = TRUE)

# ---- Distance + PCoA --------------------------------------------------------
D_sub <- as.matrix(vegdist(sub_mat_42, method = "bray"))
D_aux <- as.matrix(vegdist(aux_mat_42, method = "jaccard", binary = TRUE))

pcoa_sub <- pcoa(as.dist(D_sub))
pcoa_aux <- pcoa(as.dist(D_aux))

var_sub <- pcoa_sub$values$Relative_eig
var_aux <- pcoa_aux$values$Relative_eig

pcs_sub <- as_tibble(pcoa_sub$vectors[, 1:3]) %>%
  setNames(c("PC1","PC2","PC3")) %>%
  mutate(organism = rownames(pcoa_sub$vectors)) %>%
  left_join(all_meta, by = "organism")

pcs_aux <- as_tibble(pcoa_aux$vectors[, 1:3]) %>%
  setNames(c("PC1","PC2","PC3")) %>%
  mutate(organism = rownames(pcoa_aux$vectors)) %>%
  left_join(all_meta, by = "organism")

# ---- PERMANOVA + beta-disper -------------------------------------------------
phy_meta <- all_meta %>% filter(organism %in% rownames(sub_mat_42)) %>%
  arrange(match(organism, rownames(sub_mat_42)))

perm_sub <- adonis2(D_sub ~ Phylum, data = phy_meta,
                    permutations = 999, by = "terms")
perm_aux <- adonis2(D_aux ~ Phylum, data = phy_meta,
                    permutations = 999, by = "terms")
disp_sub <- betadisper(as.dist(D_sub), phy_meta$Phylum)
disp_aux <- betadisper(as.dist(D_aux), phy_meta$Phylum)
disp_sub_test <- permutest(disp_sub, permutations = 999)
disp_aux_test <- permutest(disp_aux, permutations = 999)

cat("\n=== PERMANOVA & dispersion ===\n")
cat(sprintf("Substrate ~ Phylum:  F = %.3f, R2 = %.3f, p = %.3f\n",
            perm_sub$F[1], perm_sub$R2[1], perm_sub$`Pr(>F)`[1]))
cat(sprintf("Auxotrophy ~ Phylum: F = %.3f, R2 = %.3f, p = %.3f\n",
            perm_aux$F[1], perm_aux$R2[1], perm_aux$`Pr(>F)`[1]))
cat(sprintf("Beta-dispersion (substrate, phylum):  F = %.3f, p = %.3f\n",
            disp_sub_test$tab$F[1], disp_sub_test$tab$`Pr(>F)`[1]))
cat(sprintf("Beta-dispersion (auxotrophy, phylum): F = %.3f, p = %.3f\n",
            disp_aux_test$tab$F[1], disp_aux_test$tab$`Pr(>F)`[1]))

# ---- Family layer prep (MST connections) ------------------------------------
# For each family with n >= 2 organisms, build the minimum spanning tree on
# Euclidean PC1-PC2 distances and return the (n-1) edges as line segments.
# n = 1 families contribute no segments.
mst_edges_for_family <- function(coords) {
  # coords: data frame with columns organism, PC1, PC2 (rows = family members)
  k <- nrow(coords)
  if (k < 2) return(NULL)
  if (k == 2) {
    return(tibble(x = coords$PC1[1], y = coords$PC2[1],
                  xend = coords$PC1[2], yend = coords$PC2[2]))
  }
  # Prim's algorithm on Euclidean distance
  d <- as.matrix(dist(coords[, c("PC1","PC2")], method = "euclidean"))
  in_tree <- rep(FALSE, k); in_tree[1] <- TRUE
  edges <- vector("list", k - 1)
  for (e in seq_len(k - 1)) {
    best_i <- NA_integer_; best_j <- NA_integer_; best_d <- Inf
    for (i in which(in_tree)) {
      for (j in which(!in_tree)) {
        if (d[i, j] < best_d) {
          best_d <- d[i, j]; best_i <- i; best_j <- j
        }
      }
    }
    edges[[e]] <- tibble(x = coords$PC1[best_i], y = coords$PC2[best_i],
                         xend = coords$PC1[best_j], yend = coords$PC2[best_j])
    in_tree[best_j] <- TRUE
  }
  bind_rows(edges)
}

prep_family_layers <- function(pcs) {
  fam_n <- pcs %>% count(Family, name = "n")
  fam_keep <- fam_n %>% filter(n >= 2) %>% pull(Family)
  
  seg_df <- pcs %>%
    filter(Family %in% fam_keep) %>%
    group_by(Family) %>%
    group_modify(~ mst_edges_for_family(.x)) %>%
    ungroup()
  
  lab_df <- pcs %>%
    filter(Family %in% fam_keep) %>%
    group_by(Family) %>%
    summarise(PC1 = mean(PC1), PC2 = mean(PC2), n = n(), .groups = "drop")
  
  list(segment = seg_df, label = lab_df,
       n_families = length(fam_keep),
       n_edges = nrow(seg_df))
}

fam_sub <- prep_family_layers(pcs_sub)
fam_aux <- prep_family_layers(pcs_aux)

cat(sprintf("\nFamily MST in substrate space: %d families, %d edges\n",
            fam_sub$n_families, fam_sub$n_edges))
cat(sprintf("Family MST in auxotrophy space: %d families, %d edges\n",
            fam_aux$n_families, fam_aux$n_edges))

# ---- Palette ----------------------------------------------------------------
class_palette <- c(
  "Bacteroidia"         = "#4F8FA8",
  "Clostridia"          = "#E07A1F",
  "Bacilli"             = "#F2CC4F",
  "Negativicutes"       = "#7A8550",
  "Methanobacteria"     = "#A04A85",
  "Coriobacteriia"      = "#9B7BAE",
  "Gammaproteobacteria" = "#6E5848",
  "Lentisphaeria"       = "#B22222",
  "Kiritimatiellia"     = "#B22222"
)

ref_label <- c(
  "Ecoli_K12_MG1655"    = "Ec K-12",
  "Ecoli_PA3"           = "Ec PA-3",
  "Btheta_VPI5482"      = "Bt VPI",
  "Btheta_KPPR3"        = "Bt KPPR",
  "Efaecalis_ATCC19433" = "Ef ATCC",
  "Efaecalis_68A"       = "Ef 68A"
)

# ---- Plot helper ------------------------------------------------------------
plot_pcoa <- function(pcs, fam_layers, var_pct, title,
                      perm_F, perm_p, disp_p) {
  pan_pts <- pcs %>% filter(role == "pan")
  ref_pts <- pcs %>% filter(role %in% c("ref_gapseq","ref_rumen")) %>%
    mutate(label = ref_label[organism])
  arg_pts <- pcs %>% filter(role == "arg_mag") %>%
    mutate(label = ifelse(organism == "MGYG000292883", "MAG-Lent", "MAG-Kirit"))
  
  ggplot() +
    # Family MST connections (n >= 2); behind everything
    geom_segment(data = fam_layers$segment,
                 aes(x = x, y = y, xend = xend, yend = yend),
                 colour = "grey45", linetype = "22",
                 linewidth = 0.45) +
    # Axis reference lines
    geom_hline(yintercept = 0, colour = "#bbb", linewidth = 0.3) +
    geom_vline(xintercept = 0, colour = "#bbb", linewidth = 0.3) +
    # Pan-Draft points
    geom_point(data = pan_pts,
               aes(PC1, PC2, fill = Class),
               shape = 21, size = 3.5, colour = "#333", stroke = 0.4,
               alpha = 0.85) +
    # Reference points
    geom_point(data = ref_pts,
               aes(PC1, PC2, fill = Class, shape = role),
               size = 4.2, colour = "#000", stroke = 1.0, alpha = 0.95) +
    geom_text_repel(data = ref_pts,
                    aes(PC1, PC2, label = label),
                    size = 2.6, colour = "#333", fontface = "italic",
                    seed = 7, max.overlaps = Inf, box.padding = 0.4) +
    # ARG-MAG points
    geom_point(data = arg_pts,
               aes(PC1, PC2),
               shape = 23, size = 5, fill = "#B22222", colour = "#000",
               stroke = 0.7) +
    geom_text_repel(data = arg_pts,
                    aes(PC1, PC2, label = label),
                    size = 3, colour = "#B22222", fontface = "bold",
                    seed = 7, max.overlaps = Inf, box.padding = 0.6) +
    # Family labels (on top, italic)
    geom_text_repel(data = fam_layers$label,
                    aes(PC1, PC2, label = Family),
                    fontface = "italic",
                    size = 2.8, colour = "grey25",
                    seed = 7, max.overlaps = Inf,
                    box.padding = 0.7, point.padding = 0.5,
                    segment.colour = "grey60",
                    segment.linetype = "dotted",
                    segment.size = 0.3) +
    scale_fill_manual(values = class_palette, name = "Class",
                      guide = guide_legend(override.aes = list(shape = 21))) +
    scale_shape_manual(values = c(ref_gapseq = 24, ref_rumen = 25),
                       labels = c(ref_gapseq = "Reference (gapseq)",
                                  ref_rumen  = "Reference (rumen)"),
                       name = "Role") +
    labs(title    = title,
         subtitle = sprintf("PERMANOVA phylum: F = %.2f, p = %.3f; betadisper p = %.3f",
                            perm_F, perm_p, disp_p),
         x = sprintf("PC1 (%.1f%%)", var_pct[1] * 100),
         y = sprintf("PC2 (%.1f%%)", var_pct[2] * 100)) +
    theme_classic(base_size = 10) +
    theme(plot.title       = element_text(size = 11, hjust = 0.5),
          plot.subtitle    = element_text(size = 8.5, colour = "#444", hjust = 0.5),
          legend.position  = "right",
          legend.key.size  = unit(0.4, "cm"),
          panel.grid.major = element_line(linewidth = 0.2, colour = "#eee"))
}

p_sub <- plot_pcoa(pcs_sub, fam_sub, var_sub,
                   "A   Substrate utilisation (Bray-Curtis on Net_Growth)",
                   perm_sub$F[1], perm_sub$`Pr(>F)`[1],
                   disp_sub_test$tab$`Pr(>F)`[1])
p_aux <- plot_pcoa(pcs_aux, fam_aux, var_aux,
                   "B   Auxotrophy patterns (Jaccard on binary)",
                   perm_aux$F[1], perm_aux$`Pr(>F)`[1],
                   disp_aux_test$tab$`Pr(>F)`[1])

fig4 <- p_sub + p_aux + plot_layout(ncol = 2, guides = "collect") &
  theme(legend.position = "right")

ggsave(file.path(out_dir, "Fig4_PCoA_niche_partitioning.png"),
       fig4, width = 14.5, height = 6.5, dpi = 1000)

cat("\nFigure 4 saved.\n")

# ---- Outputs ----------------------------------------------------------------
write.table(pcs_sub %>% mutate(matrix = "substrate"),
            file = file.path(out_dir, "Fig4_PCoA_coords.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)
write.table(pcs_aux %>% mutate(matrix = "auxotrophy"),
            file = file.path(out_dir, "Fig4_PCoA_coords.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE,
            append = TRUE, col.names = FALSE)
