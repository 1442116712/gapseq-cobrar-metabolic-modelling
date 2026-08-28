# =============================================================================
# SuppFig_PCoA_34_Gower.R
# Manuscript 4 - Supplementary niche partitioning analysis
# PCoA trained on 34 pan-Draft only; 6 references + 2 ARG-MAGs projected
# passively via Gower (1968) out-of-sample formula.
# Robustness check supporting the 42-organism main analysis (Fig 4).
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

# Force numeric and drop rows with missing compound
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

# ---- Metadata ----------------------------------------------------------------
sp_tax <- tg %>%
  group_by(Species_rep) %>%
  summarise(Domain = first(Domain),
            Phylum = first(Phylum),
            Class  = first(Class), .groups = "drop") %>%
  rename(sgb_id = Species_rep)

pan_orgs <- unique(sub_pan$organism)
pan_meta <- tibble(organism = pan_orgs) %>%
  mutate(sgb_id = sub("_pan$", "", organism)) %>%
  left_join(sp_tax, by = "sgb_id") %>%
  mutate(role = "pan")

single_meta_lookup <- tribble(
  ~organism,             ~Domain,    ~Phylum,             ~Class,                ~role,
  "Ecoli_K12_MG1655",    "Bacteria", "Pseudomonadota",    "Gammaproteobacteria", "ref_gapseq",
  "Ecoli_PA3",           "Bacteria", "Pseudomonadota",    "Gammaproteobacteria", "ref_rumen",
  "Btheta_VPI5482",      "Bacteria", "Bacteroidota",      "Bacteroidia",         "ref_gapseq",
  "Btheta_KPPR3",        "Bacteria", "Bacteroidota",      "Bacteroidia",         "ref_rumen",
  "Efaecalis_ATCC19433", "Bacteria", "Firmicutes",        "Bacilli",             "ref_gapseq",
  "Efaecalis_68A",       "Bacteria", "Firmicutes",        "Bacilli",             "ref_rumen",
  "MGYG000292883",       "Bacteria", "Lentisphaerota",    "Lentisphaeria",       "arg_mag",
  "MGYG000295553",       "Bacteria", "Verrucomicrobiota", "Kiritimatiellia",     "arg_mag")

# ---- Build matrices: train (34) vs new (8) -----------------------------------
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

sub_mat_pan <- build_mat(sub_pan, "Net_Growth", pan_orgs)
sub_mat_sin <- build_mat(sub_sin, "Net_Growth", keep_single)
miss <- setdiff(colnames(sub_mat_pan), colnames(sub_mat_sin))
if (length(miss) > 0) {
  add <- matrix(0, nrow = nrow(sub_mat_sin), ncol = length(miss),
                dimnames = list(rownames(sub_mat_sin), miss))
  sub_mat_sin <- cbind(sub_mat_sin, add)
}
sub_mat_sin <- sub_mat_sin[, colnames(sub_mat_pan), drop = FALSE]

aux_pan2 <- aux_pan %>% mutate(present = 1L)
aux_sin2 <- aux_sin %>% mutate(present = 1L)
aux_mat_pan <- build_mat(aux_pan2, "present", pan_orgs, fill = 0L, binary = TRUE)
aux_mat_sin <- build_mat(aux_sin2, "present", keep_single, fill = 0L, binary = TRUE)
miss <- setdiff(colnames(aux_mat_pan), colnames(aux_mat_sin))
if (length(miss) > 0) {
  add <- matrix(0L, nrow = nrow(aux_mat_sin), ncol = length(miss),
                dimnames = list(rownames(aux_mat_sin), miss))
  aux_mat_sin <- cbind(aux_mat_sin, add)
}
aux_mat_sin <- aux_mat_sin[, colnames(aux_mat_pan), drop = FALSE]

# ---- Train PCoA on 34 pan only -----------------------------------------------
D_sub_pan <- as.matrix(vegdist(sub_mat_pan, method = "bray"))
D_aux_pan <- as.matrix(vegdist(aux_mat_pan, method = "jaccard", binary = TRUE))

pcoa_sub <- pcoa(as.dist(D_sub_pan))
pcoa_aux <- pcoa(as.dist(D_aux_pan))

var_sub <- pcoa_sub$values$Relative_eig
var_aux <- pcoa_aux$values$Relative_eig

pan_pcs_sub <- as_tibble(pcoa_sub$vectors[, 1:3]) %>%
  setNames(c("PC1","PC2","PC3")) %>%
  mutate(organism = rownames(pcoa_sub$vectors)) %>%
  left_join(pan_meta, by = "organism")

pan_pcs_aux <- as_tibble(pcoa_aux$vectors[, 1:3]) %>%
  setNames(c("PC1","PC2","PC3")) %>%
  mutate(organism = rownames(pcoa_aux$vectors)) %>%
  left_join(pan_meta, by = "organism")

# ---- Gower (1968) out-of-sample projection -----------------------------------
project_gower <- function(mat_new, mat_train, method, binary = FALSE) {
  D_train <- as.matrix(vegdist(mat_train, method = method, binary = binary))
  D_full  <- as.matrix(vegdist(rbind(mat_train, mat_new),
                               method = method, binary = binary))
  n_train <- nrow(mat_train)
  D_cross <- D_full[(n_train + 1):nrow(D_full), 1:n_train, drop = FALSE]

  D2 <- D_train^2
  n  <- nrow(D2)
  J  <- diag(n) - matrix(1/n, n, n)
  B  <- -0.5 * J %*% D2 %*% J

  eig    <- eigen(B, symmetric = TRUE)
  ev     <- eig$values
  V      <- eig$vectors
  keep   <- which(ev > 1e-9)
  ev_k   <- ev[keep]
  V_k    <- V[, keep, drop = FALSE]

  row_means  <- rowMeans(D2)
  grand_mean <- mean(D2)
  D_cross2   <- D_cross^2
  col_new    <- rowMeans(D_cross2)
  a_new <- -0.5 * (D_cross2 -
                     matrix(row_means, nrow = nrow(D_cross2),
                            ncol = length(row_means), byrow = TRUE) -
                     col_new + grand_mean)

  pcs <- a_new %*% V_k[, 1:3, drop = FALSE] %*%
            diag(1 / sqrt(ev_k[1:3]))
  colnames(pcs) <- c("PC1", "PC2", "PC3")
  as_tibble(pcs, rownames = "organism")
}

single_pcs_sub <- project_gower(sub_mat_sin, sub_mat_pan, method = "bray") %>%
  left_join(single_meta_lookup, by = "organism")

single_pcs_aux <- project_gower(aux_mat_sin, aux_mat_pan,
                                method = "jaccard", binary = TRUE) %>%
  left_join(single_meta_lookup, by = "organism")

# ---- PERMANOVA + betadisper (on training set only) ---------------------------
phy_pan <- pan_meta %>% filter(organism %in% rownames(sub_mat_pan)) %>%
  arrange(match(organism, rownames(sub_mat_pan)))

perm_sub <- adonis2(D_sub_pan ~ Phylum, data = phy_pan,
                    permutations = 999, by = "terms")
perm_aux <- adonis2(D_aux_pan ~ Phylum, data = phy_pan,
                    permutations = 999, by = "terms")

# Homogeneity of multivariate dispersion (parallel to Fig 4 main analysis)
disp_sub <- betadisper(as.dist(D_sub_pan), phy_pan$Phylum)
disp_aux <- betadisper(as.dist(D_aux_pan), phy_pan$Phylum)
disp_sub_test <- permutest(disp_sub, permutations = 999)
disp_aux_test <- permutest(disp_aux, permutations = 999)

cat("\n=== Supplementary: 34-pan PCoA with Gower projection ===\n")
cat(sprintf("Substrate  ~ Phylum: F = %.3f, R2 = %.3f, p = %.3f\n",
            perm_sub$F[1], perm_sub$R2[1], perm_sub$`Pr(>F)`[1]))
cat(sprintf("Auxotrophy ~ Phylum: F = %.3f, R2 = %.3f, p = %.3f\n",
            perm_aux$F[1], perm_aux$R2[1], perm_aux$`Pr(>F)`[1]))
cat(sprintf("Beta-dispersion (substrate, phylum):  F = %.3f, p = %.3f\n",
            disp_sub_test$tab$F[1], disp_sub_test$tab$`Pr(>F)`[1]))
cat(sprintf("Beta-dispersion (auxotrophy, phylum): F = %.3f, p = %.3f\n",
            disp_aux_test$tab$F[1], disp_aux_test$tab$`Pr(>F)`[1]))

# ---- Palette (v4, matches Fig 4 main) ----------------------------------------
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

# ---- Plot helper -------------------------------------------------------------
plot_pcoa <- function(pan_pcs, single_pcs, var_pct, title,
                      perm_F, perm_p, disp_p) {
  ref_pts <- single_pcs %>% filter(role %in% c("ref_gapseq","ref_rumen")) %>%
    mutate(label = ref_label[organism])
  arg_pts <- single_pcs %>% filter(role == "arg_mag") %>%
    mutate(label = ifelse(organism == "MGYG000292883", "MAG-Lent", "MAG-Kirit"))

  ggplot() +
    geom_point(data = pan_pcs,
               aes(PC1, PC2, fill = Class),
               shape = 21, size = 3.5, colour = "#333", stroke = 0.4,
               alpha = 0.85) +
    geom_point(data = ref_pts,
               aes(PC1, PC2, fill = Class, shape = role),
               size = 4.2, colour = "#000", stroke = 1.0, alpha = 0.95) +
    geom_text_repel(data = ref_pts,
                    aes(PC1, PC2, label = label),
                    size = 2.6, colour = "#333", fontface = "italic",
                    seed = 7, max.overlaps = Inf, box.padding = 0.4) +
    geom_point(data = arg_pts,
               aes(PC1, PC2),
               shape = 23, size = 5, fill = "#B22222", colour = "#000",
               stroke = 0.7) +
    geom_text_repel(data = arg_pts,
                    aes(PC1, PC2, label = label),
                    size = 3, colour = "#B22222", fontface = "bold",
                    seed = 7, max.overlaps = Inf, box.padding = 0.6) +
    geom_hline(yintercept = 0, colour = "#bbb", linewidth = 0.3) +
    geom_vline(xintercept = 0, colour = "#bbb", linewidth = 0.3) +
    scale_fill_manual(values = class_palette, name = "Class",
                      guide = guide_legend(override.aes = list(shape = 21))) +
    scale_shape_manual(values = c(ref_gapseq = 24, ref_rumen = 25),
                       labels = c(ref_gapseq = "Reference (gapseq, projected)",
                                  ref_rumen  = "Reference (rumen, projected)"),
                       name = "Role") +
    labs(title    = title,
         subtitle = sprintf("PERMANOVA phylum (training): F = %.2f, p = %.3f; betadisper p = %.3f",
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

p_sub <- plot_pcoa(pan_pcs_sub, single_pcs_sub, var_sub,
                   "A   Substrate utilisation (34 pan training, Gower projection)",
                   perm_sub$F[1], perm_sub$`Pr(>F)`[1],
                   disp_sub_test$tab$`Pr(>F)`[1])
p_aux <- plot_pcoa(pan_pcs_aux, single_pcs_aux, var_aux,
                   "B   Auxotrophy patterns (34 pan training, Gower projection)",
                   perm_aux$F[1], perm_aux$`Pr(>F)`[1],
                   disp_aux_test$tab$`Pr(>F)`[1])

fig_supp <- p_sub + p_aux + plot_layout(ncol = 2, guides = "collect") &
  theme(legend.position = "right")

ggsave(file.path(out_dir, "SuppSupplementary Figure S2 (34 Gower).png"),
       fig_supp, width = 13.5, height = 6.5, dpi = 1000)

cat("\nSupplementary figure saved.\n")

write.table(bind_rows(
  pan_pcs_sub %>% mutate(matrix = "substrate", projected = "training"),
  single_pcs_sub %>% mutate(matrix = "substrate", projected = "Gower_1968"),
  pan_pcs_aux %>% mutate(matrix = "auxotrophy", projected = "training"),
  single_pcs_aux %>% mutate(matrix = "auxotrophy", projected = "Gower_1968")
), file = file.path(out_dir, "SuppFig_PCoA_34_Gower_coords.tsv"),
   sep = "\t", row.names = FALSE, quote = FALSE)
