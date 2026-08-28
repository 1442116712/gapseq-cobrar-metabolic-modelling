# =============================================================================
# SuppFig_PCoA_39_Gower.R
# Manuscript 4 - Supplementary niche partitioning analysis (variant 4)
# PCoA trained on 39 rumen-origin metabolic models:
#   34 pan-Draft + 2 ARG-MAG + 3 rumen-isolated reference.
# 3 clinical (gapseq default) reference models projected passively
# via Gower (1968).
# Directly tests whether clinical/laboratory cultured strains map to
# the same niche positions as their rumen-isolated counterparts;
# complementary to the paired clinical/rumen comparison in Fig 1.
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

# Three rumen-isolated references go into training
keep_ref_rumen    <- c("Ecoli_PA3","Btheta_KPPR3","Efaecalis_68A")
# Three clinical/gapseq references projected passively
keep_ref_clinical <- c("Ecoli_K12_MG1655","Btheta_VPI5482","Efaecalis_ATCC19433")
keep_argmag       <- c("MGYG000292883","MGYG000295553")

sub_ref_rumen    <- sub_sin %>% filter(organism %in% keep_ref_rumen)
aux_ref_rumen    <- aux_sin %>% filter(organism %in% keep_ref_rumen)
sub_ref_clinical <- sub_sin %>% filter(organism %in% keep_ref_clinical)
aux_ref_clinical <- aux_sin %>% filter(organism %in% keep_ref_clinical)
sub_arg          <- sub_sin %>% filter(organism %in% keep_argmag)
aux_arg          <- aux_sin %>% filter(organism %in% keep_argmag)

# Force numeric and drop missing-compound rows
force_num_filter <- function(df, has_growth = TRUE) {
  df <- df %>% filter(!is.na(compound))
  if (has_growth) {
    df <- df %>% mutate(Net_Growth = suppressWarnings(as.numeric(Net_Growth)))
  }
  df
}
sub_pan          <- force_num_filter(sub_pan)
sub_ref_rumen    <- force_num_filter(sub_ref_rumen)
sub_ref_clinical <- force_num_filter(sub_ref_clinical)
sub_arg          <- force_num_filter(sub_arg)
aux_pan          <- force_num_filter(aux_pan, has_growth = FALSE)
aux_ref_rumen    <- force_num_filter(aux_ref_rumen, has_growth = FALSE)
aux_ref_clinical <- force_num_filter(aux_ref_clinical, has_growth = FALSE)
aux_arg          <- force_num_filter(aux_arg, has_growth = FALSE)

floor_thresh <- function(x) ifelse(is.na(x) | x < SUB_THRESHOLD, 0, x)
sub_pan          <- sub_pan          %>% mutate(Net_Growth = floor_thresh(Net_Growth))
sub_ref_rumen    <- sub_ref_rumen    %>% mutate(Net_Growth = floor_thresh(Net_Growth))
sub_ref_clinical <- sub_ref_clinical %>% mutate(Net_Growth = floor_thresh(Net_Growth))
sub_arg          <- sub_arg          %>% mutate(Net_Growth = floor_thresh(Net_Growth))

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

ref_rumen_meta <- tribble(
  ~organism,        ~Domain,    ~Phylum,          ~Class,                ~role,
  "Ecoli_PA3",      "Bacteria", "Pseudomonadota", "Gammaproteobacteria", "ref_rumen",
  "Btheta_KPPR3",   "Bacteria", "Bacteroidota",   "Bacteroidia",         "ref_rumen",
  "Efaecalis_68A",  "Bacteria", "Firmicutes",     "Bacilli",             "ref_rumen")

ref_clinical_meta <- tribble(
  ~organism,             ~Domain,    ~Phylum,          ~Class,                ~role,
  "Ecoli_K12_MG1655",    "Bacteria", "Pseudomonadota", "Gammaproteobacteria", "ref_gapseq",
  "Btheta_VPI5482",      "Bacteria", "Bacteroidota",   "Bacteroidia",         "ref_gapseq",
  "Efaecalis_ATCC19433", "Bacteria", "Firmicutes",     "Bacilli",             "ref_gapseq")

arg_meta <- tribble(
  ~organism,         ~Domain,    ~Phylum,             ~Class,            ~role,
  "MGYG000292883",   "Bacteria", "Lentisphaerota",    "Lentisphaeria",   "arg_mag",
  "MGYG000295553",   "Bacteria", "Verrucomicrobiota", "Kiritimatiellia", "arg_mag")

train_meta <- bind_rows(
  pan_meta %>% select(organism, Domain, Phylum, Class, role),
  ref_rumen_meta,
  arg_meta
)

# ---- Build matrices ----------------------------------------------------------
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

align_cols <- function(mat_new, ref_cols, fill = 0) {
  miss <- setdiff(ref_cols, colnames(mat_new))
  if (length(miss) > 0) {
    add <- matrix(fill, nrow = nrow(mat_new), ncol = length(miss),
                  dimnames = list(rownames(mat_new), miss))
    mat_new <- cbind(mat_new, add)
  }
  mat_new[, ref_cols, drop = FALSE]
}

# 39 training matrix (34 pan + 3 rumen ref + 2 ARG-MAG): substrates
sub_train_data <- bind_rows(
  sub_pan        %>% select(organism, compound, Net_Growth),
  sub_ref_rumen  %>% select(organism, compound, Net_Growth),
  sub_arg        %>% select(organism, compound, Net_Growth)
)
train_orgs <- c(pan_orgs, keep_ref_rumen, keep_argmag)
sub_mat_train <- build_mat(sub_train_data, "Net_Growth", train_orgs)

# 3 clinical reference matrix to project
sub_mat_clinical <- build_mat(sub_ref_clinical, "Net_Growth", keep_ref_clinical)
sub_mat_clinical <- align_cols(sub_mat_clinical, colnames(sub_mat_train))

# 39 training matrix: auxotrophies (binary)
aux_train_data <- bind_rows(
  aux_pan       %>% select(organism, compound),
  aux_ref_rumen %>% select(organism, compound),
  aux_arg       %>% select(organism, compound)
) %>% mutate(present = 1L)
aux_mat_train <- build_mat(aux_train_data, "present", train_orgs,
                            fill = 0L, binary = TRUE)

aux_clinical2    <- aux_ref_clinical %>% mutate(present = 1L)
aux_mat_clinical <- build_mat(aux_clinical2, "present", keep_ref_clinical,
                               fill = 0L, binary = TRUE)
aux_mat_clinical <- align_cols(aux_mat_clinical, colnames(aux_mat_train), fill = 0L)

cat(sprintf("Training set: %d organisms x %d substrate compounds\n",
            nrow(sub_mat_train), ncol(sub_mat_train)))
cat(sprintf("Training set: %d organisms x %d auxotrophy compounds\n",
            nrow(aux_mat_train), ncol(aux_mat_train)))
cat(sprintf("Clinical reference set: %d organisms (projected)\n",
            nrow(sub_mat_clinical)))

# ---- Train PCoA on 39 -------------------------------------------------------
D_sub_train <- as.matrix(vegdist(sub_mat_train, method = "bray"))
D_aux_train <- as.matrix(vegdist(aux_mat_train, method = "jaccard", binary = TRUE))

pcoa_sub <- pcoa(as.dist(D_sub_train))
pcoa_aux <- pcoa(as.dist(D_aux_train))

var_sub <- pcoa_sub$values$Relative_eig
var_aux <- pcoa_aux$values$Relative_eig

train_pcs_sub <- as_tibble(pcoa_sub$vectors[, 1:3]) %>%
  setNames(c("PC1","PC2","PC3")) %>%
  mutate(organism = rownames(pcoa_sub$vectors)) %>%
  left_join(train_meta, by = "organism")

train_pcs_aux <- as_tibble(pcoa_aux$vectors[, 1:3]) %>%
  setNames(c("PC1","PC2","PC3")) %>%
  mutate(organism = rownames(pcoa_aux$vectors)) %>%
  left_join(train_meta, by = "organism")

# ---- Gower (1968) projection for 3 clinical refs ----------------------------
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

clinical_pcs_sub <- project_gower(sub_mat_clinical, sub_mat_train,
                                  method = "bray") %>%
  left_join(ref_clinical_meta, by = "organism")

clinical_pcs_aux <- project_gower(aux_mat_clinical, aux_mat_train,
                                  method = "jaccard", binary = TRUE) %>%
  left_join(ref_clinical_meta, by = "organism")

# ---- PERMANOVA + betadisper on 39 training set -------------------------------
phy_train <- train_meta %>%
  filter(organism %in% rownames(sub_mat_train)) %>%
  arrange(match(organism, rownames(sub_mat_train)))

perm_sub <- adonis2(D_sub_train ~ Phylum, data = phy_train,
                    permutations = 999, by = "terms")
perm_aux <- adonis2(D_aux_train ~ Phylum, data = phy_train,
                    permutations = 999, by = "terms")

disp_sub <- betadisper(as.dist(D_sub_train), phy_train$Phylum)
disp_aux <- betadisper(as.dist(D_aux_train), phy_train$Phylum)
disp_sub_test <- permutest(disp_sub, permutations = 999)
disp_aux_test <- permutest(disp_aux, permutations = 999)

cat("\n=== Supplementary (39-train rumen-origin + 3 clinical projection) ===\n")
cat(sprintf("Substrate  ~ Phylum: F = %.3f, R2 = %.3f, p = %.3f\n",
            perm_sub$F[1], perm_sub$R2[1], perm_sub$`Pr(>F)`[1]))
cat(sprintf("Auxotrophy ~ Phylum: F = %.3f, R2 = %.3f, p = %.3f\n",
            perm_aux$F[1], perm_aux$R2[1], perm_aux$`Pr(>F)`[1]))
cat(sprintf("Beta-dispersion (substrate, phylum):  F = %.3f, p = %.3f\n",
            disp_sub_test$tab$F[1], disp_sub_test$tab$`Pr(>F)`[1]))
cat(sprintf("Beta-dispersion (auxotrophy, phylum): F = %.3f, p = %.3f\n",
            disp_aux_test$tab$F[1], disp_aux_test$tab$`Pr(>F)`[1]))

# ---- Concordance metric: distance between clinical and rumen ref pairs ------
# In the trained space, each clinical ref (projected) has a counterpart
# rumen ref (training). Compute Euclidean distance in PC1-PC2 space.
cat("\n=== Clinical (projected) vs paired rumen (training) PC1-PC2 distance ===\n")
pairs <- list(
  c("Ecoli_K12_MG1655",    "Ecoli_PA3",     "E. coli"),
  c("Btheta_VPI5482",      "Btheta_KPPR3",  "B. thetaiotaomicron"),
  c("Efaecalis_ATCC19433", "Efaecalis_68A", "E. faecalis"))
for (pr in pairs) {
  for (mat_lbl in c("substrate","auxotrophy")) {
    tbl_clin  <- if (mat_lbl == "substrate") clinical_pcs_sub else clinical_pcs_aux
    tbl_train <- if (mat_lbl == "substrate") train_pcs_sub    else train_pcs_aux
    a <- tbl_clin  %>% filter(organism == pr[1])
    b <- tbl_train %>% filter(organism == pr[2])
    if (nrow(a) == 1 && nrow(b) == 1) {
      d <- sqrt((a$PC1 - b$PC1)^2 + (a$PC2 - b$PC2)^2)
      cat(sprintf("  %s (%s): PC1-PC2 distance = %.4f\n", pr[3], mat_lbl, d))
    }
  }
}

# ---- Palette (v4) -----------------------------------------------------------
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
plot_pcoa <- function(train_pcs, clinical_pcs, var_pct, title,
                      perm_F, perm_p, disp_p) {
  pan_pts <- train_pcs %>% filter(role == "pan")
  rumen_pts <- train_pcs %>% filter(role == "ref_rumen") %>%
    mutate(label = ref_label[organism])
  arg_pts <- train_pcs %>% filter(role == "arg_mag") %>%
    mutate(label = ifelse(organism == "MGYG000292883", "MAG-Lent", "MAG-Kirit"))
  clin <- clinical_pcs %>%
    mutate(label = ref_label[organism])

  # Pair-connecting segments (rumen training -> clinical projected)
  seg_data <- bind_rows(
    lapply(list(c("Ecoli_K12_MG1655","Ecoli_PA3"),
                c("Btheta_VPI5482","Btheta_KPPR3"),
                c("Efaecalis_ATCC19433","Efaecalis_68A")),
           function(pr) {
             a <- clinical_pcs %>% filter(organism == pr[1])
             b <- train_pcs    %>% filter(organism == pr[2])
             if (nrow(a) == 1 && nrow(b) == 1) {
               tibble(x = a$PC1, y = a$PC2, xend = b$PC1, yend = b$PC2)
             } else NULL
           }))

  ggplot() +
    # Pair connectors (drawn first, behind points)
    geom_segment(data = seg_data,
                 aes(x = x, y = y, xend = xend, yend = yend),
                 colour = "#999", linewidth = 0.4, linetype = "22") +
    # Pan-Draft training
    geom_point(data = pan_pts,
               aes(PC1, PC2, fill = Class),
               shape = 21, size = 3.5, colour = "#333", stroke = 0.4,
               alpha = 0.85) +
    # Rumen reference training (real coords, inverted triangle, filled)
    geom_point(data = rumen_pts,
               aes(PC1, PC2, fill = Class),
               shape = 25, size = 4.2, colour = "#000", stroke = 1.0,
               alpha = 0.95) +
    geom_text_repel(data = rumen_pts,
                    aes(PC1, PC2, label = label),
                    size = 2.6, colour = "#333", fontface = "italic",
                    seed = 7, max.overlaps = Inf, box.padding = 0.4) +
    # ARG-MAG training (real coords, red filled diamond)
    geom_point(data = arg_pts,
               aes(PC1, PC2),
               shape = 23, size = 5, fill = "#B22222", colour = "#000",
               stroke = 0.7) +
    geom_text_repel(data = arg_pts,
                    aes(PC1, PC2, label = label),
                    size = 3, colour = "#B22222", fontface = "bold",
                    seed = 7, max.overlaps = Inf, box.padding = 0.6) +
    # Clinical reference projected (upward triangle, filled)
    geom_point(data = clin,
               aes(PC1, PC2, fill = Class),
               shape = 24, size = 4.2, colour = "#000", stroke = 1.0,
               alpha = 0.95) +
    geom_text_repel(data = clin,
                    aes(PC1, PC2, label = label),
                    size = 2.6, colour = "#333", fontface = "italic",
                    seed = 7, max.overlaps = Inf, box.padding = 0.4) +
    geom_hline(yintercept = 0, colour = "#bbb", linewidth = 0.3) +
    geom_vline(xintercept = 0, colour = "#bbb", linewidth = 0.3) +
    scale_fill_manual(values = class_palette, name = "Class",
                      guide = guide_legend(override.aes = list(shape = 21))) +
    labs(title    = title,
         subtitle = sprintf("PERMANOVA phylum (39 training): F = %.2f, p = %.3f; betadisper p = %.3f",
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

p_sub <- plot_pcoa(train_pcs_sub, clinical_pcs_sub, var_sub,
                   "A   Substrate utilisation (39 rumen-origin training, clinical refs projected)",
                   perm_sub$F[1], perm_sub$`Pr(>F)`[1],
                   disp_sub_test$tab$`Pr(>F)`[1])
p_aux <- plot_pcoa(train_pcs_aux, clinical_pcs_aux, var_aux,
                   "B   Auxotrophy patterns (39 rumen-origin training, clinical refs projected)",
                   perm_aux$F[1], perm_aux$`Pr(>F)`[1],
                   disp_aux_test$tab$`Pr(>F)`[1])

fig_supp <- p_sub + p_aux + plot_layout(ncol = 2, guides = "collect") &
  theme(legend.position = "right")

ggsave(file.path(out_dir, "SuppSupplementary Figure S2 (39 Gower).png"),
       fig_supp, width = 13.5, height = 6.5, dpi = 1000)

cat("\nFigure saved.\n")

# ---- Outputs -----------------------------------------------------------------
write.table(bind_rows(
  train_pcs_sub    %>% mutate(matrix = "substrate",  projected = "training"),
  clinical_pcs_sub %>% mutate(matrix = "substrate",  projected = "Gower_1968"),
  train_pcs_aux    %>% mutate(matrix = "auxotrophy", projected = "training"),
  clinical_pcs_aux %>% mutate(matrix = "auxotrophy", projected = "Gower_1968")
), file = file.path(out_dir, "SuppFig_PCoA_39_Gower_coords.tsv"),
   sep = "\t", row.names = FALSE, quote = FALSE)
