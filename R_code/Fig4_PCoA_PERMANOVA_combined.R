# =============================================================================
# Fig4_PCoA_PERMANOVA_combined.R
# Manuscript 4 - Niche partitioning by PCoA + dataset-scope sensitivity sweep.
#
# Combined 4-panel figure:
#   A. Substrate utilisation PCoA (Bray-Curtis on Net_Growth).
#   B. Auxotrophy PCoA (Jaccard on binary essential factors).
#   C. PERMANOVA F-statistic across five dataset scopes.
#   D. Betadisper p-value across five dataset scopes.
#
# Layout:
#   - All four panel titles have identical font size and are centred.
#   - Panel tags (A/B/C/D) are auto-placed at the top-left of each panel
#     via patchwork's tag system.
#
# Output: figures/Fig4_PCoA_PERMANOVA_combined.png
# =============================================================================
suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(vegan)
  library(ape)
  library(ggplot2)
  library(patchwork)
  library(ggrepel)
  library(scales)
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

sub_pan <- sub_pan %>% filter(!is.na(compound))
sub_sin <- sub_sin %>% filter(!is.na(compound))
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

all_meta <- bind_rows(
  pan_meta %>% select(organism, Domain, Phylum, Class, role),
  single_meta_lookup
)

# ---- Build 42-organism matrices ----------------------------------------------
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
  sub_sin %>% select(organism, compound, Net_Growth)
)
sub_mat_42 <- build_mat(sub_combined, "Net_Growth", all42)

aux_combined <- bind_rows(
  aux_pan %>% select(organism, compound),
  aux_sin %>% select(organism, compound)
) %>% mutate(present = 1L)
aux_mat_42 <- build_mat(aux_combined, "present", all42, fill = 0L, binary = TRUE)

# ---- Distance + PCoA ---------------------------------------------------------
D_sub <- as.matrix(vegdist(sub_mat_42, method = "bray"))
D_aux <- as.matrix(vegdist(aux_mat_42, method = "jaccard", binary = TRUE))

pcoa_sub <- pcoa(as.dist(D_sub))
pcoa_aux <- pcoa(as.dist(D_aux))
var_sub  <- pcoa_sub$values$Relative_eig
var_aux  <- pcoa_aux$values$Relative_eig

pcs_sub <- as_tibble(pcoa_sub$vectors[, 1:3]) %>%
  setNames(c("PC1","PC2","PC3")) %>%
  mutate(organism = rownames(pcoa_sub$vectors)) %>%
  left_join(all_meta, by = "organism")
pcs_aux <- as_tibble(pcoa_aux$vectors[, 1:3]) %>%
  setNames(c("PC1","PC2","PC3")) %>%
  mutate(organism = rownames(pcoa_aux$vectors)) %>%
  left_join(all_meta, by = "organism")

# ---- PERMANOVA on the 42 set (for panel A/B subtitles) -----------------------
phy_meta <- all_meta %>% filter(organism %in% rownames(sub_mat_42)) %>%
  arrange(match(organism, rownames(sub_mat_42)))
perm_sub <- adonis2(D_sub ~ Phylum, data = phy_meta,
                    permutations = 999, by = "terms")
perm_aux <- adonis2(D_aux ~ Phylum, data = phy_meta,
                    permutations = 999, by = "terms")

# ---- v4 palette --------------------------------------------------------------
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

# ---- Plot helper for PCoA (panels A/B) ---------------------------------------
plot_pcoa <- function(pcs, var_pct, title_text, perm_F, perm_p) {
  pan_pts <- pcs %>% filter(role == "pan")
  ref_pts <- pcs %>% filter(role %in% c("ref_gapseq","ref_rumen")) %>%
    mutate(label = ref_label[organism])
  arg_pts <- pcs %>% filter(role == "arg_mag") %>%
    mutate(label = ifelse(organism == "MGYG000292883", "MAG-Lent", "MAG-Kirit"))

  ggplot() +
    geom_point(data = pan_pts,
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
                       labels = c(ref_gapseq = "Reference (gapseq)",
                                  ref_rumen  = "Reference (rumen)"),
                       name = "Role") +
    labs(title    = title_text,
         subtitle = sprintf("PERMANOVA phylum: F = %.2f, p = %.3f",
                            perm_F, perm_p),
         x = sprintf("PC1 (%.1f%%)", var_pct[1] * 100),
         y = sprintf("PC2 (%.1f%%)", var_pct[2] * 100))
}

p_a <- plot_pcoa(pcs_sub, var_sub,
                 "Substrate utilisation (Bray-Curtis on Net_Growth)",
                 perm_sub$F[1], perm_sub$`Pr(>F)`[1])
p_b <- plot_pcoa(pcs_aux, var_aux,
                 "Auxotrophy patterns (Jaccard on binary)",
                 perm_aux$F[1], perm_aux$`Pr(>F)`[1])

# ---- PERMANOVA sweep data (panels C/D) ---------------------------------------
perm <- tribble(
  ~scope,             ~n,  ~metric,       ~test_type,    ~stat,  ~p,    ~R2,
  "All 42",           42L, "Substrate",   "PERMANOVA",   3.509,  0.001, 0.376,
  "All 42",           42L, "Auxotrophy",  "PERMANOVA",   4.155,  0.002, 0.416,
  "All 42",           42L, "Substrate",   "betadisper",  4.864,  0.019, NA_real_,
  "All 42",           42L, "Auxotrophy",  "betadisper",  4.481,  0.003, NA_real_,
  "No ARG-MAGs",      40L, "Substrate",   "PERMANOVA",   3.884,  0.001, 0.307,
  "No ARG-MAGs",      40L, "Auxotrophy",  "PERMANOVA",   5.562,  0.001, 0.389,
  "No ARG-MAGs",      40L, "Substrate",   "betadisper",  4.858,  0.032, NA_real_,
  "No ARG-MAGs",      40L, "Auxotrophy",  "betadisper",  4.840,  0.005, NA_real_,
  "No clinical refs", 39L, "Substrate",   "PERMANOVA",   3.127,  0.001, 0.370,
  "No clinical refs", 39L, "Auxotrophy",  "PERMANOVA",   3.789,  0.002, 0.415,
  "No clinical refs", 39L, "Substrate",   "betadisper",  3.953,  0.062, NA_real_,
  "No clinical refs", 39L, "Auxotrophy",  "betadisper",  4.368,  0.016, NA_real_,
  "No refs",          36L, "Substrate",   "PERMANOVA",   3.495,  0.001, 0.368,
  "No refs",          36L, "Auxotrophy",  "PERMANOVA",   4.370,  0.001, 0.421,
  "No refs",          36L, "Substrate",   "betadisper",  3.775,  0.091, NA_real_,
  "No refs",          36L, "Auxotrophy",  "betadisper",  4.872,  0.018, NA_real_,
  "Pan only",         34L, "Substrate",   "PERMANOVA",   4.057,  0.001, 0.289,
  "Pan only",         34L, "Auxotrophy",  "PERMANOVA",   6.403,  0.001, 0.390,
  "Pan only",         34L, "Substrate",   "betadisper",  3.110,  0.067, NA_real_,
  "Pan only",         34L, "Auxotrophy",  "betadisper",  5.124,  0.022, NA_real_
)
scope_levels <- c("All 42", "No ARG-MAGs", "No clinical refs", "No refs", "Pan only")
n_levels     <- c(42, 40, 39, 36, 34)
perm <- perm %>%
  mutate(
    scope_label = factor(sprintf("%s\n(n = %d)", scope, n),
                         levels = sprintf("%s\n(n = %d)", scope_levels, n_levels)),
    metric      = factor(metric, levels = c("Substrate", "Auxotrophy"))
  )
metric_cols <- c("Substrate" = "#2A9D8F", "Auxotrophy" = "#9B5DE5")

# Panel C: PERMANOVA F
perm_F <- perm %>% filter(test_type == "PERMANOVA")
p_c <- ggplot(perm_F, aes(x = scope_label, y = stat, fill = metric)) +
  geom_col(position = position_dodge(width = 0.75), width = 0.7) +
  geom_text(aes(label = sprintf("F = %.2f\nR² = %.3f", stat, R2)),
            position = position_dodge(width = 0.75),
            vjust = -0.2, size = 2.7, lineheight = 0.9) +
  scale_fill_manual(values = metric_cols, name = NULL) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.22)),
                     limits = c(0, 8)) +
  labs(title = "PERMANOVA on phylum across five dataset scopes",
       subtitle = "All PERMANOVA p ≤ 0.002 (significant on every scope and metric).",
       x = NULL, y = "PERMANOVA F-statistic")

# Panel D: betadisper p
perm_disp <- perm %>% filter(test_type == "betadisper")
p_d <- ggplot(perm_disp, aes(x = scope_label, y = p, fill = metric)) +
  geom_col(position = position_dodge(width = 0.75), width = 0.7) +
  geom_hline(yintercept = 0.05, linetype = "dashed", colour = "#264653") +
  annotate("text", x = 6, y = 0.052,
           label = "p = 0.05", hjust = 1, vjust = 0,
           size = 3, colour = "#264653", fontface = "italic") +
  geom_text(aes(label = sprintf("p = %.3f", p)),
            position = position_dodge(width = 0.75),
            vjust = -0.4, size = 2.9) +
  scale_fill_manual(values = metric_cols, name = NULL) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.22)),
                     limits = c(0, 0.12),
                     breaks = c(0, 0.025, 0.05, 0.075, 0.10)) +
  labs(title = "Betadisper on phylum across five dataset scopes",
       subtitle = "Bars below the dashed line (p > 0.05) indicate homogeneous within-phylum dispersion.",
       x = NULL, y = "Betadisper p-value")

# ---- Compose 2 x 2 with patchwork tags ---------------------------------------
# - Tags A/B/C/D placed at top-left of each panel by patchwork.
# - Title is centred. Subtitle is centred and smaller.
# - All four panel titles have the same size.
TITLE_SIZE    <- 13
SUBTITLE_SIZE <- 10
TAG_SIZE      <- 18

final_fig <- (p_a + p_b) / (p_c + p_d) +
  plot_annotation(tag_levels = "A") +
  plot_layout(guides = "collect", heights = c(1, 1)) &
  theme_classic(base_size = 11) &
  theme(
    plot.title         = element_text(face = "bold", size = TITLE_SIZE, hjust = 0.5),
    plot.subtitle      = element_text(size = SUBTITLE_SIZE, colour = "grey30", hjust = 0.5),
    plot.tag           = element_text(face = "bold", size = TAG_SIZE, hjust = 0),
    plot.tag.position  = c(0.01, 0.99),
    legend.position    = "right",
    legend.key.size    = unit(0.4, "cm"),
    axis.text.x        = element_text(size = 9),
    panel.grid.major   = element_line(linewidth = 0.2, colour = "#eee")
  )

out_png <- file.path(out_dir, "Fig4_PCoA_PERMANOVA_combined.png")
ggsave(out_png, final_fig, width = 16, height = 13, dpi = 1000)
cat(sprintf(">>> Wrote %s\n", out_png))

# ---- Console summary ---------------------------------------------------------
cat(sprintf("\nSubstrate PERMANOVA (42): F = %.3f, R² = %.3f, p = %.3f\n",
            perm_sub$F[1], perm_sub$R2[1], perm_sub$`Pr(>F)`[1]))
cat(sprintf("Auxotrophy PERMANOVA (42): F = %.3f, R² = %.3f, p = %.3f\n",
            perm_aux$F[1], perm_aux$R2[1], perm_aux$`Pr(>F)`[1]))
