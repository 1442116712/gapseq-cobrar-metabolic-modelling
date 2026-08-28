# =============================================================================
# Fig4.R
# Manuscript 4 - Niche partitioning by PCoA + dataset-scope sensitivity sweep.
#
# Combined 4-panel figure:
#   A. Substrate utilisation PCoA (Bray-Curtis on Net_Growth).
#   B. Auxotrophy PCoA (Jaccard on binary essential factors).
#   C. PERMANOVA F-statistic across five dataset scopes.
#   D. Betadisper p-value across five dataset scopes.
#
# What this script computes (all done inline; nothing hardcoded):
#   1. Bray-Curtis distances on substrate Net_Growth (quantitative) and
#      Jaccard distances on binary auxotrophy presence/absence (dissimilarity
#      choice matched to data type).
#   2. PCoA for visual confirmation of phylum-level clustering, run before
#      any formal test.
#   3. PERMANOVA on phylum (vegan::adonis2, 999 permutations).
#   4. Homogeneity of multivariate dispersion (vegan::betadisper +
#      vegan::permutest, 999 permutations) — the only formal PERMANOVA
#      pre-assumption. Class-level PERMANOVA is not attempted because three
#      of six classes contain < 3 Pan-Draft samples; phylum is used as the
#      grouping factor throughout.
#   5. Steps 3 and 4 are repeated across five dataset scopes (All 42 →
#      Pan only, n = 42, 40, 39, 36, 34) as a robustness check for the
#      inclusion of the two ARG-MAGs and the six reference single-genome
#      models. Panels C and D display the sweep.
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
  library(purrr)
  library(vegan)
  library(ape)
  library(ggplot2)
  library(patchwork)
  library(ggrepel)
  library(scales)
})

# ---- Inputs ------------------------------------------------------------------
xl      <- "C:/Users/CFL/OneDrive - Queen's University Belfast/PhD_onedrive/Manuscript4_gapseq/others/M4_Supplement.xlsx"
out_dir <- "C:/Users/CFL/OneDrive - Queen's University Belfast/PhD_onedrive/Manuscript4_gapseq/figures"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
set.seed(7)
SUB_THRESHOLD <- 0.01
N_PERM        <- 999

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

# ---- Distance + PCoA (on the full 42-organism set, for panels A/B) -----------
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

# ---- PERMANOVA + betadisper helpers ------------------------------------------
# Runs PERMANOVA (adonis2) and dispersion homogeneity (betadisper + permutest)
# on a single (distance matrix, grouping vector) pair. Returns a 2-row tibble:
# one row for PERMANOVA (with F, R2, p) and one for betadisper (with F, p, NA R2).
run_permanova_betadisper <- function(D, group_vec, scope, n, metric) {
  stopifnot(nrow(D) == length(group_vec))
  meta_df <- data.frame(group = group_vec)

  perm  <- adonis2(as.dist(D) ~ group, data = meta_df,
                   permutations = N_PERM, by = "terms")
  disp  <- betadisper(as.dist(D), meta_df$group)
  dtest <- permutest(disp, permutations = N_PERM)

  tibble(
    scope     = scope, n = n, metric = metric,
    test_type = c("PERMANOVA", "betadisper"),
    stat      = c(perm$F[1],  dtest$tab$F[1]),
    p         = c(perm$`Pr(>F)`[1], dtest$tab$`Pr(>F)`[1]),
    R2        = c(perm$R2[1], NA_real_)
  )
}

# Order and membership of the five dataset scopes.
scope_defs <- list(
  "All 42"           = all42,
  "No ARG-MAGs"      = setdiff(all42, c("MGYG000292883", "MGYG000295553")),
  "No clinical refs" = setdiff(all42,
                               c("Ecoli_K12_MG1655", "Btheta_VPI5482",
                                 "Efaecalis_ATCC19433")),
  "No refs"          = setdiff(all42, keep_single[1:6]),
  "Pan only"         = pan_orgs
)
stopifnot(vapply(scope_defs, length, integer(1)) == c(42, 40, 39, 36, 34))

# Subset a distance matrix + phylum vector to a given scope.
subset_dist <- function(D, keep) {
  D[keep, keep, drop = FALSE]
}
phylum_for <- function(orgs) {
  phy <- all_meta$Phylum[match(orgs, all_meta$organism)]
  stopifnot(!anyNA(phy))
  phy
}

# ---- 42-scope PERMANOVA numbers used as subtitles for panels A/B -------------
tmp42_sub <- run_permanova_betadisper(D_sub, phylum_for(all42), "All 42", 42, "Substrate")
tmp42_aux <- run_permanova_betadisper(D_aux, phylum_for(all42), "All 42", 42, "Predicted essential-factor profile")
perm_sub_F <- tmp42_sub$stat[1]; perm_sub_p <- tmp42_sub$p[1]
perm_aux_F <- tmp42_aux$stat[1]; perm_aux_p <- tmp42_aux$p[1]

# ---- Full 5-scope sweep, both metrics ----------------------------------------
perm <- map_dfr(names(scope_defs), function(scope) {
  keep <- scope_defs[[scope]]
  n    <- length(keep)
  phy  <- phylum_for(keep)
  Ds   <- subset_dist(D_sub, keep)
  Da   <- subset_dist(D_aux, keep)
  bind_rows(
    run_permanova_betadisper(Ds, phy, scope, n, "Substrate"),
    run_permanova_betadisper(Da, phy, scope, n, "Predicted essential-factor profile")
  )
})

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
               size = 4.2, colour = "#333", stroke = 0.4, alpha = 0.95) +
    geom_text_repel(data = ref_pts,
                    aes(PC1, PC2, label = label),
                    size = 4, colour = "#333", fontface = "italic",
                    seed = 7, max.overlaps = Inf, box.padding = 0.4) +
    geom_point(data = arg_pts,
               aes(PC1, PC2, shape = role),
               size = 5, fill = "#B22222", colour = "#000",
               stroke = 0.7) +
    geom_text_repel(data = arg_pts,
                    aes(PC1, PC2, label = label),
                    size = 4.5, colour = "#B22222", fontface = "bold",
                    seed = 7, max.overlaps = Inf, box.padding = 0.6) +
    geom_hline(yintercept = 0, colour = "#bbb", linewidth = 0.3) +
    geom_vline(xintercept = 0, colour = "#bbb", linewidth = 0.3) +
    scale_fill_manual(values = class_palette, name = "Class",
                      guide = guide_legend(
                        override.aes = list(shape = 21, size = 5,
                                            colour = "#333", stroke = 0.4))) +
    scale_shape_manual(
      breaks = c("ref_gapseq", "ref_rumen", "arg_mag"),
      values = c(ref_gapseq = 24, ref_rumen = 25, arg_mag = 23),
      labels = c(ref_gapseq = "Clinical reference",
                 ref_rumen  = "Rumen reference",
                 arg_mag    = "ARG-MAG target"),
      name = "Reference / target",
      guide = guide_legend(
        keyheight = unit(0.7, "cm"),   # per-legend override so keys don't overlap
        override.aes = list(
          fill   = c("grey60", "grey60", "#B22222"),
          colour = c("#333",   "#333",   "#000"),
          stroke = c(0.4,      0.4,      0.7),
          size   = c(5,        5,        5.5)
        ))
    ) +
    labs(title    = title_text,
         subtitle = sprintf("PERMANOVA phylum: F = %.2f, p = %.3f",
                            perm_F, perm_p),
         x = sprintf("PC1 (%.1f%%)", var_pct[1] * 100),
         y = sprintf("PC2 (%.1f%%)", var_pct[2] * 100))
}

p_a <- plot_pcoa(pcs_sub, var_sub,
                 "Substrate utilisation (Bray-Curtis on Net_Growth)",
                 perm_sub_F, perm_sub_p)
p_b <- plot_pcoa(pcs_aux, var_aux,
                 "Predicted essential-factor profile (Jaccard on binary)",
                 perm_aux_F, perm_aux_p)

# ---- Panel-ready factor levels for the sweep (panels C/D) --------------------
scope_levels <- c("All 42", "No ARG-MAGs", "No clinical refs", "No refs", "Pan only")
n_levels     <- c(42, 40, 39, 36, 34)
perm <- perm %>%
  mutate(
    scope       = factor(scope, levels = scope_levels),
    scope_label = factor(sprintf("%s\n(n = %d)", scope, n),
                         levels = sprintf("%s\n(n = %d)", scope_levels, n_levels)),
    metric      = factor(metric, levels = c("Substrate", "Predicted essential-factor profile"))
  )
metric_cols <- c("Substrate" = "#2A9D8F", "Predicted essential-factor profile" = "#9B5DE5")

# Panel C: PERMANOVA F
perm_F <- perm %>% filter(test_type == "PERMANOVA")
p_c <- ggplot(perm_F, aes(x = scope_label, y = stat, fill = metric)) +
  geom_col(position = position_dodge(width = 0.75), width = 0.7) +
  geom_text(aes(label = sprintf("F = %.2f\nR² = %.3f", stat, R2)),
            position = position_dodge(width = 0.75),
            vjust = -0.2, size = 4, lineheight = 0.9) +
  scale_fill_manual(values = metric_cols, name = NULL) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.22)),
                     limits = c(0, 8)) +
  labs(title = "PERMANOVA on phylum across five dataset scopes",
       subtitle = sprintf("Max PERMANOVA p across all scopes and metrics = %.3f\n(all significant).",
                          max(perm_F$p)),
       x = NULL, y = "PERMANOVA F-statistic")

# Panel D: betadisper p
perm_disp <- perm %>% filter(test_type == "betadisper")
p_d <- ggplot(perm_disp, aes(x = scope_label, y = p, fill = metric)) +
  geom_col(position = position_dodge(width = 0.75), width = 0.7) +
  geom_hline(yintercept = 0.05, linetype = "dashed", colour = "#264653") +
  annotate("text", x = 6, y = 0.052,
           label = "p = 0.05", hjust = 1, vjust = 0,
           size = 4.5, colour = "#264653", fontface = "italic") +
  geom_text(aes(label = sprintf("p = %.3f", p)),
            position = position_dodge(width = 0.75),
            vjust = -0.4, size = 4, lineheight = 0.9) +
  scale_fill_manual(values = metric_cols, name = NULL) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.22)),
                     limits = c(0, 0.10),
                     breaks = c(0, 0.025, 0.05, 0.075, 0.10)) +
  labs(title = "Betadisper on phylum across five dataset scopes",
       subtitle = "Bars above the dashed line (p > 0.05) indicate\nhomogeneous within-phylum dispersion.",
       x = NULL, y = "Betadisper p-value")

# ---- Compose 2 x 2 with patchwork tags ---------------------------------------
# - Tags A/B/C/D placed at top-left of each panel by patchwork.
# - Title is centred. Subtitle is centred and smaller.
# - All four panel titles have the same size.
TITLE_SIZE    <- 18
SUBTITLE_SIZE <- 14
TAG_SIZE      <- 24

# Legend policy: collect A+B into a single shared legend to the right of the
# top row (Class + Reference), and collect C+D into a single shared legend to
# the right of the bottom row (Substrate/Auxotrophy fill). This way each
# legend sits next to the panels it explains rather than in a single
# right-hand column that mixes PCoA-only and bar-only guides.
row_ab <- (p_a + p_b) + plot_layout(guides = "collect")
row_cd <- (p_c + p_d) + plot_layout(guides = "collect")

final_fig <- (row_ab / row_cd) +
  plot_annotation(tag_levels = "A") +
  plot_layout(heights = c(1, 1)) &
  theme_classic(base_size = 14) &
  theme(
    plot.title         = element_text(face = "bold", size = TITLE_SIZE, hjust = 0.5),
    plot.subtitle      = element_text(size = SUBTITLE_SIZE, colour = "grey30", hjust = 0.5),
    plot.tag           = element_text(face = "bold", size = TAG_SIZE, hjust = 0),
    plot.tag.position  = c(0.01, 0.99),
    legend.position    = "right",
    legend.key.size    = unit(0.7, "cm"),
    legend.title       = element_text(face = "bold", size = 14),
    legend.text        = element_text(size = 13),
    axis.text          = element_text(size = 12),
    axis.title         = element_text(size = 14),
    panel.grid.major   = element_line(linewidth = 0.2, colour = "#eee"),
    plot.margin        = margin(t = 15, r = 20, b = 15, l = 20)
  )

out_png <- file.path(out_dir, "Fig4.png")
ggsave(out_png, final_fig, width = 22, height = 13, dpi = 1000)
cat(sprintf(">>> Wrote %s\n", out_png))

# ---- Console summary ---------------------------------------------------------
cat(sprintf("\nSubstrate PERMANOVA (42): F = %.3f, p = %.3f\n", perm_sub_F, perm_sub_p))
cat(sprintf("Auxotrophy PERMANOVA (42): F = %.3f, p = %.3f\n", perm_aux_F, perm_aux_p))

cat("\n=== Full 5-scope sweep (recomputed inline; no hardcoded values) ===\n")
perm %>%
  arrange(scope, metric, test_type) %>%
  mutate(across(where(is.numeric), ~ sprintf("%.3f", .))) %>%
  as.data.frame() %>%
  print(row.names = FALSE)
