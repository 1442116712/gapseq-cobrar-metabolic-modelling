# =============================================================
# Figure 1: Cow rumen catalogue overview and SGB selection
#
# Four panels (2 x 2):
#   A. Catalogue class distribution (top 10 classes by n_MAGs)
#   B. MAGs-per-SGB distribution with >=15 threshold highlighted
#   C. Class breakdown of the 34 selected SGBs
#   D. Per-SGB MAG quality (mean completeness vs mean N50)
#
# Colour design (v2):
#   * One palette used throughout: bar/point COLOUR encodes PHYLUM.
#     Same colour means same phylum in every panel.
#   * In panel D, class is encoded by point SHAPE (6 shapes) so the
#     three Firmicutes classes remain distinguishable while sharing
#     the Firmicutes red.
#   * In panels A and C, class names appear on the Y axis;
#     no class-level colour is used to avoid clashing with phylum.
#
# Inputs:
#   M4_Supplement.xlsx
#     - cow-rumen-MAGs   full MGnify cow rumen catalogue v1.0.1
#     - target-SGBs      34 selected SGBs + 2 ARG-MAGs, MAG-level
#
# Output: figures/Fig1_dataset_overview.png
# =============================================================

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
  library(scales)
})

# ---------- Configuration ----------
xlsx_path <- "C:/Users/CFL/OneDrive - Queen's University Belfast/PhD_onedrive/Manuscript4/M4_Supplement.xlsx"
out_dir   <- "C:/Users/CFL/OneDrive - Queen's University Belfast/PhD_onedrive/Manuscript4/figures"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

THRESHOLD <- 15  # MAGs per SGB threshold for Pan-Draft

# ---------- Read data ----------
catalog <- read_excel(xlsx_path, sheet = "cow-rumen-MAGs")
target  <- read_excel(xlsx_path, sheet = "target-SGBs")

# Defensive: empty / NA Phylum -> "Other"
catalog <- catalog %>%
  mutate(Phylum = ifelse(is.na(Phylum) | Phylum == "", "Other", Phylum))

# ---------- Colours and shapes ----------
# ONE phylum palette is used throughout A, C, D.
phylum_colors <- c(
  "Firmicutes"        = "#E76F51",  # red    (Clostridia, Bacilli, Negativicutes)
  "Bacteroidota"      = "#2A9D8F",  # teal   (Bacteroidia)
  "Methanobacteriota" = "#9B5DE5",  # purple (Methanobacteria)
  "Actinobacteriota"  = "#F4A261",  # orange (Coriobacteriia)
  "Cyanobacteria"     = "#7FB069",  # green  (Vampirovibrionia, catalogue only)
  "Verrucomicrobiota" = "#B8B27D",  # olive  (Kiritimatiellia)
  "Proteobacteria"    = "#5C8AE8",  # blue   (Alphaproteobacteria)
  "Spirochaetota"     = "#D4A373",  # tan
  "Other"             = "#9C9C9C"   # grey   (fallback)
)

# Six classes appear in the 34 selected SGBs. Shape encodes class in panel D.
class_shapes <- c(
  "Clostridia"      = 16,  # filled circle
  "Bacilli"         = 17,  # filled triangle
  "Negativicutes"   = 15,  # filled square
  "Bacteroidia"     = 18,  # filled diamond
  "Coriobacteriia"  = 25,  # filled triangle down
  "Methanobacteria" = 8    # asterisk-like
)

# Helper to make Phylum a factor ordered by palette
factor_phylum <- function(x)
  factor(x, levels = intersect(names(phylum_colors), unique(x)))

# ---------- Aggregations ----------
catalog_sgb <- catalog %>%
  group_by(Species_rep, Class, Phylum) %>%
  summarise(n_MAGs = n(), .groups = "drop")

selected_sgbs <- catalog_sgb %>% filter(n_MAGs >= THRESHOLD) %>% pull(Species_rep)
stopifnot(length(selected_sgbs) == 34)

selected_summary <- target %>%
  filter(Species_rep %in% selected_sgbs) %>%
  group_by(Species_rep, Class, Phylum) %>%
  summarise(
    n_MAGs            = n(),
    mean_completeness = mean(Completeness, na.rm = TRUE),
    mean_n50          = mean(n50_contigs,  na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    Phylum = factor_phylum(Phylum),
    Class  = factor(Class, levels = names(class_shapes))
  )
stopifnot(nrow(selected_summary) == 34)

# ---------- Panel A: catalogue class distribution ----------
catalog_class <- catalog %>%
  group_by(Class, Phylum) %>%
  summarise(n_MAGs = n(),
            n_SGBs = n_distinct(Species_rep),
            .groups = "drop") %>%
  arrange(desc(n_MAGs)) %>%
  slice_head(n = 10) %>%
  mutate(Phylum = factor_phylum(Phylum))

p_a <- ggplot(catalog_class,
              aes(x = reorder(Class, n_MAGs), y = n_MAGs, fill = Phylum)) +
  geom_col() +
  geom_text(aes(label = sprintf("%d MAGs (%d SGBs)", n_MAGs, n_SGBs)),
            hjust = -0.05, size = 2.8) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.5))) +
  scale_fill_manual(values = phylum_colors, name = "Phylum",
                    drop = FALSE) +
  coord_flip() +
  labs(title = "A  Catalogue class distribution",
       subtitle = "MGnify cow rumen v1.0.1: 5,578 MAGs in 2,729 SGBs (top 10 classes)",
       x = NULL, y = "Number of MAGs") +
  theme_classic(base_size = 11) +
  theme(legend.position = "right",
        plot.title    = element_text(face = "bold", size = 12),
        plot.subtitle = element_text(size = 9, colour = "grey30"))

# ---------- Panel B: MAGs-per-SGB distribution ----------
breaks <- c(0.5, 1.5, 4.5, 9.5, 14.5, 29.5, 99.5, 1000)
labels <- c("1", "2-4", "5-9", "10-14", "15-29", "30-99", "100+")
bin_data <- catalog_sgb %>%
  mutate(bin = cut(n_MAGs, breaks = breaks, labels = labels, right = TRUE)) %>%
  count(bin) %>%
  mutate(selected = bin %in% c("15-29", "30-99", "100+"))

p_b <- ggplot(bin_data, aes(x = bin, y = n, fill = selected)) +
  geom_col() +
  geom_text(aes(label = n), vjust = -0.4, size = 3.2) +
  geom_vline(xintercept = 4.5, linetype = "dashed", colour = "#264653") +
  scale_fill_manual(values = c("FALSE" = "#cccccc", "TRUE" = "#264653"),
                    labels = c("FALSE" = "Excluded", "TRUE" = "Selected (n = 34)"),
                    name = NULL) +
  scale_y_log10(expand = expansion(mult = c(0, 0.25)),
                labels = comma_format()) +
  labs(title = "B  Distribution of MAG count per SGB",
       subtitle = expression("Pan-Draft accuracy threshold: " >= "15 MAGs per SGB"),
       x = "MAGs per SGB", y = "Number of SGBs (log10)") +
  theme_classic(base_size = 11) +
  theme(legend.position = "top",
        plot.title    = element_text(face = "bold", size = 12),
        plot.subtitle = element_text(size = 9, colour = "grey30"))

# ---------- Panel C: 34 selected SGBs by class (bars by phylum) ----------
selected_class <- selected_summary %>%
  count(Class, Phylum, name = "n_SGBs") %>%
  arrange(desc(n_SGBs))

p_c <- ggplot(selected_class,
              aes(x = reorder(Class, n_SGBs), y = n_SGBs, fill = Phylum)) +
  geom_col() +
  geom_text(aes(label = sprintf("n = %d", n_SGBs)),
            hjust = -0.15, size = 3.2, fontface = "bold") +
  scale_fill_manual(values = phylum_colors, name = "Phylum",
                    drop = FALSE) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.20)),
                     breaks = c(0, 5, 10, 15)) +
  coord_flip() +
  labs(title = "C  Selected SGBs by class",
       subtitle = "34 SGBs span six classes across four phyla",
       x = NULL, y = "Number of SGBs") +
  theme_classic(base_size = 11) +
  theme(legend.position = "right",
        plot.title    = element_text(face = "bold", size = 12),
        plot.subtitle = element_text(size = 9, colour = "grey30"))

# ---------- Panel D: per-SGB MAG quality ----------
# Color = Phylum (consistent with A and C). Shape = Class (six classes).
p_d <- ggplot(selected_summary,
              aes(x = mean_completeness, y = mean_n50 / 1000,
                  colour = Phylum, shape = Class, size = n_MAGs)) +
  geom_point(alpha = 0.85, stroke = 0.4) +
  scale_colour_manual(values = phylum_colors, name = "Phylum",
                      drop = FALSE) +
  scale_shape_manual(values = class_shapes, name = "Class") +
  scale_size_continuous(range = c(2.5, 8), breaks = c(15, 30, 60, 84),
                        name = "MAGs per SGB") +
  guides(
    colour = guide_legend(order = 1, override.aes = list(shape = 16, size = 3.5)),
    shape  = guide_legend(order = 2, override.aes = list(colour = "grey30", size = 3.5)),
    size   = guide_legend(order = 3)
  ) +
  labs(title = "D  Per-SGB MAG quality",
       subtitle = "Each point = one SGB. Colour = phylum, shape = class.",
       x = "Mean completeness (%)",
       y = "Mean contig N50 (kb)") +
  theme_classic(base_size = 11) +
  theme(legend.position = "right",
        legend.box      = "vertical",
        plot.title      = element_text(face = "bold", size = 12),
        plot.subtitle   = element_text(size = 9, colour = "grey30"))

# ---------- Combine ----------
fig1 <- (p_a | p_b) / (p_c | p_d) +
  plot_layout(heights = c(1, 1))

out_png <- file.path(out_dir, "Fig1_dataset_overview.png")
ggsave(out_png, fig1, width = 14, height = 10, dpi = 300)
cat(sprintf(">>> Wrote %s\n", out_png))

# ---------- Console diagnostics ----------
cat("\n=== Summary of selected dataset ===\n")
cat(sprintf("  Total catalogue: %d MAGs in %d SGBs\n",
            nrow(catalog), n_distinct(catalog$Species_rep)))
cat(sprintf("  Selected SGBs (>=15 MAGs): %d\n", length(selected_sgbs)))
cat(sprintf("  MAGs in selected SGBs: %d\n", sum(selected_summary$n_MAGs)))
cat("  Class breakdown:\n")
print(selected_class)
cat("\n  Quality:\n")
cat(sprintf("    Mean completeness across SGBs: median = %.2f%% (range %.2f-%.2f%%)\n",
            median(selected_summary$mean_completeness),
            min(selected_summary$mean_completeness),
            max(selected_summary$mean_completeness)))
cat(sprintf("    Mean N50          across SGBs: median = %.0f bp (range %.0f-%.0f bp)\n",
            median(selected_summary$mean_n50),
            min(selected_summary$mean_n50),
            max(selected_summary$mean_n50)))
