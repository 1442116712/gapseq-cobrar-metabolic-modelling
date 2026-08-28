# =============================================================================
# Fig_trace_sensitivity.R
# Manuscript 4 - Trace background growth sensitivity to essential-factor
# uptake bound. 42 organisms tested across 4 lower bounds (-0.01, -0.1, -1, -10
# mmol/gDW/h). Calibrates the choice of -0.1 used in the main pipeline.
#
# Panel A: 42 trajectories on log-log axes; QC bands at <0.01, 0.01-0.1, >0.1.
#          Reference (gapseq) and reference (rumen) annotated; ARG-MAGs marked
#          as red diamonds. Insensitive (max/min < 1.5x) shown solid;
#          sensitive shown dashed.
# Panel B: Sensitive subset only (n = ~18), horizontal bars at 4 thresholds,
#          ordered ARG-MAG -> reference -> pan, descending max trace_bg.
# Right column: role legend, sensitivity legend, and group-composition box.
# =============================================================================

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(patchwork)
  library(ggrepel)
})

# ---- Inputs ------------------------------------------------------------------
xl      <- "C:/Users/CFL/OneDrive - Queen's University Belfast/PhD_onedrive/Manuscript4/M4_Supplement.xlsx"
out_dir <- "C:/Users/CFL/OneDrive - Queen's University Belfast/PhD_onedrive/Manuscript4/figures"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

INSENS_FOLD <- 1.5   # max/min ratio threshold for "insensitive" classification

# ---- Read sensitivity sheet --------------------------------------------------
ts <- read_excel(xl, sheet = "trace_sensitivity")

# Expected columns: SGB, trace_0p01, trace_0p1, trace_1, trace_10
stopifnot(all(c("SGB","trace_0p01","trace_0p1","trace_1","trace_10")
              %in% colnames(ts)))

# Force numeric (some cells may be read as character if Excel contains
# "NA" strings, blanks, or stray non-numeric tokens)
ts <- ts %>%
  mutate(across(c(trace_0p01, trace_0p1, trace_1, trace_10),
                ~ suppressWarnings(as.numeric(.))))

# ---- Classify and annotate ---------------------------------------------------
classify_sensitivity <- function(v01, v1, v_1, v_10, fold = INSENS_FOLD) {
  vals <- c(v01, v1, v_1, v_10)
  vals <- vals[!is.na(vals)]
  if (length(vals) < 2) return("unknown")
  if (max(vals) / min(vals) < fold) "insensitive" else "sensitive"
}

ts <- ts %>%
  rowwise() %>%
  mutate(group = classify_sensitivity(trace_0p01, trace_0p1, trace_1, trace_10)) %>%
  ungroup()

clinical_set <- c("Ecoli_K12_MG1655","Btheta_VPI5482","Efaecalis_ATCC19433")
rumen_set    <- c("Ecoli_PA3","Btheta_KPPR3","Efaecalis_68A")
argmag_set   <- c("MGYG000292883","MGYG000295553")

ts <- ts %>%
  mutate(role = case_when(
    SGB %in% clinical_set ~ "ref_gapseq",
    SGB %in% rumen_set    ~ "ref_rumen",
    SGB %in% argmag_set   ~ "argmag",
    TRUE                  ~ "pan"))

display_label <- function(s) {
  case_when(
    s == "Btheta_KPPR3"        ~ "B. theta KPPR-3",
    s == "Btheta_VPI5482"      ~ "B. theta VPI-5482",
    s == "Ecoli_K12_MG1655"    ~ "E. coli K-12",
    s == "Ecoli_PA3"           ~ "E. coli PA-3",
    s == "Efaecalis_68A"       ~ "E. faecalis 68A",
    s == "Efaecalis_ATCC19433" ~ "E. faecalis ATCC19433",
    s == "MGYG000292883"       ~ "MAG-Lent (M292883)",
    s == "MGYG000295553"       ~ "MAG-Kirit (M295553)",
    TRUE ~ sub("_pan$", "", s))
}
ts <- ts %>% mutate(display = display_label(SGB))

# ---- Long form for trajectories ---------------------------------------------
ts_long <- ts %>%
  pivot_longer(cols = c(trace_0p01, trace_0p1, trace_1, trace_10),
               names_to = "threshold", values_to = "trace_bg") %>%
  mutate(threshold_value = case_when(
    threshold == "trace_0p01" ~ 0.01,
    threshold == "trace_0p1"  ~ 0.1,
    threshold == "trace_1"    ~ 1,
    threshold == "trace_10"   ~ 10),
    threshold_label = case_when(
      threshold == "trace_0p01" ~ "-0.01",
      threshold == "trace_0p1"  ~ "-0.1",
      threshold == "trace_1"    ~ "-1",
      threshold == "trace_10"   ~ "-10")) %>%
  filter(!is.na(trace_bg))

# ---- Palettes ----------------------------------------------------------------
ROLE_PAN <- "#2E86AB"
ROLE_RG  <- "#D4A017"
ROLE_RR  <- "#2E7D32"
ROLE_ARG <- "#B22222"

role_pal <- c(pan = ROLE_PAN, ref_gapseq = ROLE_RG,
              ref_rumen = ROLE_RR, argmag = ROLE_ARG)
group_lty <- c(insensitive = "solid", sensitive = "22")
threshold_pal <- c(`-0.01` = "#A8C8E1", `-0.1` = "#5BA0D0",
                   `-1` = "#2E5F8A",  `-10` = "#0F2A4F")

# =============================================================================
# Panel A - Trajectories
# =============================================================================
xs <- c(0.01, 0.1, 1, 10)

# Bring ARG-MAGs and references to top z-order by drawing them last
draw_order <- c("pan", "ref_gapseq", "ref_rumen", "argmag")
ts_long_ord <- ts_long %>%
  mutate(role = factor(role, levels = draw_order)) %>%
  arrange(role, SGB)

# Build label coordinates: each annotated organism uses its trace_10 value
label_df <- ts %>%
  filter(role != "pan") %>%
  mutate(x_end = 10, y_end = trace_10) %>%
  filter(!is.na(y_end))

p_a <- ggplot() +
  # QC bands
  annotate("rect", xmin = 0.005, xmax = 60, ymin = 0,    ymax = 0.01,
           fill = "#E8F4D9", alpha = 0.5) +
  annotate("rect", xmin = 0.005, xmax = 60, ymin = 0.01, ymax = 0.1,
           fill = "#FFF3D6", alpha = 0.5) +
  annotate("rect", xmin = 0.005, xmax = 60, ymin = 0.1,  ymax = 5,
           fill = "#FCE4E4", alpha = 0.5) +
  annotate("text", x = 0.0075, y = 0.003, label = "pass\n<0.01",
           fontface = "bold", size = 2.6, colour = "#5D8C2E") +
  annotate("text", x = 0.0075, y = 0.03,  label = "borderline\n0.01-0.1",
           fontface = "bold", size = 2.6, colour = "#B8860B") +
  annotate("text", x = 0.0075, y = 0.5,   label = "fail\n>0.1",
           fontface = "bold", size = 2.6, colour = "#A53A3A") +
  # Trajectories - aes(colour = role) drives legend; aes(linetype = group) drives sensitivity legend
  geom_line(data = ts_long_ord %>% filter(role == "pan"),
            aes(threshold_value, trace_bg, group = SGB,
                colour = role, linetype = group),
            linewidth = 0.4, alpha = 0.55) +
  geom_point(data = ts_long_ord %>% filter(role == "pan"),
             aes(threshold_value, trace_bg, colour = role),
             size = 1.0, alpha = 0.55, show.legend = FALSE) +
  geom_line(data = ts_long_ord %>% filter(role == "ref_gapseq"),
            aes(threshold_value, trace_bg, group = SGB,
                colour = role, linetype = group),
            linewidth = 0.9, alpha = 0.95) +
  geom_point(data = ts_long_ord %>% filter(role == "ref_gapseq"),
             aes(threshold_value, trace_bg, colour = role),
             size = 1.8, alpha = 0.95, show.legend = FALSE) +
  geom_line(data = ts_long_ord %>% filter(role == "ref_rumen"),
            aes(threshold_value, trace_bg, group = SGB,
                colour = role, linetype = group),
            linewidth = 0.9, alpha = 0.95) +
  geom_point(data = ts_long_ord %>% filter(role == "ref_rumen"),
             aes(threshold_value, trace_bg, colour = role),
             size = 1.8, alpha = 0.95, show.legend = FALSE) +
  geom_line(data = ts_long_ord %>% filter(role == "argmag"),
            aes(threshold_value, trace_bg, group = SGB,
                colour = role, linetype = group),
            linewidth = 1.0, alpha = 1.0) +
  geom_point(data = ts_long_ord %>% filter(role == "argmag"),
             aes(threshold_value, trace_bg, colour = role),
             size = 2.4, shape = 18, show.legend = FALSE) +
  # Annotation labels
  geom_text_repel(data = label_df,
                  aes(x_end, y_end, label = display, colour = role),
                  size = 2.6, fontface = "italic",
                  nudge_x = 4, direction = "y", segment.size = 0.25,
                  segment.colour = "grey60", min.segment.length = 0,
                  max.overlaps = Inf, seed = 7,
                  show.legend = FALSE) +
  scale_colour_manual(values = role_pal,
                      breaks = c("pan","ref_gapseq","ref_rumen","argmag"),
                      labels = c("pan-Draft (n=34)",
                                 "Reference, gapseq (n=3)",
                                 "Reference, rumen (n=3)",
                                 "ARG-MAG (n=2)"),
                      name = "Role") +
  scale_linetype_manual(values = group_lty,
                        breaks = c("insensitive","sensitive"),
                        labels = c("Insensitive (max/min < 1.5x)",
                                   "Sensitive"),
                        name = "Sensitivity") +
  scale_x_log10(breaks = xs, labels = c("-0.01","-0.1","-1","-10"),
                limits = c(0.005, 60), expand = c(0, 0)) +
  scale_y_log10(limits = c(5e-4, 5)) +
  labs(title = "A   Trace background sensitivity (n = 42)",
       x = expression(paste("Essential-factor exchange lower bound (mmol gDW"^-1, " h"^-1, ")")),
       y = expression(paste("Trace background growth (h"^-1, ")"))) +
  guides(colour   = guide_legend(override.aes = list(linewidth = 1.2,
                                                     alpha = 1)),
         linetype = guide_legend(override.aes = list(colour = "#666",
                                                     linewidth = 0.9))) +
  theme_classic(base_size = 10) +
  theme(plot.title       = element_text(size = 13, hjust = 0.5,
                                        face = "bold"),
        legend.position  = "right",
        legend.key.size  = unit(0.5, "cm"),
        legend.key.width = unit(1.2, "cm"),
        legend.spacing.y = unit(0.2, "cm"),
        panel.grid.major = element_line(linewidth = 0.2, colour = "#eee"))

# =============================================================================
# Panel B - Sensitive subset bar chart
# =============================================================================
sens <- ts %>% filter(group == "sensitive") %>%
  mutate(role_order = case_when(role == "argmag" ~ 0,
                                role %in% c("ref_gapseq","ref_rumen") ~ 1,
                                TRUE ~ 2),
         max_trace = pmax(trace_0p01, trace_0p1, trace_1, trace_10, na.rm = TRUE)) %>%
  arrange(role_order, desc(max_trace))

sens_long <- sens %>%
  mutate(display = factor(display, levels = rev(display))) %>%
  pivot_longer(cols = c(trace_0p01, trace_0p1, trace_1, trace_10),
               names_to = "threshold", values_to = "trace_bg") %>%
  mutate(threshold_label = case_when(
    threshold == "trace_0p01" ~ "-0.01",
    threshold == "trace_0p1"  ~ "-0.1",
    threshold == "trace_1"    ~ "-1",
    threshold == "trace_10"   ~ "-10"),
    threshold_label = factor(threshold_label,
                             levels = c("-0.01","-0.1","-1","-10")))

# Y-axis label colour
ylab_colour <- sens %>%
  mutate(col = case_when(role == "argmag"     ~ ROLE_ARG,
                         role == "ref_gapseq" ~ ROLE_RG,
                         role == "ref_rumen"  ~ ROLE_RR,
                         TRUE                 ~ "#222")) %>%
  pull(col)

p_b <- ggplot(sens_long,
              aes(x = trace_bg, y = display, fill = threshold_label)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.65,
           colour = "white", linewidth = 0.2) +
  geom_vline(xintercept = 0.01, linetype = "dashed",
             colour = "#5D8C2E", linewidth = 0.4) +
  annotate("text", x = 0.01, y = nrow(sens) + 1.6,
           label = "QC threshold (0.01)", colour = "#5D8C2E",
           hjust = 0.5, vjust = 0, size = 2.7, fontface = "bold") +
  scale_fill_manual(values = threshold_pal,
                    name = expression(paste("ess.-factor lb (mmol gDW"^-1, " h"^-1, ")"))) +
  scale_x_continuous(trans = scales::pseudo_log_trans(sigma = 0.001),
                     breaks = c(0, 0.001, 0.01, 0.1, 1),
                     labels = c("0","0.001","0.01","0.1","1"),
                     limits = c(0, 3)) +
  coord_cartesian(ylim = c(0.5, nrow(sens) + 1.0), clip = "off") +
  labs(title = sprintf("B   Sensitive subset (n = %d)", nrow(sens)),
       x = expression(paste("Trace background growth (h"^-1, ")")),
       y = NULL) +
  theme_classic(base_size = 10) +
  theme(plot.title       = element_text(size = 13, hjust = 0.5,
                                        face = "bold",
                                        margin = margin(b = 14)),
        axis.text.y      = element_text(colour = rev(ylab_colour),
                                        face = "bold"),
        legend.position  = "right",
        legend.key.size  = unit(0.4, "cm"),
        plot.margin      = margin(t = 14, r = 8, b = 8, l = 8),
        panel.grid.major.x = element_line(linewidth = 0.2, colour = "#eee"))

# =============================================================================
# Counts text for caption (printed to console for use in figure caption)
# =============================================================================
n_pan_ins <- sum(ts$role == "pan"        & ts$group == "insensitive")
n_pan_sen <- sum(ts$role == "pan"        & ts$group == "sensitive")
n_ref_ins <- sum(ts$role %in% c("ref_gapseq","ref_rumen") & ts$group == "insensitive")
n_ref_sen <- sum(ts$role %in% c("ref_gapseq","ref_rumen") & ts$group == "sensitive")
n_arg_ins <- sum(ts$role == "argmag"     & ts$group == "insensitive")
n_arg_sen <- sum(ts$role == "argmag"     & ts$group == "sensitive")

cat("\n=== Group composition ===\n")
cat(sprintf("Insensitive (n=%d):  pan %d/34, ref %d/6, ARG-MAG %d/2\n",
            n_pan_ins + n_ref_ins + n_arg_ins,
            n_pan_ins, n_ref_ins, n_arg_ins))
cat(sprintf("Sensitive   (n=%d):  pan %d/34, ref %d/6, ARG-MAG %d/2\n",
            n_pan_sen + n_ref_sen + n_arg_sen,
            n_pan_sen, n_ref_sen, n_arg_sen))

# =============================================================================
# Final figure
# =============================================================================
fig <- p_a / p_b + plot_layout(heights = c(1.0, 1.4))

ggsave(file.path(out_dir, "Supplementary Figure S4.png"),
       fig, width = 12, height = 11, dpi = 1000)

cat("\nFigure saved.\n")

# Persist the classified table
write.table(ts, file = file.path(out_dir, "trace_sensitivity_classified.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)
