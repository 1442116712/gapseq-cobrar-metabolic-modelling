# =============================================================
# Fig_substrate_assay.R
#
# Full six-panel Figure 6 for the v4.1 SH-panel experiments, drawn from
# S12e_Validation_result:
#
#   A. SH_0542 on M2 (5 conditions)     B. SH_0542 on FAB (6-condition factorial)
#   C. SH_0717 on M2 (5 conditions)     D. SH_0760 on M2 (5 conditions)
#   E. SH_0717 FAB rescue               F. SH_0760 FAB rescue
#      (FAB alone / FAB + strain-specific cocktail / M2 positive control; 0,24,48,72 h)
#
# Panels A-D use the upper (substrate-assay) block of S12e; panels E-F use the
# lower (FAB-rescue) block. Substrate rows are mapped to protocol conditions by
# block-row index rather than by string label.
# =============================================================

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
})

# ---------- Paths ----------
xlsx_path <- "C:/Users/CFL/OneDrive - Queen's University Belfast/PhD_onedrive/Manuscript4_gapseq/M4_Supplementary_Tables.xlsx"
out_dir   <- "C:/Users/CFL/OneDrive - Queen's University Belfast/PhD_onedrive/Manuscript4_gapseq/figures"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ---------- Colour palette (indexed by internal short condition code) ----------
cond_colours <- c(
  "Baseline"           = "#95A5A6",  # grey  - relabelled per panel to "M2" or "FAB"
  "+ maltose"          = "#5B9BD5",  # blue
  "+ strain-specific"  = "#F2C94C",  # yellow
  "+ combined"         = "#C9622B",  # orange (M2 public + specific)
  "+ thiamine"         = "#9B59B6",  # purple (FAB only)
  "+ urea (negative)"  = "#2E8B57",  # green - distinct from grey baseline
  "+ all three"        = "#B22222"   # red (FAB only)
)

cond_levels_all <- names(cond_colours)

# ---------- Read S12e once and split into substrate (A-D) / rescue (E-F) blocks ----------
# The rescue block was appended below the substrate block and has its own header
# row ("OD600/time(h) | 0h | 24h | 48h | 72h"). Reading the whole sheet in one go
# mixes the two blocks and makes some time columns character, so we split first
# on the header rows and coerce each block's time columns to numeric.
time_pts <- c(0:15, 22, 23, 24, 48, 72)
allc <- read_excel(xlsx_path, sheet = "S12e_Validation_result", col_names = FALSE)
hdr_rows <- which(allc[[1]] == "OD600/time(h)")
main_hdr <- hdr_rows[1]
resc_hdr <- hdr_rows[length(hdr_rows)]

raw <- allc[(main_hdr + 1):(resc_hdr - 1), 1:(2 + length(time_pts))]
colnames(raw) <- c("strain_tag", "condition_raw", paste0("t", time_pts))
raw <- raw %>%
  mutate(across(starts_with("t"), ~ suppressWarnings(as.numeric(.)))) %>%
  mutate(strain_tag = ifelse(is.na(strain_tag) | strain_tag == "", NA, strain_tag)) %>%
  tidyr::fill(strain_tag, .direction = "down") %>%
  filter(strain_tag %in% c("SH542", "SH717", "SH760"))

sh542_rows <- which(raw$strain_tag == "SH542")
b_first <- sh542_rows[which(startsWith(raw$condition_raw[sh542_rows], "FAB"))[1]]
raw$experiment <- ifelse(seq_len(nrow(raw)) >= b_first & raw$strain_tag == "SH542", "B", "A")
raw$block_row <- ave(seq_len(nrow(raw)), raw$strain_tag, raw$experiment, FUN = seq_along)

assign_condition_A <- function(i) {
  if (i <= 3)  return("Baseline")
  if (i <= 6)  return("+ maltose")
  if (i <= 9)  return("+ strain-specific")
  if (i <= 12) return("+ combined")
  if (i <= 15) return("+ urea (negative)")
  NA_character_
}
assign_condition_B <- function(i) {
  if (i <= 3)  return("Baseline")
  if (i <= 6)  return("+ maltose")
  if (i <= 9)  return("+ strain-specific")
  if (i <= 12) return("+ thiamine")
  if (i <= 15) return("+ all three")
  if (i <= 18) return("+ urea (negative)")
  NA_character_
}

raw$condition <- mapply(function(exp, i) {
  if (exp == "A") assign_condition_A(i) else assign_condition_B(i)
}, raw$experiment, raw$block_row)

data <- raw %>%
  filter(strain_tag %in% c("SH542", "SH717", "SH760"), !is.na(condition))

long <- data %>%
  select(strain_tag, experiment, condition, starts_with("t")) %>%
  pivot_longer(cols = starts_with("t"), names_to = "timepoint", values_to = "od") %>%
  mutate(time_h = as.numeric(sub("^t", "", timepoint)), od = as.numeric(od)) %>%
  filter(!is.na(od))

summ <- long %>%
  group_by(strain_tag, experiment, condition, time_h) %>%
  summarise(mean_od = mean(od), sd_od = sd(od), n = n(), .groups = "drop") %>%
  mutate(ymin = pmax(mean_od - sd_od, 0),
         ymax = mean_od + sd_od)

# ---------- FAB-rescue block (E-F): rows below the rescue header ----------
resc <- allc[(resc_hdr + 1):nrow(allc), 1:6]
colnames(resc) <- c("strain_tag", "condition", "t0", "t24", "t48", "t72")
resc <- resc %>%
  mutate(across(starts_with("t"), ~ suppressWarnings(as.numeric(.)))) %>%
  mutate(strain_tag = ifelse(is.na(strain_tag) | strain_tag == "", NA, strain_tag)) %>%
  tidyr::fill(strain_tag, .direction = "down") %>%
  filter(!is.na(condition) & condition != "")

resc_long <- resc %>%
  pivot_longer(cols = starts_with("t"), names_to = "timepoint", values_to = "od") %>%
  mutate(time_h = as.numeric(sub("^t", "", timepoint)), od = as.numeric(od)) %>%
  filter(!is.na(od))

resc_summ <- resc_long %>%
  group_by(strain_tag, condition, time_h) %>%
  summarise(mean_od = mean(od), sd_od = sd(od), n = n(), .groups = "drop") %>%
  mutate(ymin = pmax(mean_od - sd_od, 0),
         ymax = mean_od + sd_od)

# ---------- Helper: growth-curve panel (substrate assays A-D) ----------
plot_curves <- function(df, label_map, title = "", subtitle = "",
                        ymax_val = NULL, legend_ncol = 3) {
  df$condition <- factor(df$condition, levels = names(label_map))
  df$condition <- droplevels(df$condition)
  present <- levels(df$condition)
  present_labels <- label_map[present]

  p <- ggplot(df, aes(x = time_h, y = mean_od,
                      colour = condition, fill = condition,
                      group = condition)) +
    geom_errorbar(aes(ymin = ymin, ymax = ymax), width = 0.9, linewidth = 0.45, alpha = 0.95) +
    geom_line(linewidth = 1.1) +
    geom_point(size = 1.9, shape = 21, stroke = 0.4, colour = "grey20") +
    scale_colour_manual(values = cond_colours[present], labels = present_labels, name = NULL) +
    scale_fill_manual(values   = cond_colours[present], labels = present_labels, name = NULL) +
    scale_x_continuous(breaks = c(0:15, 24, 48, 72), expand = expansion(mult = c(0.01, 0.02))) +
    guides(colour = guide_legend(ncol = legend_ncol, byrow = FALSE,
                                 label.position = "right",
                                 label.hjust = 0),
           fill   = guide_legend(ncol = legend_ncol, byrow = FALSE,
                                 label.position = "right",
                                 label.hjust = 0)) +
    labs(title = title, subtitle = subtitle,
         x = "Time (h)",
         y = expression(paste("Blank-corrected OD"[600]))) +
    theme_classic(base_size = 14) +
    theme(legend.position = "bottom",
          legend.text = element_text(size = 12, lineheight = 1.0, hjust = 0),
          legend.key.width = unit(1.1, "cm"),
          legend.key.height = unit(0.9, "cm"),
          legend.spacing.x = unit(0.4, "cm"),
          legend.spacing.y = unit(0.15, "cm"),
          legend.box.margin = margin(t = 4, r = 4, b = 4, l = 4),
          axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
          plot.title = element_text(face = "bold", size = 15),
          plot.subtitle = element_text(size = 11, colour = "grey30"),
          plot.margin = margin(10, 12, 10, 12))
  if (!is.null(ymax_val)) p <- p + coord_cartesian(ylim = c(0, ymax_val))
  p
}

# ---------- Helper: FAB-rescue panel (E-F), same visual grammar ----------
resc_colours <- c("FAB" = "#95A5A6", "FAB+cocktail" = "#B22222", "M2" = "#2E8B57")
plot_rescue <- function(df, title = "", subtitle = "3 conditions \u00d7 n = 3",
                        cocktail_label = "+ cocktail") {
  lab_map <- c("FAB" = "FAB alone",
               "FAB+cocktail" = cocktail_label,
               "M2" = "M2 (positive control)")
  df$condition <- factor(df$condition, levels = names(lab_map))
  df$condition <- droplevels(df$condition)
  present <- levels(df$condition)

  ggplot(df, aes(x = time_h, y = mean_od,
                 colour = condition, fill = condition, group = condition)) +
    geom_errorbar(aes(ymin = ymin, ymax = ymax), width = 1.6, linewidth = 0.5, alpha = 0.95) +
    geom_line(linewidth = 1.1) +
    geom_point(size = 2.1, shape = 21, stroke = 0.4, colour = "grey20") +
    scale_colour_manual(values = resc_colours[present], labels = lab_map[present], name = NULL) +
    scale_fill_manual(values   = resc_colours[present], labels = lab_map[present], name = NULL) +
    scale_x_continuous(breaks = c(0, 24, 48, 72), limits = c(0, 72),
                       expand = expansion(mult = c(0, 0.02))) +
    coord_cartesian(ylim = c(0, 1.0)) +
    guides(colour = guide_legend(ncol = 3), fill = guide_legend(ncol = 3)) +
    labs(title = title, subtitle = subtitle,
         x = "Time (h)",
         y = expression(paste("Blank-corrected OD"[600]))) +
    theme_classic(base_size = 14) +
    theme(legend.position = "bottom",
          legend.text = element_text(size = 12, lineheight = 1.0, hjust = 0),
          legend.key.width = unit(1.1, "cm"),
          legend.key.height = unit(0.9, "cm"),
          legend.spacing.x = unit(0.4, "cm"),
          legend.box.margin = margin(t = 4, r = 4, b = 4, l = 4),
          plot.title = element_text(face = "bold", size = 15),
          plot.subtitle = element_text(size = 11, colour = "grey30"),
          plot.margin = margin(10, 12, 10, 12))
}

# ===============================================
# COMPOSITE FIGURE: 3x2 grid
#   Row 1: A. SH_0542 M2      | B. SH_0542 FAB
#   Row 2: C. SH_0717 M2      | D. SH_0760 M2
#   Row 3: E. SH_0717 rescue  | F. SH_0760 rescue
# ===============================================
sh542_m2  <- summ %>% filter(strain_tag == "SH542", experiment == "A")
sh542_fab <- summ %>% filter(strain_tag == "SH542", experiment == "B")
sh717     <- summ %>% filter(strain_tag == "SH717", experiment == "A")
sh760     <- summ %>% filter(strain_tag == "SH760", experiment == "A")

# Common y-max for the two SH_0542 panels so M2 vs FAB is visually comparable
ymax_542 <- max(c(sh542_m2$ymax, sh542_fab$ymax), na.rm = TRUE) * 1.05

# Panel A: SH_0542 on M2
labels_A <- c(
  "Baseline"          = "M2",
  "+ maltose"         = "+ 0.1% (w/v) maltose",
  "+ strain-specific" = "+ 0.1% (w/v) galactose",
  "+ combined"        = "+ 0.1% (w/v) maltose\n+ 0.1% (w/v) galactose",
  "+ urea (negative)" = "+ 0.1% (w/v) urea"
)
p_A <- plot_curves(sh542_m2, label_map = labels_A,
  title = "A. SH_0542 (Coprobacillaceae) on M2",
  subtitle = "5 conditions \u00d7 n = 3", ymax_val = ymax_542)

# Panel B: SH_0542 on FAB
labels_B <- c(
  "Baseline"          = "FAB",
  "+ maltose"         = "+ 0.1% (w/v) maltose",
  "+ strain-specific" = "+ 0.1% (w/v) galactose",
  "+ thiamine"        = "+ 5 \u00b5M thiamine\u00b7HCl",
  "+ urea (negative)" = "+ 0.1% (w/v) urea",
  "+ all three"       = "+ 0.1% (w/v) maltose\n+ 0.1% (w/v) galactose\n+ 5 \u00b5M thiamine\u00b7HCl"
)
p_B <- plot_curves(sh542_fab, label_map = labels_B,
  title = "B. SH_0542 (Coprobacillaceae) on FAB",
  subtitle = "6-condition factorial \u00d7 n = 3", ymax_val = ymax_542)

# Panel C: SH_0717 on M2
labels_C <- c(
  "Baseline"          = "M2",
  "+ maltose"         = "+ 0.1% (w/v) maltose",
  "+ strain-specific" = "+ 0.1% (w/v) xylose",
  "+ combined"        = "+ 0.1% (w/v) maltose\n+ 0.1% (w/v) xylose",
  "+ urea (negative)" = "+ 0.1% (w/v) urea"
)
p_C <- plot_curves(sh717, label_map = labels_C,
  title = "C. SH_0717 (Ruminococcaceae) on M2",
  subtitle = "5 conditions \u00d7 n = 3")

# Panel D: SH_0760 on M2
labels_D <- c(
  "Baseline"          = "M2",
  "+ maltose"         = "+ 0.1% (w/v) maltose",
  "+ strain-specific" = "+ 0.1% (w/v) cellobiose",
  "+ combined"        = "+ 0.1% (w/v) maltose\n+ 0.1% (w/v) cellobiose",
  "+ urea (negative)" = "+ 0.1% (w/v) urea"
)
p_D <- plot_curves(sh760, label_map = labels_D,
  title = "D. SH_0760 (Erysipelotrichaceae) on M2",
  subtitle = "5 conditions \u00d7 n = 3")

# Panel E: SH_0717 FAB rescue
p_E <- plot_rescue(resc_summ %>% filter(strain_tag == "SH717"),
  title = "E. SH_0717 (Ruminococcaceae) FAB rescue",
  cocktail_label = "+ cocktail (metals, vitamins, chorismate)")

# Panel F: SH_0760 FAB rescue
p_F <- plot_rescue(resc_summ %>% filter(strain_tag == "SH760"),
  title = "F. SH_0760 (Erysipelotrichaceae) FAB rescue",
  cocktail_label = "+ cocktail (metals, vitamins, D-glucose)")

# ---------- Compose 3x2 grid ----------
fig <- (p_A | p_B) / (p_C | p_D) / (p_E | p_F) + plot_layout(heights = c(1, 1, 1))

# Full-page figure: wide panels so the hourly x-axis, labels and legend are legible.
ggsave(file.path(out_dir, "Fig_substrate_assay.png"), fig,
       width = 18, height = 22, units = "in", dpi = 400, bg = "white")
cat("Saved: Fig_substrate_assay.png (panels A-F, full page)\n")

# ---------- Also save summary CSVs ----------
write.csv(summ,      file.path(out_dir, "Fig_substrate_assay_summary.csv"), row.names = FALSE)
write.csv(resc_summ, file.path(out_dir, "Fig_substrate_rescue_summary.csv"), row.names = FALSE)
cat("Summary CSVs saved.\n")
