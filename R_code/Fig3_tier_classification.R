# =============================================================
# Fig 3 — Three-axis tier classification of inferred essential
# factors across the 42-model dataset.
#
# Inputs (single Excel file with 4 sheets):
#   S10a_Auxotrophy_tier     -- per-(model x compound) tier
#   S2a_Per_model_summary    -- model_id, model_type, category
#   S1b_Selected_SGBs        -- Species_rep -> Class mapping
#   (read sheets with skip = 3 because of title/desc/blank rows)
#
# Output: figures/Fig3_tier_classification.png
#
# Tiers (from S10a column 'tier'):
#   A_hard              -> "Tier A" (robust auxotroph)
#   B_gap_fill_rescued  -> "Tier B" (gap-fill bridged synthesis)
#   C_network_coupled   -> "Tier C" (condition-dependent)
#   D_anomaly           -> "Tier D" (rare edge case)
# =============================================================

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(forcats)
  library(ggnewscale)   # second fill scale (class colour strip)
})

# ---------- Paths (edit if needed) ----------
xlsx_path <- "C:/Users/CFL/OneDrive - Queen's University Belfast/PhD_onedrive/Manuscript4/M4_Supplementary_Tables.xlsx"
out_path  <- "C:/Users/CFL/OneDrive - Queen's University Belfast/PhD_onedrive/Manuscript4/figures/Fig3_tier_classification.png"

# ---------- Load sheets ----------
tier <- read_excel(xlsx_path, sheet = "S10a_Auxotrophy_tier", skip = 3) |>
  filter(in_main_dataset == TRUE)
s2a  <- read_excel(xlsx_path, sheet = "S2a_Per_model_summary", skip = 3)
s1b  <- read_excel(xlsx_path, sheet = "S1b_Selected_SGBs",     skip = 3)

# ---------- Resolve metadata per organism ----------
meta_lookup <- function(o) {
  hit <- match(o, s2a$model_id)
  if (!is.na(hit)) return(list(type = s2a$model_type[hit], cat = s2a$category[hit]))
  if (endsWith(o, "_pan")) {
    hit <- match(sub("_pan$", "", o), s2a$model_id)
    if (!is.na(hit)) return(list(type = s2a$model_type[hit], cat = s2a$category[hit]))
  }
  list(type = NA_character_, cat = NA_character_)
}
class_lookup <- function(o) {
  if (o == "MGYG000292883") return("Lentisphaeria")
  if (o == "MGYG000295553") return("Kiritimatiellia")
  if (endsWith(o, "_pan")) {
    sgb <- sub("_pan$", "", o)
    hit <- match(sgb, s1b$Species_rep)
    if (!is.na(hit)) return(s1b$Class[hit])
  }
  NA_character_
}

# ---------- Per-organism tier counts ----------
per_org <- tier |>
  count(organism, tier) |>
  pivot_wider(names_from = tier, values_from = n, values_fill = 0)

# Ensure all four tier columns exist (some organisms may have 0 of some)
for (col in c("A_hard","B_gap_fill_rescued","C_network_coupled","D_anomaly")) {
  if (!col %in% names(per_org)) per_org[[col]] <- 0
}

per_org <- per_org |>
  rowwise() |>
  mutate(category = meta_lookup(organism)$cat,
         class    = class_lookup(organism)) |>
  ungroup() |>
  mutate(group = case_when(
    grepl("Reference", category) ~ "Reference",
    grepl("ARG-MAG",   category) ~ "ARG-MAG",
    TRUE                         ~ "Pan-Draft"
  )) |>
  mutate(group = factor(group, levels = c("Reference","ARG-MAG","Pan-Draft"))) |>
  arrange(group, class, organism) |>
  mutate(organism = fct_inorder(organism),
         n_essential = A_hard + C_network_coupled + D_anomaly,
         n_total     = A_hard + B_gap_fill_rescued + C_network_coupled + D_anomaly,
         x_idx       = as.integer(organism))

# Long format for stacking
long <- per_org |>
  select(organism, x_idx,
         A_hard, B_gap_fill_rescued, C_network_coupled, D_anomaly) |>
  pivot_longer(c(A_hard, B_gap_fill_rescued, C_network_coupled, D_anomaly),
               names_to = "tier_raw", values_to = "n") |>
  mutate(tier = factor(case_when(
    tier_raw == "A_hard"             ~ "Tier A",
    tier_raw == "B_gap_fill_rescued" ~ "Tier B",
    tier_raw == "C_network_coupled"  ~ "Tier C",
    tier_raw == "D_anomaly"          ~ "Tier D"
  ), levels = c("Tier A","Tier B","Tier C","Tier D")))

# ---------- Colours ----------
tier_colors <- c("Tier A" = "#1d5a30",
                 "Tier B" = "#C9C9C9",
                 "Tier C" = "#B9D2B6",
                 "Tier D" = "#C04545")

class_colors <- c(
  Bacteroidia     = "#1F78B4", Clostridia    = "#E31A1C",
  Bacilli         = "#33A02C", Methanobacteria = "#FF7F00",
  Negativicutes   = "#6A3D9A", Coriobacteriia = "#B15928",
  Lentisphaeria   = "#666666", Kiritimatiellia = "#999999"
)

# ---------- Layout helpers ----------
N     <- nrow(per_org)
ref_n <- sum(per_org$group == "Reference")
arg_n <- sum(per_org$group == "ARG-MAG")
ymax  <- max(per_org$n_total) * 1.25

# Class strip y-band (below x-axis)
strip_top <- -ymax * 0.005
strip_bot <- -ymax * 0.030

# ---------- Plot ----------
p <- ggplot(long, aes(x = organism, y = n, fill = tier)) +
  geom_col(width = 0.8, colour = "white", linewidth = 0.2) +
  geom_text(data = per_org,
            aes(x = organism, y = n_total + ymax * 0.015,
                label = n_essential),
            inherit.aes = FALSE, size = 2.4,
            colour = "#1d5a30", fontface = "bold") +
  scale_fill_manual(values = tier_colors, name = "Tier") +
  geom_vline(xintercept = c(ref_n + 0.5, ref_n + arg_n + 0.5),
             linetype = "dashed", colour = "gray60", linewidth = 0.4) +
  # Section headers
  annotate("text", x = (1 + ref_n) / 2,                          y = ymax * 0.98,
           label = "Cultured references", fontface = "bold", size = 3.6) +
  annotate("text", x = (ref_n + 1 + ref_n + arg_n) / 2,           y = ymax * 1.02,
           label = "ARG-MAG", fontface = "bold", size = 3.6) +
  annotate("text", x = (ref_n + arg_n + 1 + N) / 2,               y = ymax * 0.98,
           label = "Pan-Draft species-representative",
           fontface = "bold", size = 3.6) +
  # Second fill scale: class colour strip
  ggnewscale::new_scale_fill() +
  geom_rect(data = per_org |> filter(!is.na(class)),
            aes(xmin = x_idx - 0.4, xmax = x_idx + 0.4,
                ymin = strip_bot,   ymax = strip_top,
                fill = class),
            inherit.aes = FALSE) +
  scale_fill_manual(values = class_colors, name = "Class", drop = FALSE) +
  scale_y_continuous(expand = expansion(mult = c(0.0, 0.02)),
                     limits = c(strip_bot * 1.2, ymax)) +
  coord_cartesian(clip = "off") +
  labs(x = NULL, y = "Number of organic compounds per model") +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor   = element_blank(),
    axis.text.x  = element_text(angle = 70, hjust = 1, size = 7.5),
    axis.text.y  = element_text(size = 8),
    plot.margin  = margin(t = 12, r = 30, b = 10, l = 10),

    # ===== Change legend placement here =====
    # A) outside, right (default — clean academic look):
    legend.position = "right",
    legend.box      = "vertical",

    # B) above plot, horizontal:
    # legend.position = "top",
    # legend.direction = "horizontal",
    # legend.box      = "vertical",

    # C) below plot, horizontal:
    # legend.position = "bottom",
    # legend.direction = "horizontal",
    # legend.box      = "vertical",

    legend.title = element_text(face = "bold"),
    legend.key.height = unit(0.8, "lines"),
    legend.spacing.y  = unit(0.2, "lines")
  )

ggsave(out_path, p, width = 14.5, height = 6.8, dpi = 200)
cat("Saved:", out_path, "\n")
