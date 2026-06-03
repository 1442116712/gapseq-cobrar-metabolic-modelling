# =============================================================
# Figure 6 — In-silico FAB rescue pre-validation (DSM 11370)
#
# Three panels in a single horizontal row (1 × 3 layout):
#   B (left)   — Cocktail leave-one-out (7 horizontal bars, 3 tiers)
#                — reveals Co²⁺ is the sole rescue-critical component
#   C (middle) — Key conditions (6 growth bars including M2-no-RF arms)
#                — confirms the LOO finding by adding Co²⁺ alone
#                  to FAB and to M2-no-RF
#   D (right)  — Co²⁺ titration on FAB-only (log–log dose–response)
#                — quantifies the rescue threshold (~10⁻⁶ mmol gDW⁻¹ h⁻¹)
#
# Output: figures/Fig6_FAB_rescue.png  (18 × 5.5 inches, 1000 dpi)
# Combine with Fig6_invitro_4cards.png (Panel A, 2 × 2 cards) and
# Fig6E_B12_pathway.png (Panel E) for the full Figure 6.
# =============================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(patchwork)
  library(scales)
})

# ---------- Output ----------
out_dir  <- "C:/Users/CFL/OneDrive - Queen's University Belfast/PhD_onedrive/Manuscript4/figures"
out_path <- file.path(out_dir, "Fig6_FAB_rescue.png")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ---------- Palette ----------
PAL <- c(
  critical   = "#B22222",   # red — rescue-critical / negative
  rate_limit = "#D4A017",   # gold — rate-limiting cofactors / FAB + Co²⁺
  redundant  = "#888888",   # grey — redundant in FAB
  positive   = "#2E7D32",   # green — positive control / rescued
  fab        = "#555555",   # dark grey — FAB-only (no growth)
  m2         = "#888888",   # medium grey — M2-no-rumen-only (no growth)
  cocktail   = "#1F4E79",   # navy — full cocktail / SOP
  m2_co      = "#C46A39"    # orange-brown — M2-no-rumen + Co²⁺ (partial rescue)
)

# =============================================================
# Panel B — Cocktail leave-one-out (7 components, 3 tiers)
# Establishes that Co²⁺ is the SOLE rescue-critical component
# of the 7-part cocktail — provides the rationale for testing
# Co²⁺ alone in Panel C.
# =============================================================
panel_B <- data.frame(
  component = c("Co²⁺", "Mn²⁺", "Cu²⁺", "Fe²⁺",
                "Zn²⁺", "Ca²⁺", "Thiamin"),
  growth    = c(0, 0.137, 0.137, 0.187, 0.274, 0.274, 0.274),
  tier      = c("critical", "rate_limit", "rate_limit", "rate_limit",
                "redundant", "redundant", "redundant"),
  stringsAsFactors = FALSE
)
panel_B <- panel_B[order(panel_B$growth), ]
panel_B$component <- factor(panel_B$component, levels = panel_B$component)
panel_B$tier_lab <- factor(c(critical = "Critical (rescue fails)",
                              rate_limit = "Rate-limiting",
                              redundant = "Redundant")[panel_B$tier],
                            levels = c("Critical (rescue fails)",
                                       "Rate-limiting", "Redundant"))

# Pre-compute label x-positions: nudge the "0" label right so it doesn't overlap the axis
panel_B$text_x <- ifelse(panel_B$growth < 0.001, 0.012, panel_B$growth)

p_B <- ggplot(panel_B, aes(x = growth, y = component, fill = tier_lab)) +
  geom_col(width = 0.7, colour = "white", linewidth = 1.0) +
  geom_text(aes(x = text_x,
                label = ifelse(growth < 0.001, "0", sprintf("%.3f", growth))),
            hjust = -0.15, size = 3.1, colour = "grey15") +
  geom_vline(xintercept = 0.274, linetype = "dashed",
             colour = PAL[["cocktail"]], linewidth = 0.5) +
  # Place "full cocktail = 0.274" ABOVE the topmost bar (Thiamin at y=7)
  annotate("text", x = 0.20, y = 7.55,
           label = "full cocktail = 0.274",
           colour = PAL[["cocktail"]], size = 3.0, hjust = 0, vjust = 0) +
  scale_fill_manual(values = c("Critical (rescue fails)" = PAL[["critical"]],
                                "Rate-limiting"          = PAL[["rate_limit"]],
                                "Redundant"              = PAL[["redundant"]]),
                    name = NULL) +
  scale_x_continuous(limits = c(0, 0.34), expand = c(0, 0)) +
  scale_y_discrete(expand = expansion(add = c(0.5, 1.0))) +
  coord_cartesian(clip = "off") +
  labs(x = expression(paste("Growth when component dropped (h"^-1, ")")),
       y = NULL) +
  theme_classic(base_size = 11) +
  theme(axis.text.y = element_text(size = 10, colour = "grey20"),
        axis.text.x = element_text(size = 10, colour = "grey20"),
        axis.title.x = element_text(size = 11),
        axis.line   = element_line(colour = "grey50"),
        axis.ticks  = element_line(colour = "grey50"),
        legend.position = c(0.80, 0.18),
        legend.background = element_rect(fill = "white",
                                         colour = "grey80",
                                         linewidth = 0.3),
        legend.key.size = unit(0.42, "cm"),
        legend.text = element_text(size = 8.5),
        legend.margin = margin(2, 4, 2, 4))

# =============================================================
# Panel C — Key conditions (6 growth bars)
# Validates Panel B's LOO finding by directly testing Co²⁺ alone:
#   positive control → two "fails" (FAB, M2-no-RF) →
#   full rescue (FAB+cocktail) → two single-Co²⁺ arms
#   (FAB+Co, M2+Co)
# =============================================================
panel_C <- data.frame(
  condition = c("gapseq min\n(pos ctrl)",
                "FAB\nonly",
                "M2 (no RF)\nonly",
                "FAB + full\ncocktail",
                "FAB +\nCo²⁺ only",
                "M2 (no RF) +\nCo²⁺ only"),
  growth    = c(0.204, 0, 0, 0.274, 0.137, 0.034),
  group     = c("positive", "fab", "m2", "cocktail", "rate_limit", "m2_co"),
  stringsAsFactors = FALSE
)
panel_C$condition <- factor(panel_C$condition, levels = panel_C$condition)
panel_C$lbl <- ifelse(panel_C$growth == 0,
                      "0\n(no growth)", sprintf("%.3f", panel_C$growth))
panel_C$lbl_y <- ifelse(panel_C$growth == 0, 0.012, panel_C$growth + 0.015)

p_C <- ggplot(panel_C, aes(x = condition, y = growth, fill = group)) +
  geom_col(width = 0.55, colour = "white", linewidth = 1.0) +
  geom_text(aes(y = lbl_y, label = lbl),
            colour = "grey15", size = 2.8, lineheight = 0.85) +
  geom_hline(yintercept = 1e-4, linetype = "dashed",
             colour = PAL[["critical"]], linewidth = 0.5) +
  annotate("label", x = 0.55, y = 0.30,
           label = "growth threshold = 10⁻⁴ h⁻¹",
           colour = PAL[["critical"]], fill = "white",
           label.size = 0.4, size = 2.8, hjust = 0,
           label.padding = unit(0.2, "lines"),
           label.r = unit(0.15, "lines")) +
  scale_fill_manual(values = PAL, guide = "none") +
  scale_y_continuous(limits = c(0, 0.33), expand = c(0, 0)) +
  labs(x = NULL, y = expression(paste("Growth rate (h"^-1, ")"))) +
  theme_classic(base_size = 11) +
  theme(axis.text.x = element_text(size = 7.5, colour = "grey20",
                                   lineheight = 0.9),
        axis.text.y = element_text(size = 10, colour = "grey20"),
        axis.title.y = element_text(size = 11),
        axis.line   = element_line(colour = "grey50"),
        axis.ticks  = element_line(colour = "grey50"))

# =============================================================
# Panel D — Co²⁺ titration on FAB-only (log–log dose–response)
# =============================================================
panel_D <- data.frame(
  flux   = c(1e-8, 1e-7, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1, 10),
  growth = c(1.37e-6, 1.37e-5, 1.37e-4, 1.37e-3, 1.37e-2,
             0.137, 0.137, 0.137, 0.137, 0.137),
  stringsAsFactors = FALSE
)
panel_D$status <- factor(ifelse(panel_D$growth > 1e-4, "rescued", "below"),
                          levels = c("rescued", "below"))

p_D <- ggplot(panel_D, aes(x = flux, y = growth)) +
  # Blood Co²⁺ band
  annotate("rect", xmin = 1e-6, xmax = 2e-6, ymin = 1e-7, ymax = 1.5,
           fill = PAL[["critical"]], alpha = 0.16) +
  annotate("text", x = sqrt(1e-6 * 2e-6), y = 0.5,
           label = "5 % horse blood\nCo²⁺ range",
           colour = PAL[["critical"]], size = 3.0,
           fontface = "italic", lineheight = 0.85) +
  # Rescue threshold
  geom_hline(yintercept = 1e-4, linetype = "dashed",
             colour = "grey55", linewidth = 0.5) +
  annotate("text", x = 8e-4, y = 1.6e-4,
           label = "rescue threshold 10⁻⁴ h⁻¹",
           colour = "grey40", size = 3.0, hjust = 0, vjust = 0) +
  # SOP cocktail line
  geom_vline(xintercept = 0.02, linetype = "dotted",
             colour = PAL[["cocktail"]], linewidth = 0.5) +
  annotate("text", x = 0.025, y = 1.5e-7,
           label = "SOP cocktail\nCo²⁺ (10 µM)",
           colour = PAL[["cocktail"]], size = 3.0,
           hjust = 0, lineheight = 0.85) +
  # Saturation line
  geom_hline(yintercept = 0.137, linetype = "dashed",
             colour = "grey55", linewidth = 0.4) +
  annotate("text", x = 0.5, y = 0.07,
           label = "saturation 0.137 h⁻¹",
           colour = "grey40", size = 3.0, hjust = 0.5, vjust = 1) +
  # Data
  geom_line(colour = PAL[["cocktail"]], linewidth = 0.5, alpha = 0.7) +
  geom_point(aes(colour = status), size = 2.8, stroke = 1.2) +
  scale_colour_manual(values = c(rescued = PAL[["positive"]],
                                  below   = PAL[["critical"]]),
                      labels = c(rescued = "rescued (growth > 10⁻⁴)",
                                 below   = "below threshold"),
                      name = NULL) +
  scale_x_log10(limits = c(3e-9, 30),
                breaks = 10^seq(-8, 1, 1),
                labels = trans_format("log10", math_format(10^.x))) +
  scale_y_log10(limits = c(1e-7, 1.5),
                breaks = 10^seq(-7, 0, 1),
                labels = trans_format("log10", math_format(10^.x))) +
  labs(x = expression(paste("Co"^"2+", " uptake flux (mmol gDW"^-1, " h"^-1, ")")),
       y = expression(paste("Growth rate (h"^-1, ")"))) +
  theme_classic(base_size = 11) +
  theme(axis.text   = element_text(size = 10, colour = "grey20"),
        axis.title  = element_text(size = 11),
        axis.line   = element_line(colour = "grey50"),
        axis.ticks  = element_line(colour = "grey50"),
        # Legend INSIDE the panel at top — avoids clipping by the figure boundary
        legend.position    = "top",
        legend.direction   = "horizontal",
        legend.justification = "right",
        legend.background  = element_rect(fill = "white",
                                          colour = "grey80",
                                          linewidth = 0.3),
        legend.key.size    = unit(0.4, "cm"),
        legend.text        = element_text(size = 9),
        legend.margin      = margin(2, 4, 2, 4),
        plot.margin        = margin(t = 8, r = 8, b = 8, l = 8))

# =============================================================
# Compose 1 × 3 horizontal layout
# Order:  B (LOO)  ->  C (key conditions)  ->  D (titration)
# Rationale: LOO discovers Co²⁺ as critical (B); key conditions
# panel then confirms by adding Co²⁺ alone (C); titration
# quantifies the threshold (D).
# =============================================================
fig <- p_B | p_C | p_D

ggsave(out_path, fig,
       width = 18, height = 5.5, units = "in",
       dpi = 1000, bg = "white")

cat(sprintf("Wrote: %s\n", out_path))
