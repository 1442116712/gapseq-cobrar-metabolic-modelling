# =============================================================
# Figure 6E — B12 biosynthesis pathway in the DSM 11370 model
#
# Standalone template for the bottom-right panel of Figure 6
# (cobalamin / B12 biosynthesis schematic).
#
# Layout:
#   Main horizontal flow (left → right):
#     Co²⁺(ext)  →  Co²⁺(cyto)  →  adenosyl-cobyric acid
#                                  →  adenosyl-cobinamide
#                                  →  cob(I)alamin
#                                  →  AdoCbl (active B12)
#                                  →  biomass
#   Top branch (blocked):
#     External B12 (blood / yeast extract)
#         ⤬→ AdoCbl     (NO transporter in the model)
#
# Reaction IDs sit above each arrow. Node fills group concepts:
#   yellow = cobalt ion pools
#   blue   = cobalamin biosynthesis intermediates
#   green  = active product / biomass
#   red    = blocked external pool (no exchange in the model)
#
# Output: figures/Fig6E_B12_pathway.png  (12 × 6 inches, 600 dpi)
# =============================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(tibble)
  library(dplyr)
})

# ---------- Paths ----------
out_dir  <- "C:/Users/CFL/OneDrive - Queen's University Belfast/PhD_onedrive/Manuscript4/figures"
out_path <- file.path(out_dir, "Fig6E_B12_pathway.png")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ---------- Palette ----------
node_fill <- c(
  ion          = "#FFE082",   # cobalt pools (yellow)
  intermediate = "#BBDEFB",   # cobalamin biosynthesis intermediates (blue)
  product      = "#A5D6A7",   # AdoCbl, active B12 (green)
  biomass      = "#C5E1A5",   # biomass node (light green)
  blocked      = "#FFCDD2"    # external B12 with no transporter (red)
)
arrow_color_ok      <- "grey40"
arrow_color_blocked <- "#B22222"

# ---------- Nodes ----------
# x increases left → right; y = 0 is the main flow row;
# y > 0 puts a node above the main flow.
nodes <- tribble(
  ~id,           ~x,    ~y,   ~label,                                  ~type,
  "Co_ext",      0.0,   0.0,  "Co²⁺\n(medium / blood)",      "ion",
  "Co_in",       2.4,   0.0,  "Co²⁺\n(cytoplasm)",           "ion",
  "cobyric",     4.8,   0.0,  "adenosyl-\ncobyric acid",               "intermediate",
  "cobinamide",  7.2,   0.0,  "adenosyl-\ncobinamide",                 "intermediate",
  "cobI",        9.6,   0.0,  "cob(I)alamin",                          "intermediate",
  "AdoCbl",     12.0,   0.0,  "AdoCbl\n(active B12)",                  "product",
  "biomass",    14.4,   0.0,  "biomass",                               "biomass",
  # B12_ext positioned DIRECTLY ABOVE AdoCbl so the blocked arrow is a clean vertical line
  "B12_ext",    12.0,   2.4,  "external B12\n(blood / yeast extract)", "blocked"
)
nodes <- nodes %>%
  mutate(type = factor(type, levels = names(node_fill)))

# ---------- Edges ----------
edges <- tribble(
  ~from,         ~to,           ~label,                       ~ok,
  "Co_ext",      "Co_in",       "uptake",                     TRUE,
  "Co_in",       "cobyric",     "Co insertion\n→ corrin",TRUE,
  "cobyric",     "cobinamide",  "rxn04384\nrxn05054",         TRUE,
  "cobinamide",  "cobI",        "rxn05029",                   TRUE,
  "cobI",        "AdoCbl",      "rxn08194\nrxn10091",         TRUE,
  "AdoCbl",      "biomass",     "biomass\ndemand",            TRUE,
  "B12_ext",     "AdoCbl",      "NO transporter\nin the model", FALSE
)

# Resolve coordinates for each edge endpoint
edges <- edges %>%
  left_join(nodes %>% transmute(from = id, x1 = x, y1 = y), by = "from") %>%
  left_join(nodes %>% transmute(to   = id, x2 = x, y2 = y), by = "to") %>%
  mutate(midx = (x1 + x2) / 2,
         midy = (y1 + y2) / 2)

# Working arrows: same fixed half-length, centred at the midpoint between
# consecutive node centres. This makes all gray arrows identical visually
# regardless of how wide the label boxes are.
ARROW_HALF_LEN <- 0.55
edges <- edges %>%
  mutate(
    ax1 = ifelse(ok, midx - ARROW_HALF_LEN, x1),
    ay1 = ifelse(ok, y1,                   y1 - 0.55),
    ax2 = ifelse(ok, midx + ARROW_HALF_LEN, x2),
    ay2 = ifelse(ok, y2,                   y2 + 0.45)
  )

# ---------- Plot ----------
p <- ggplot() +
  # 1. Successful arrows (grey solid) — all uniform length
  geom_segment(
    data = edges %>% filter(ok),
    aes(x = ax1, y = ay1, xend = ax2, yend = ay2),
    arrow = arrow(length = unit(3, "mm"), type = "closed"),
    colour = arrow_color_ok, linewidth = 0.6
  ) +
  # 2. Blocked arrow (red dashed) — clean vertical line from B12_ext down to AdoCbl
  geom_segment(
    data = edges %>% filter(!ok),
    aes(x = ax1, y = ay1, xend = ax2, yend = ay2),
    arrow = arrow(length = unit(4, "mm"), type = "closed"),
    colour = arrow_color_blocked, linewidth = 0.7, linetype = "dashed"
  ) +
  # 3. X mark sitting on the blocked vertical line, mid-way down
  geom_point(
    data = edges %>% filter(!ok),
    aes(x = midx, y = midy),
    shape = 4, size = 11, colour = arrow_color_blocked, stroke = 2.2
  ) +
  # 4. Edge labels — working arrows (above arrow)
  geom_text(
    data = edges %>% filter(ok, label != ""),
    aes(x = midx, y = midy + 0.5, label = label),
    size = 2.8, colour = "grey30", lineheight = 0.85
  ) +
  # 5. Edge label — blocked arrow (bold red, placed to the RIGHT of X mark)
  geom_text(
    data = edges %>% filter(!ok),
    aes(x = midx + 0.6, y = midy, label = label),
    size = 3.0, colour = arrow_color_blocked, fontface = "bold",
    lineheight = 0.85, hjust = 0
  ) +
  # 6. Node labels
  geom_label(
    data = nodes,
    aes(x = x, y = y, label = label, fill = type),
    size = 2.9, label.size = 0.4, lineheight = 0.9,
    label.padding = unit(0.32, "lines"),
    label.r       = unit(0.18, "lines"),
    colour = "grey15", fontface = "bold"
  ) +
  scale_fill_manual(values = node_fill, guide = "none") +
  # 7. Footer note (compact, two lines close to the flow)
  annotate("text", x = 7.0, y = -1.0,
           label = paste0(
             "5 cobalamin biosynthesis reactions present in the reconstruction; ",
             "0 cobalamin exchange reactions encoded"),
           size = 3.0, colour = "grey25", fontface = "italic") +
  annotate("text", x = 7.0, y = -1.4,
           label = paste0(
             "→ AdoCbl produced exclusively by de novo synthesis ",
             "from imported Co²⁺"),
           size = 3.0, colour = "grey25", fontface = "italic") +
  coord_cartesian(xlim = c(-1.1, 15.1), ylim = c(-1.7, 3.0)) +
  theme_void()

# ---------- Save ----------
ggsave(out_path, p,
       width = 12, height = 5, dpi = 600, bg = "white")

cat(sprintf("Wrote: %s\n", out_path))
