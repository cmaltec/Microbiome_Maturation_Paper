#________________________________________
# - File Name : palette.R
# - Purpose : Shared colour palette and theme_paper2() ggplot theme
# - Created By : Christian Maltecca
# - e-mail: cmaltec@ncsu.edu
#      / \       _-'
#     /|  \-''- _ /
# ' { |          \
#   /             \
#  /      "o.  |o }
#     |           \
#                 ',
#      \_         __\
#        ''-_    \.//
#          / '-____'
#         /
#        '
#     -'
# ╔╗╔╔═╗╔═╗╦ ╦
# ║║║║  ╚═╗║ ║
# ╝╚╝╚═╝╚═╝╚═╝
# Mon May 25 10:09:00 2026 ------------------------------

pal <- list(
  SEN      = "#C2785C",
  TOL      = "#5B8FA8",
  neutral  = "#4A4A4A",
  c2       = "#D4A84B",   # litter / c² component
  h2       = "#234A72",   # additive / h² — must stay distinct from TOL
  ref_line = "#9A8F82"    # reference / identity lines
)

line_colors <- c(SEN = pal$SEN, TOL = pal$TOL)
line_fills  <- c(SEN = pal$SEN, TOL = pal$TOL)


# Canonical ggplot2 theme for all paper 2 figures.
# Soft grey panel + white gridlines, thin grey axis lines, no panel border.
# Title/subtitle flush left; panel-letter tags via plot_annotation(tag_levels).
# Base font size 10 by default.
theme_paper2 <- function(base_size = 10) {
  theme_bw(base_size = base_size) %+replace%
    theme(
      # Panel surface
      panel.background = element_rect(fill = "grey98", color = NA),
      panel.border     = element_blank(),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color = "white", linewidth = 0.55),
      # Axes
      axis.line    = element_line(color = "grey75", linewidth = 0.4),
      axis.ticks   = element_line(color = "grey75", linewidth = 0.4),
      axis.title.x = element_text(margin = margin(t = 8)),
      axis.title.y = element_text(margin = margin(r = 8), angle = 90),
      # Strips (for facets)
      strip.background = element_rect(fill = "#F5F0EB", color = NA),
      strip.text       = element_text(face = "bold", color = "grey15"),
      # Legend
      legend.background = element_rect(fill = alpha("white", 0.85), color = NA),
      legend.key        = element_rect(fill = "transparent", color = NA),
      # Titles & subtitles — flush left, panel-letter via tags
      plot.title          = element_text(size = base_size + 1, face = "bold",
                                         color = "grey15", hjust = 0,
                                         margin = margin(b = 4)),
      plot.subtitle       = element_text(size = base_size - 1.5, color = "grey25",
                                         hjust = 0, margin = margin(b = 8)),
      plot.title.position = "plot",
      # Panel letter tags (A, B, ...) via plot_annotation(tag_levels = "A")
      plot.tag          = element_text(face = "bold", size = base_size + 3,
                                       color = "grey15"),
      # Plot margins
      plot.margin = margin(t = 8, r = 12, b = 8, l = 12)
    )
}
