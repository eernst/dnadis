# common.R — shared setup for assembly and comparison report templates
#
# Sourced at the top of each Rmd's setup chunk.  Provides:
#   - Package loading (base set; callers may p_load additional packages)
#   - Font detection and showtext setup
#   - Shared theming constants (base_font_pt, axis_theme, apply_legend_theme)
#   - Okabe-Ito palette, oi_rgba(), bar_cell()
#   - Classification colors / labels / levels
#   - Helper functions: save_pdf, girafe_styled, gt_theme, fmt_mb
#   - setup_subgenomes(): compute subgenome color palettes from ref_lengths

# ---------------------------------------------------------------------------
# Packages
# ---------------------------------------------------------------------------
if (!requireNamespace("pacman", quietly = TRUE))
  install.packages("pacman", repos = "https://cloud.r-project.org")
library(pacman)
pacman::p_load(readr, dplyr, stringr, tidyr, tibble,
               ggplot2, ggnewscale, patchwork, grid,
               ggiraph, htmlwidgets, ggrepel,
               ggokabeito, colorspace,
               gt, gtExtras, scales,
               showtext, sysfonts, knitr)

# ---------------------------------------------------------------------------
# Font setup
# ---------------------------------------------------------------------------
.resolve_font <- function(family) {
  slug <- tolower(gsub(" ", "", family))
  path <- systemfonts::match_fonts(family)$path
  if (nzchar(path) && grepl(slug, tolower(gsub("[- ]", "", basename(path)))))
    return(path)
  NULL
}
base_family <- "sans"
for (.fam in c("Helvetica Neue", "Helvetica", "Liberation Sans")) {
  .r <- .resolve_font(.fam)
  if (!is.null(.r)) {
    sysfonts::font_add(.fam,
      regular    = .r,
      bold       = systemfonts::match_fonts(.fam, weight = "bold")$path,
      italic     = systemfonts::match_fonts(.fam, italic = TRUE)$path,
      bolditalic = systemfonts::match_fonts(.fam, italic = TRUE, weight = "bold")$path
    )
    base_family <- .fam
    break
  }
}
showtext_auto()

# ---------------------------------------------------------------------------
# Shared theming
# ---------------------------------------------------------------------------
base_font_pt <- 8
axis_title_margin_pt <- 4
axis_text_margin_pt <- 4

axis_theme <- theme(
  axis.title.x = element_text(size = base_font_pt, family = base_family,
                               margin = margin(t = axis_title_margin_pt)),
  axis.title.y = element_text(size = base_font_pt, family = base_family,
                               margin = margin(r = axis_title_margin_pt)),
  axis.text.x = element_text(size = base_font_pt, family = base_family,
                               margin = margin(t = axis_text_margin_pt)),
  axis.text.y = element_text(size = base_font_pt, family = base_family)
)

apply_legend_theme <- function(p, text_pt = base_font_pt, key_pt = base_font_pt,
                               tight = TRUE) {
  base <- theme(
    legend.title = element_text(size = text_pt, family = base_family),
    legend.text  = element_text(size = text_pt, family = base_family),
    legend.key.height = grid::unit(key_pt, "pt"),
    legend.key.width  = grid::unit(key_pt, "pt"),
    legend.direction = "horizontal",
    legend.box = "horizontal"
  )
  if (tight) {
    base <- base + theme(
      legend.box.margin = margin(0, 0, 0, 0),
      legend.margin = margin(0, 0, 0, 0),
      legend.spacing.y = grid::unit(0, "pt"),
      legend.spacing.x = grid::unit(2, "pt")
    )
  }
  p + base
}

# ---------------------------------------------------------------------------
# Okabe-Ito palette and classification colors
# ---------------------------------------------------------------------------
oi <- palette_okabe_ito()

oi_rgba <- function(col, alpha = 0.25) {
  v <- col2rgb(col)
  sprintf("rgba(%d,%d,%d,%.2f)", v[1], v[2], v[3], alpha)
}

bar_cell <- function(pct, col, label, alpha = 0.25) {
  paste0(
    '<div style="position:relative;padding:2px 6px;white-space:nowrap;border-radius:2px;">',
    sprintf('<div style="position:absolute;inset:0;width:%.1f%%;background-color:%s;border-radius:2px;z-index:0;"></div>', pct, oi_rgba(col, alpha)),
    '<span style="position:relative;z-index:1;">', label, '</span></div>'
  )
}

classification_colors <- c(
  "chrom_assigned"     = oi[5],
  "chrom_unassigned"   = oi[2],
  "organelle_complete" = oi[3],
  "organelle_debris"   = oi[4],
  "rDNA"               = oi[1],
  "contaminant"        = oi[7],
  "chrom_debris"       = lighten(oi[5], amount = 0.5),
  "debris"             = oi[9],
  "unclassified"       = oi[8]
)

classification_labels <- c(
  "chrom_assigned"     = "Chromosome (assigned)",
  "chrom_unassigned"   = "Chromosome (unassigned)",
  "organelle_complete" = "Organelle (complete)",
  "organelle_debris"   = "Organelle (debris)",
  "rDNA"               = "rDNA",
  "contaminant"        = "Contaminant",
  "chrom_debris"       = "Chromosome debris",
  "debris"             = "Debris",
  "unclassified"       = "Unclassified"
)

clf_labels_short <- c(
  "chrom_assigned"     = "Chr assigned",
  "chrom_unassigned"   = "Chr unassigned",
  "organelle_complete" = "Organelle",
  "organelle_debris"   = "Organelle debris",
  "rDNA"               = "rDNA",
  "contaminant"        = "Contaminant",
  "chrom_debris"       = "Chr debris",
  "debris"             = "Debris",
  "unclassified"       = "Unclassified"
)

classification_levels <- c(
  "chrom_assigned", "chrom_unassigned", "organelle_complete", "organelle_debris",
  "rDNA", "contaminant", "chrom_debris", "debris", "unclassified"
)

# ---------------------------------------------------------------------------
# Helper functions
# ---------------------------------------------------------------------------
save_pdf <- function(path, plot, width = 7.2, height = 5, dpi = 300) {
  dev <- if (capabilities("cairo")) cairo_pdf else "pdf"
  invisible(ggsave(path, plot, width = width, height = height, dpi = dpi, device = dev))
}

girafe_styled <- function(ggobj, width_svg = 10, height_svg = 6.9) {
  g <- ggiraph::girafe(ggobj = ggobj, width_svg = width_svg, height_svg = height_svg)
  ggiraph::girafe_options(g,
    ggiraph::opts_tooltip(css = paste0(
      "font-family: ", base_family, "; font-size: 9pt; ",
      "background: white; padding: 6px; border: 1px solid #666;"
    )),
    ggiraph::opts_hover(css = "opacity: 0.8;"),
    ggiraph::opts_toolbar(saveaspng = FALSE)
  )
}

gt_theme <- function(gt_obj) {
  gt_obj %>%
    tab_options(
      table.font.size = px(11),
      column_labels.font.weight = "bold",
      column_labels.border.top.style = "none",
      column_labels.border.bottom.width = px(2),
      column_labels.border.bottom.color = "#333",
      table.border.top.style = "none",
      table.border.bottom.style = "none",
      table_body.border.bottom.color = "#999",
      data_row.padding = px(4),
      heading.border.bottom.style = "none"
    ) %>%
    opt_row_striping(row_striping = FALSE)
}

fmt_mb <- function(x) sprintf("%.1f", x / 1e6)

# ---------------------------------------------------------------------------
# Subgenome detection and color setup
# ---------------------------------------------------------------------------
setup_subgenomes <- function(ref_lengths) {
  sg_levels <- ref_lengths %>%
    distinct(subgenome) %>% pull(subgenome) %>% as.character()
  sg_levels <- sg_levels[grepl("^[A-Z]$", sg_levels)]
  sg_levels <- sort(sg_levels)
  n_sg <- length(sg_levels)
  if (n_sg > 8) { sg_levels <- sg_levels[1:8]; n_sg <- 8 }
  has_subgenomes <- (n_sg >= 2)

  pal_sg_dark  <- c(oi[5], oi[1], oi[3], oi[6], oi[7], oi[4], oi[2], oi[8])
  pal_sg_light <- lighten(pal_sg_dark, amount = 0.5)
  plot_levels  <- if (has_subgenomes) sg_levels else c("G")

  sg_cols_dark  <- setNames(pal_sg_dark[seq_along(plot_levels)],  plot_levels)
  sg_cols_light <- setNames(pal_sg_light[seq_along(plot_levels)], plot_levels)
  col_dark  <- c(sg_cols_dark,  "NA" = pal_sg_dark[1])
  col_light <- c(sg_cols_light, "NA" = pal_sg_light[1])

  list(
    sg_levels = sg_levels,
    n_sg = n_sg,
    has_subgenomes = has_subgenomes,
    pal_sg_dark = pal_sg_dark,
    pal_sg_light = pal_sg_light,
    plot_levels = plot_levels,
    sg_cols_dark = sg_cols_dark,
    sg_cols_light = sg_cols_light,
    col_dark = col_dark,
    col_light = col_light
  )
}
