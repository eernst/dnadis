#!/usr/bin/env Rscript
# Classification summary bar chart for final_finalizer
# Waterfall design: top bar shows full assembly, bottom bar zooms into non-chromosome categories

if (!requireNamespace("pacman", quietly = TRUE)) {
  install.packages("pacman", repos = "https://cloud.r-project.org")
}
library(pacman)
pacman::p_load(
  readr, dplyr, stringr, ggplot2, tibble, tidyr, scales,
  ggiraph, htmlwidgets, ggrepel, ggokabeito, showtext, sysfonts
)

# Register Liberation Sans for consistent cross-platform rendering
sysfonts::font_add("Liberation Sans",
  regular    = "/usr/share/fonts/liberation-sans/LiberationSans-Regular.ttf",
  bold       = "/usr/share/fonts/liberation-sans/LiberationSans-Bold.ttf",
  italic     = "/usr/share/fonts/liberation-sans/LiberationSans-Italic.ttf",
  bolditalic = "/usr/share/fonts/liberation-sans/LiberationSans-BoldItalic.ttf"
)
showtext_auto()

# Placeholders - replaced by Python
summary_file   <- "__SUMMARY__"
out_pdf        <- "__OUTPDF__"
out_html       <- "__OUTHTML__"
plot_html      <- as.logical("__PLOTHTML__")
plot_suffix    <- "__SUFFIX__"
assembly_name  <- "__ASMNAME__"
reference_name <- "__REFNAME__"

# Font size hierarchy
base_family <- "Liberation Sans"
bar_label_font_pt <- 8    # top bar label (Chr assigned)
label_font_pt <- 6.5      # bottom bar labels
axis_font_pt <- 7         # axis and annotation text
annot_font_pt <- 7        # "Remaining X%" annotation

# Read contig summary
df <- read_tsv(summary_file, show_col_types = FALSE,
               col_types = cols(
                 original_name = col_character(),
                 classification = col_character(),
                 length = col_double()
               ))

# Check if we have any data
if (nrow(df) == 0) {
  message("No contigs in summary. Skipping classification bar plot.")
  quit(status = 0)
}

# Define classification order and colors
classification_levels <- c(
  "chrom_assigned",
  "chrom_unassigned",
  "organelle_complete",
  "organelle_debris",
  "rdna",
  "contaminant",
  "chrom_debris",
  "debris",
  "unclassified"
)

classification_labels <- c(
  "chrom_assigned" = "Chr assigned",
  "chrom_unassigned" = "Chr unassigned",
  "organelle_complete" = "Organelle",
  "organelle_debris" = "Organelle debris",
  "rdna" = "rDNA",
  "contaminant" = "Contaminant",
  "chrom_debris" = "Chr debris",
  "debris" = "Debris",
  "unclassified" = "Unclassified"
)

# Okabe-Ito colorblind-friendly palette
oi <- palette_okabe_ito()
classification_colors <- c(
  "chrom_assigned"     = oi[5],  # blue
  "chrom_unassigned"   = oi[7],  # reddish purple
  "organelle_complete" = oi[3],  # bluish green
  "organelle_debris"   = oi[4],  # yellow
  "rdna"               = oi[1],  # orange
  "contaminant"        = oi[6],  # vermillion
  "chrom_debris"       = oi[2],  # sky blue
  "debris"             = oi[9],  # gray
  "unclassified"       = oi[8]   # dark gray (black)
)

# Aggregate by classification
df_agg <- df %>%
  mutate(
    classification = factor(classification, levels = classification_levels),
    length_mb = length / 1e6
  ) %>%
  group_by(classification) %>%
  summarise(
    n_contigs = n(),
    total_mb = sum(length_mb, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  filter(!is.na(classification))

# Calculate overall totals
assembly_total_contigs <- sum(df_agg$n_contigs)
assembly_total_mb <- sum(df_agg$total_mb)

# Calculate percentages for full assembly
# Order: chrom_assigned first, then remaining categories by desc(total_mb)
# This matches the bottom bar order so connector ribbons don't cross
df_agg <- df_agg %>%
  mutate(
    pct_contigs = n_contigs / assembly_total_contigs * 100,
    pct_mb = total_mb / assembly_total_mb * 100
  ) %>%
  arrange(classification != "chrom_assigned", desc(total_mb)) %>%
  mutate(
    xmax = cumsum(pct_mb),
    xmin = xmax - pct_mb,
    x_center = (xmin + xmax) / 2
  )

# Build subtitle: "929.2 Mb · n=455 contigs", optionally prefixed with assembly name
subtitle_parts <- c(
  sprintf("%.1f Mb", assembly_total_mb),
  paste0("n=", assembly_total_contigs, " contigs")
)
subtitle_text <- paste(subtitle_parts, collapse = " \u00b7 ")
if (nzchar(assembly_name)) {
  subtitle_text <- paste0(assembly_name, " \u00b7 ", subtitle_text)
}

# ============================================
# TOP BAR: All categories with proper colors, only label chr_assigned
# ============================================
bar1_ymin <- 0.665
bar1_ymax <- 0.90

# Top bar uses all categories from df_agg with their proper colors
df_bar1 <- df_agg %>%
  filter(pct_mb > 0) %>%
  mutate(
    tooltip = paste0(
      classification_labels[as.character(classification)],
      "\n", sprintf("%.1f Mb (%.1f%%)", total_mb, pct_mb),
      "\nn=", n_contigs, " (", sprintf("%.1f%%", pct_contigs), ")"
    )
  )

# Chr assigned data (for label and connector)
df_chr <- df_agg %>%
  filter(classification == "chrom_assigned" & pct_mb > 0)

# Other categories (for bottom bar zoom)
df_other <- df_agg %>%
  filter(classification != "chrom_assigned" & pct_mb > 0)

if (nrow(df_chr) > 0) {
  chr_xmax <- df_chr$xmax[1]
} else {
  chr_xmax <- 0
}

# Chr assigned label
df_chr_label <- df_chr %>%
  mutate(
    label_text = paste0(
      classification_labels[as.character(classification)], "\n",
      sprintf("%.1f Mb (%.1f%%)", total_mb, pct_mb), "\n",
      "n=", n_contigs, " (", sprintf("%.1f%%", pct_contigs), ")"
    )
  )

# ============================================
# BOTTOM BAR: Zoom into non-chr_assigned categories
# ============================================
bar2_ymin <- 0.1925
bar2_ymax <- 0.2675

if (nrow(df_other) > 0) {
  other_total_mb <- sum(df_other$total_mb)
  other_pct <- sum(df_other$pct_mb)

  df_bar2 <- df_other %>%
    mutate(
      pct_within_other = total_mb / other_total_mb * 100
    ) %>%
    arrange(desc(total_mb)) %>%
    mutate(
      xmax2 = cumsum(pct_within_other),
      xmin2 = xmax2 - pct_within_other,
      x_center2 = (xmin2 + xmax2) / 2
    ) %>%
    mutate(
      tooltip = paste0(
        classification_labels[as.character(classification)],
        "\n", sprintf("%.1f Mb (%.1f%%)", total_mb, pct_mb), " of assembly",
        "\nn=", n_contigs, " (", sprintf("%.1f%%", pct_contigs), " of assembly)"
      ),
      label_text = paste0(
        classification_labels[as.character(classification)], "\n",
        sprintf("%.1f Mb (%.1f%%)", total_mb, pct_mb), "\n",
        "n=", n_contigs, " (", sprintf("%.1f%%", pct_contigs), ")"
      )
    )

  # Per-category connector ribbons from top bar to bottom bar
  # Each ribbon connects a segment's position in the top bar to its position in the bottom bar
  y_top <- bar1_ymin
  y_bot <- bar2_ymax
  n_curve_pts <- 20
  t_vals <- seq(0, 1, length.out = n_curve_pts)

  connector_polys <- bind_rows(lapply(seq_len(nrow(df_bar2)), function(i) {
    row <- df_bar2[i, ]
    cls <- as.character(row$classification)

    # Top bar position (from df_agg coordinates)
    top_left  <- row$xmin
    top_right <- row$xmax

    # Bottom bar position (zoomed coordinates)
    bot_left  <- row$xmin2
    bot_right <- row$xmax2

    # Smooth S-curve (ease-in-out) so edges leave top bar and arrive at bottom bar horizontally
    s <- 0.5 - 0.5 * cos(pi * t_vals)  # 0→1 sigmoid
    left_x  <- top_left  + s * (bot_left  - top_left)
    right_x <- top_right + s * (bot_right - top_right)
    y_seq   <- y_top + t_vals * (y_bot - y_top)

    tibble(
      x = c(left_x, rev(right_x)),
      y = c(y_seq, rev(y_seq)),
      classification = cls
    )
  }))
} else {
  df_bar2 <- tibble()
  other_pct <- 0
  connector_polys <- tibble()
}

# ============================================
# Shared ggrepel parameters for bottom bar labels
# ============================================
repel_params <- list(
  size = label_font_pt / .pt,
  family = base_family,
  lineheight = 0.95,
  direction = "x",
  nudge_y = -0.5,
  segment.color = "grey70",
  segment.size = 0.25,
  segment.curvature = -1e-20,
  segment.shape = 1,
  segment.ncp = 1,
  segment.angle = 90,
  segment.inflect = TRUE,
  min.segment.length = 0,
  box.padding = 0.25,
  point.padding = 0.75,
  point.size = 2,
  force = 1.5,
  force_pull = 0,
  max.overlaps = Inf,
  ylim = c(-Inf, bar2_ymin - 0.02),
  seed = 42
)

# ============================================
# BUILD PLOT
# ============================================
plot_height <- 3.0
y_min <- -0.35
y_max <- 1.0

p <- ggplot() +
  # Connector ribbons (draw first so bars overlay edges)
  {
    if (nrow(df_bar2) > 0) {
      geom_polygon(
        data = connector_polys,
        aes(x = x, y = y, fill = classification, group = classification),
        alpha = 0.25,
        color = NA
      )
    }
  } +
  # TOP BAR segments (all categories with proper colors, white borders)
  geom_rect(
    data = df_bar1,
    aes(xmin = xmin, xmax = xmax, ymin = bar1_ymin, ymax = bar1_ymax,
        fill = classification),
    color = "white",
    linewidth = 0.3
  ) +
  # Chr assigned label (inside the bar)
  {
    if (nrow(df_chr_label) > 0) {
      geom_text(
        data = df_chr_label,
        aes(x = x_center, y = (bar1_ymin + bar1_ymax) / 2, label = label_text),
        size = bar_label_font_pt / .pt,
        family = base_family,
        lineheight = 0.95,
        color = "white",
        fontface = "plain"
      )
    }
  } +
  # Tick marks at 25%, 50%, 75%
  annotate("segment", x = c(25, 50, 75), xend = c(25, 50, 75),
           y = bar1_ymax, yend = bar1_ymax + 0.015,
           color = "gray80", linewidth = 0.3) +
  # Axis labels at top
  annotate("text", x = 0, y = bar1_ymax + 0.03, label = "0%",
           size = axis_font_pt / .pt, family = base_family, hjust = 0) +
  annotate("text", x = 100, y = bar1_ymax + 0.03, label = "100%",
           size = axis_font_pt / .pt, family = base_family, hjust = 1) +
  # "Remaining X%" annotation (left-aligned)
  {
    if (nrow(df_bar2) > 0) {
      annotate("text", x = 0, y = bar2_ymax + 0.08,
               label = paste0("Remaining ", sprintf("%.1f%%", other_pct), " of assembly:"),
               size = annot_font_pt / .pt, family = base_family, fontface = "italic",
               color = "grey30", hjust = 0)
    }
  } +
  # BOTTOM BAR: colored segments
  {
    if (nrow(df_bar2) > 0) {
      geom_rect(
        data = df_bar2,
        aes(xmin = xmin2, xmax = xmax2, ymin = bar2_ymin, ymax = bar2_ymax,
            fill = classification),
        color = NA
      )
    }
  } +
  # Labels with ggrepel
  {
    if (nrow(df_bar2) > 0) {
      do.call(ggrepel::geom_text_repel, c(
        list(
          data = df_bar2,
          mapping = aes(x = x_center2, y = bar2_ymin, label = label_text)
        ),
        repel_params
      ))
    }
  } +
  scale_fill_manual(
    values = classification_colors,
    labels = classification_labels,
    drop = FALSE
  ) +
  scale_x_continuous(
    expand = c(0.01, 0.01),
    limits = c(0, 100)
  ) +
  scale_y_continuous(
    limits = c(y_min, y_max),
    expand = c(0, 0)
  ) +
  coord_cartesian(clip = "off") +
  theme_void(base_family = base_family) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 11, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = 9, color = "grey40"),
    legend.position = "none",
    plot.margin = margin(t = 10, r = 15, b = 15, l = 15)
  ) +
  labs(
    title = "Contig classification summary",
    subtitle = subtitle_text
  )

# Save PDF
ggsave(out_pdf, plot = p, width = 7.2, height = plot_height, units = "in", dpi = 300, device = cairo_pdf)
message("Classification summary bar plot saved to: ", out_pdf)

# Interactive HTML version
if (plot_html) {
  p_html <- ggplot() +
    # Connector ribbons
    {
      if (nrow(df_bar2) > 0) {
        geom_polygon(
          data = connector_polys,
          aes(x = x, y = y, fill = classification, group = classification),
          alpha = 0.25,
          color = NA
        )
      }
    } +
    # TOP BAR segments (interactive, all categories, white borders)
    ggiraph::geom_rect_interactive(
      data = df_bar1,
      aes(xmin = xmin, xmax = xmax, ymin = bar1_ymin, ymax = bar1_ymax,
          fill = classification, tooltip = tooltip),
      color = "white",
      linewidth = 0.3
    ) +
    # Chr assigned label
    {
      if (nrow(df_chr_label) > 0) {
        geom_text(
          data = df_chr_label,
          aes(x = x_center, y = (bar1_ymin + bar1_ymax) / 2, label = label_text),
          size = bar_label_font_pt / .pt,
          family = base_family,
          lineheight = 0.95,
          color = "white",
          fontface = "plain"
        )
      }
    } +
    # Tick marks at 25%, 50%, 75%
    annotate("segment", x = c(25, 50, 75), xend = c(25, 50, 75),
             y = bar1_ymax, yend = bar1_ymax + 0.015,
             color = "gray80", linewidth = 0.3) +
    # Axis labels
    annotate("text", x = 0, y = bar1_ymax + 0.03, label = "0%",
             size = axis_font_pt / .pt, family = base_family, hjust = 0) +
    annotate("text", x = 100, y = bar1_ymax + 0.03, label = "100%",
             size = axis_font_pt / .pt, family = base_family, hjust = 1) +
    # Remaining annotation (left-aligned)
    {
      if (nrow(df_bar2) > 0) {
        annotate("text", x = 0, y = bar2_ymax + 0.08,
                 label = paste0("Remaining ", sprintf("%.1f%%", other_pct), " of assembly:"),
                 size = annot_font_pt / .pt, family = base_family, fontface = "italic",
                 color = "grey30", hjust = 0)
      }
    } +
    # BOTTOM BAR (interactive)
    {
      if (nrow(df_bar2) > 0) {
        ggiraph::geom_rect_interactive(
          data = df_bar2,
          aes(xmin = xmin2, xmax = xmax2, ymin = bar2_ymin, ymax = bar2_ymax,
              fill = classification, tooltip = tooltip),
          color = NA
        )
      }
    } +
    # Labels with ggrepel (same params as PDF)
    {
      if (nrow(df_bar2) > 0) {
        do.call(ggrepel::geom_text_repel, c(
          list(
            data = df_bar2,
            mapping = aes(x = x_center2, y = bar2_ymin, label = label_text)
          ),
          repel_params
        ))
      }
    } +
    scale_fill_manual(
      values = classification_colors,
      labels = classification_labels,
      drop = FALSE
    ) +
    scale_x_continuous(
      expand = c(0.01, 0.01),
      limits = c(0, 100)
    ) +
    scale_y_continuous(
      limits = c(y_min, y_max),
      expand = c(0, 0)
    ) +
    coord_cartesian(clip = "off") +
    theme_void(base_family = base_family) +
    theme(
      plot.title = element_text(hjust = 0.5, size = 11, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, size = 9, color = "grey40"),
      legend.position = "none",
      plot.margin = margin(t = 10, r = 15, b = 15, l = 15)
    ) +
    labs(
      title = "Contig classification summary",
      subtitle = subtitle_text
    )

  girafe_obj <- ggiraph::girafe(ggobj = p_html, width_svg = 7.2, height_svg = plot_height)
  girafe_obj <- ggiraph::girafe_options(
    girafe_obj,
    opts_tooltip(css = paste0(
      "font-family: ", base_family, "; ",
      "font-size: 9pt; ",
      "background: white; padding: 6px; border: 1px solid #666;"
    )),
    opts_hover(css = "opacity: 0.8; stroke: black; stroke-width: 1;")
  )
  htmlwidgets::saveWidget(girafe_obj, out_html, selfcontained = TRUE)
  message("Interactive HTML saved to: ", out_html)
}
