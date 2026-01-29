#!/usr/bin/env Rscript
# Classification summary contigs visualization for final_finalizer
# Shows individual contigs as shapes (pills for linear, rings for circular) faceted by classification

if (!requireNamespace("pacman", quietly = TRUE)) {
  install.packages("pacman", repos = "https://cloud.r-project.org")
}
library(pacman)
pacman::p_load(
  readr, dplyr, stringr, ggplot2, tibble, tidyr, scales,
  ggiraph, htmlwidgets, ggforce
)

# Placeholders - replaced by Python
summary_file <- "__SUMMARY__"
out_pdf      <- "__OUTPDF__"
out_html     <- "__OUTHTML__"
plot_html    <- as.logical("__PLOTHTML__")
plot_suffix  <- "__SUFFIX__"

base_family <- "Helvetica"
base_font_pt <- 8
facet_title_pt <- 9

# Read contig summary
df <- read_tsv(summary_file, show_col_types = FALSE,
               col_types = cols(
                 original_name = col_character(),
                 classification = col_character(),
                 length = col_double()
               ))

# Check if we have any data
if (nrow(df) == 0) {
  message("No contigs in summary. Skipping classification contigs plot.")
  quit(status = 0)
}

# Define classification order and colors (same as bar chart)
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

classification_colors <- c(
  "chrom_assigned" = "#2166AC",
  "chrom_unassigned" = "#92C5DE",
  "organelle_complete" = "#1B9E77",
  "organelle_debris" = "#B2DF8A",
  "rdna" = "#D95F02",
  "contaminant" = "#E31A1C",
  "chrom_debris" = "#FDBF6F",
  "debris" = "#FB9A99",
  "unclassified" = "#CCCCCC"
)

# Parse topology from original_name suffix (l = linear pill, c = circular ring)
df_plot <- df %>%
  mutate(
    classification = factor(classification, levels = classification_levels),
    length_mb = length / 1e6,
    # Extract topology: "l" = linear (pill), "c" = circular (ring)
    # Default to linear if no suffix
    topology = case_when(
      str_detect(original_name, "l$") ~ "linear",
      str_detect(original_name, "c$") ~ "circular",
      TRUE ~ "linear"
    ),
    tooltip = paste0(
      original_name, "\n",
      classification_labels[as.character(classification)], "\n",
      sprintf("%.3f Mb", length_mb)
    )
  ) %>%
  filter(!is.na(classification))

# Calculate stats per classification for facet labels
df_facet_stats <- df_plot %>%
  group_by(classification) %>%
  summarise(
    n_contigs = n(),
    total_mb = sum(length_mb, na.rm = TRUE),
    .groups = "drop"
  )

# Calculate overall totals for percentages
total_contigs <- sum(df_facet_stats$n_contigs)
total_mb <- sum(df_facet_stats$total_mb)

df_facet_stats <- df_facet_stats %>%
  mutate(
    pct_contigs = n_contigs / total_contigs * 100,
    pct_mb = total_mb / total_mb * 100,
    facet_label = paste0(
      classification_labels[as.character(classification)], "\n",
      n_contigs, " (", sprintf("%.1f%%", pct_contigs), ")  ",
      sprintf("%.1f Mb (%.1f%%)", total_mb, pct_mb)
    )
  )

# Arrange contigs in rows within each facet
# Sort by length descending, then pack into rows
df_layout <- df_plot %>%
  left_join(df_facet_stats %>% select(classification, facet_label), by = "classification") %>%
  group_by(classification) %>%
  arrange(desc(length_mb), .by_group = TRUE) %>%
  mutate(
    # Simple row packing: arrange bottom-to-top, aligned to left baseline
    row_num = row_number(),
    y_pos = row_num * 0.4  # vertical spacing between rows
  ) %>%
  ungroup()

# Calculate max x for each facet (for axis limits)
max_x <- max(df_layout$length_mb, na.rm = TRUE)
max_y_per_facet <- df_layout %>%
  group_by(classification) %>%
  summarise(max_y = max(y_pos), .groups = "drop")

# Pill/ring parameters
pill_height <- 0.25
ring_thickness <- 0.08  # thickness of circular ring

# Static plot
p <- ggplot(df_layout) +
  # Linear contigs as rounded rectangles (pills)
  ggforce::geom_ellipse(
    data = df_layout %>% filter(topology == "linear"),
    aes(x0 = length_mb / 2, y0 = y_pos,
        a = length_mb / 2, b = pill_height / 2,
        angle = 0,
        fill = classification),
    color = NA,
    alpha = 0.85
  ) +
  # Circular contigs as rings
  # Outer circle
  ggforce::geom_circle(
    data = df_layout %>% filter(topology == "circular"),
    aes(x0 = length_mb / 2, y0 = y_pos,
        r = pill_height / 2,
        fill = classification),
    color = NA,
    alpha = 0.85
  ) +
  # Inner circle (white to create ring effect)
  ggforce::geom_circle(
    data = df_layout %>% filter(topology == "circular"),
    aes(x0 = length_mb / 2, y0 = y_pos,
        r = pill_height / 2 - ring_thickness),
    fill = "white",
    color = NA
  ) +
  scale_fill_manual(
    values = classification_colors,
    labels = classification_labels,
    drop = FALSE
  ) +
  scale_x_continuous(
    expand = expansion(mult = c(0.02, 0.05)),
    labels = function(x) sprintf("%.1f", x)
  ) +
  facet_wrap(
    ~ facet_label,
    ncol = 1,
    scales = "free_y",
    strip.position = "top"
  ) +
  theme_minimal(base_family = base_family, base_size = base_font_pt) +
  theme(
    plot.title = element_text(hjust = 0.5, size = base_font_pt + 2, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = base_font_pt, color = "grey40"),
    strip.text = element_text(
      size = facet_title_pt,
      hjust = 0,
      face = "plain",
      family = base_family
    ),
    strip.background = element_rect(fill = "grey95", color = NA),
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.title.x = element_text(size = base_font_pt),
    axis.text.x = element_text(size = base_font_pt - 1),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_line(color = "grey90", linewidth = 0.3),
    panel.spacing.y = unit(0.3, "lines"),
    legend.position = "none",
    plot.margin = margin(t = 10, r = 10, b = 10, l = 10)
  ) +
  labs(
    title = paste0("Contig classification (", plot_suffix, ")"),
    subtitle = paste0(total_contigs, " contigs, ", sprintf("%.1f Mb", total_mb)),
    x = "Length (Mb)"
  )

# Save PDF using cairo_pdf with max 7.2" width
ggsave(out_pdf, plot = p, width = 7.2, height = 10, units = "in", dpi = 300, device = cairo_pdf)
message("Classification contigs plot saved to: ", out_pdf)

# Interactive HTML version
# Use simplified geoms for compatibility with ggiraph
if (plot_html) {
  p_html <- ggplot(df_layout) +
    # Linear contigs as rectangles with rounded ends (approximated with rect)
    ggiraph::geom_rect_interactive(
      data = df_layout %>% filter(topology == "linear"),
      aes(xmin = 0, xmax = length_mb,
          ymin = y_pos - pill_height / 2, ymax = y_pos + pill_height / 2,
          fill = classification,
          tooltip = tooltip),
      color = NA,
      alpha = 0.85
    ) +
    # Circular contigs as filled points
    ggiraph::geom_point_interactive(
      data = df_layout %>% filter(topology == "circular"),
      aes(x = length_mb / 2, y = y_pos,
          fill = classification,
          tooltip = tooltip),
      shape = 21,  # filled circle with border
      size = 4,
      color = "white",
      stroke = 1.5,
      alpha = 0.85
    ) +
    scale_fill_manual(
      values = classification_colors,
      labels = classification_labels,
      drop = FALSE
    ) +
    scale_x_continuous(
      expand = expansion(mult = c(0.02, 0.05)),
      labels = function(x) sprintf("%.1f", x)
    ) +
    facet_wrap(
      ~ facet_label,
      ncol = 1,
      scales = "free_y",
      strip.position = "top"
    ) +
    theme_minimal(base_family = base_family, base_size = base_font_pt) +
    theme(
      plot.title = element_text(hjust = 0.5, size = base_font_pt + 2, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, size = base_font_pt, color = "grey40"),
      strip.text = element_text(
        size = facet_title_pt,
        hjust = 0,
        face = "plain",
        family = base_family
      ),
      strip.background = element_rect(fill = "grey95", color = NA),
      axis.title.y = element_blank(),
      axis.text.y = element_blank(),
      axis.title.x = element_text(size = base_font_pt),
      axis.text.x = element_text(size = base_font_pt - 1),
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_line(color = "grey90", linewidth = 0.3),
      panel.spacing.y = unit(0.3, "lines"),
      legend.position = "none",
      plot.margin = margin(t = 10, r = 10, b = 10, l = 10)
    ) +
    labs(
      title = paste0("Contig classification (", plot_suffix, ")"),
      subtitle = paste0(total_contigs, " contigs, ", sprintf("%.1f Mb", total_mb)),
      x = "Length (Mb)"
    )

  girafe_obj <- ggiraph::girafe(ggobj = p_html, width_svg = 7.2, height_svg = 10)
  girafe_obj <- ggiraph::girafe_options(
    girafe_obj,
    opts_tooltip(css = paste0(
      "font-family: ", base_family, "; ",
      "font-size: ", base_font_pt, "pt; ",
      "background: white; padding: 6px; border: 1px solid #666;"
    )),
    opts_hover(css = "opacity: 1.0; stroke: black; stroke-width: 1;")
  )
  htmlwidgets::saveWidget(girafe_obj, out_html, selfcontained = TRUE)
  message("Interactive HTML saved to: ", out_html)
}
