#!/usr/bin/env Rscript
# Classification summary bar chart for final_finalizer
# Shows proportions of contigs across classification categories as a single horizontal stacked bar

if (!requireNamespace("pacman", quietly = TRUE)) {
  install.packages("pacman", repos = "https://cloud.r-project.org")
}
library(pacman)
pacman::p_load(
  readr, dplyr, stringr, ggplot2, tibble, tidyr, scales,
  ggiraph, htmlwidgets
)

# Placeholders - replaced by Python
summary_file <- "__SUMMARY__"
out_pdf      <- "__OUTPDF__"
out_html     <- "__OUTHTML__"
plot_html    <- as.logical("__PLOTHTML__")
plot_suffix  <- "__SUFFIX__"

base_family <- "Helvetica"
base_font_pt <- 10

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
# Order: chromosome assignments first, then special features, then problems
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

# Colorblind-safe palette with distinct colors for each category
classification_colors <- c(
  "chrom_assigned" = "#2166AC",      # Blue (primary success)
  "chrom_unassigned" = "#92C5DE",    # Light blue (partial success)
  "organelle_complete" = "#1B9E77",  # Teal (organelle)
  "organelle_debris" = "#B2DF8A",    # Light green (partial organelle)
  "rdna" = "#D95F02",                # Orange (rDNA)
  "contaminant" = "#E31A1C",         # Red (contaminant)
  "chrom_debris" = "#FDBF6F",        # Light orange (chromosome debris)
  "debris" = "#FB9A99",              # Light red (general debris)
  "unclassified" = "#CCCCCC"         # Gray (unknown)
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
    total_bp = sum(length, na.rm = TRUE),
    total_mb = sum(length_mb, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  filter(!is.na(classification))

# Calculate overall totals
total_contigs <- sum(df_agg$n_contigs)
total_mb <- sum(df_agg$total_mb)

# Calculate percentages
df_agg <- df_agg %>%
  mutate(
    pct_contigs = n_contigs / total_contigs * 100,
    pct_mb = total_mb / total_mb * 100,
    # Create position for cumulative bar
    pos_mb = cumsum(pct_mb) - pct_mb / 2
  )

# Create labels for below the bar
df_labels <- df_agg %>%
  mutate(
    label_text = paste0(
      classification_labels[as.character(classification)], "\n",
      n_contigs, " (", sprintf("%.1f%%", pct_contigs), ")\n",
      sprintf("%.1f Mb (%.1f%%)", total_mb, pct_mb)
    ),
    tooltip = paste0(
      classification_labels[as.character(classification)],
      "\nContigs: ", n_contigs, " (", sprintf("%.1f%%", pct_contigs), ")",
      "\nTotal: ", sprintf("%.1f Mb (%.1f%%)", total_mb, pct_mb)
    )
  )

# Static plot
p <- ggplot(df_agg, aes(x = pct_mb, y = 1, fill = classification)) +
  geom_col(width = 0.3, color = NA, position = position_stack(reverse = FALSE)) +
  # Labels below the bar
  geom_text(
    data = df_labels,
    aes(x = pos_mb, y = 0.5, label = label_text),
    size = base_font_pt / .pt * 0.85,
    family = base_family,
    lineheight = 0.85,
    vjust = 1,
    inherit.aes = FALSE
  ) +
  scale_fill_manual(
    values = classification_colors,
    labels = classification_labels,
    drop = FALSE
  ) +
  scale_x_continuous(
    expand = c(0, 0),
    labels = function(x) paste0(x, "%")
  ) +
  scale_y_continuous(expand = expansion(mult = c(0.3, 0.3))) +
  coord_cartesian(clip = "off") +
  theme_minimal(base_family = base_family, base_size = base_font_pt) +
  theme(
    plot.title = element_text(hjust = 0.5, size = base_font_pt + 2, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = base_font_pt, color = "grey40"),
    axis.title = element_blank(),
    axis.text.y = element_blank(),
    axis.text.x = element_text(size = base_font_pt),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_line(color = "grey90", linewidth = 0.3),
    legend.position = "none",
    plot.margin = margin(t = 10, r = 10, b = 80, l = 10)
  ) +
  labs(
    title = paste0("Contig classification summary (", plot_suffix, ")"),
    subtitle = paste0(total_contigs, " contigs, ", sprintf("%.1f Mb", total_mb))
  )

# Save PDF using cairo_pdf
ggsave(out_pdf, plot = p, width = 7.2, height = 4, units = "in", dpi = 300, device = cairo_pdf)
message("Classification summary bar plot saved to: ", out_pdf)

# Interactive HTML version
if (plot_html) {
  p_html <- ggplot(df_agg, aes(x = pct_mb, y = 1, fill = classification)) +
    ggiraph::geom_col_interactive(
      aes(tooltip = df_labels$tooltip),
      width = 0.3,
      color = NA,
      position = position_stack(reverse = FALSE)
    ) +
    # Labels below the bar (non-interactive)
    geom_text(
      data = df_labels,
      aes(x = pos_mb, y = 0.5, label = label_text),
      size = base_font_pt / .pt * 0.85,
      family = base_family,
      lineheight = 0.85,
      vjust = 1,
      inherit.aes = FALSE
    ) +
    scale_fill_manual(
      values = classification_colors,
      labels = classification_labels,
      drop = FALSE
    ) +
    scale_x_continuous(
      expand = c(0, 0),
      labels = function(x) paste0(x, "%")
    ) +
    scale_y_continuous(expand = expansion(mult = c(0.3, 0.3))) +
    coord_cartesian(clip = "off") +
    theme_minimal(base_family = base_family, base_size = base_font_pt) +
    theme(
      plot.title = element_text(hjust = 0.5, size = base_font_pt + 2, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, size = base_font_pt, color = "grey40"),
      axis.title = element_blank(),
      axis.text.y = element_blank(),
      axis.text.x = element_text(size = base_font_pt),
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_line(color = "grey90", linewidth = 0.3),
      legend.position = "none",
      plot.margin = margin(t = 10, r = 10, b = 80, l = 10)
    ) +
    labs(
      title = paste0("Contig classification summary (", plot_suffix, ")"),
      subtitle = paste0(total_contigs, " contigs, ", sprintf("%.1f Mb", total_mb))
    )

  girafe_obj <- ggiraph::girafe(ggobj = p_html, width_svg = 7.2, height_svg = 4)
  girafe_obj <- ggiraph::girafe_options(
    girafe_obj,
    opts_tooltip(css = paste0(
      "font-family: ", base_family, "; ",
      "font-size: ", base_font_pt, "pt; ",
      "background: white; padding: 6px; border: 1px solid #666;"
    )),
    opts_hover(css = "opacity: 0.8;")
  )
  htmlwidgets::saveWidget(girafe_obj, out_html, selfcontained = TRUE)
  message("Interactive HTML saved to: ", out_html)
}
