#!/usr/bin/env Rscript

if (!requireNamespace("pacman", quietly = TRUE)) {
  install.packages("pacman", repos = "https://cloud.r-project.org")
}
library(pacman)
pacman::p_load(
  readr, dplyr, stringr, ggplot2, tibble, tidyr, patchwork,
  grid, ggiraph, htmlwidgets, ggokabeito, colorspace, showtext, sysfonts
)

# Register Liberation Sans for consistent cross-platform rendering
sysfonts::font_add("Liberation Sans",
  regular    = "/usr/share/fonts/liberation-sans/LiberationSans-Regular.ttf",
  bold       = "/usr/share/fonts/liberation-sans/LiberationSans-Bold.ttf",
  italic     = "/usr/share/fonts/liberation-sans/LiberationSans-Italic.ttf",
  bolditalic = "/usr/share/fonts/liberation-sans/LiberationSans-BoldItalic.ttf"
)
showtext_auto()

summary_file <- "__SUMMARY__"
out_pdf      <- "__OUTPDF__"
out_html     <- "__OUTHTML__"
plot_html    <- as.logical("__PLOTHTML__")
plot_title_suffix <- "__SUFFIX__"
base_family <- "Liberation Sans"
base_font_pt <- 8
axis_title_margin_pt <- 4
axis_text_margin_pt <- 4

axis_theme <- theme(
  axis.title.x = element_text(size = base_font_pt, family = base_family, margin = margin(t = axis_title_margin_pt)),
  axis.title.y = element_text(size = base_font_pt, family = base_family, margin = margin(r = axis_title_margin_pt)),
  axis.text.x = element_text(size = base_font_pt, family = base_family, margin = margin(t = axis_text_margin_pt)),
  axis.text.y = element_text(size = base_font_pt, family = base_family)
)

apply_legend_theme <- function(p, text_pt = base_font_pt, key_pt = base_font_pt, tight = TRUE) {
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

# Read summary data
df <- read_tsv(summary_file, show_col_types = FALSE)

# Check if depth data is available
has_depth <- "depth_mean" %in% names(df) && any(!is.na(df$depth_mean))

if (!has_depth) {
  message("No depth data available in summary file. Skipping depth plot generation.")
  quit(status = 0)
}

# Classification order for consistent display
# Each classification is followed by its corresponding debris category
# organelle_debris is split into chrC_debris and chrM_debris when assigned_ref_id is available
classification_order <- c(
  "chrom_assigned",
  "chrom_debris",
  "chrom_unassigned",
  "chrC",
  "chrC_debris",
  "chrM",
  "chrM_debris",
  "organelle_complete",
  "organelle_debris",
  "rDNA",
  "contaminant",
  "debris",
  "unclassified"
)

# Okabe-Ito colorblind-friendly palette
# Each debris category uses a lightened shade of its parent classification
oi <- palette_okabe_ito()
classification_colors <- c(
  "chrom_assigned"     = oi[5],                      # blue
  "chrom_debris"       = lighten(oi[5], amount = 0.5),
  "chrom_unassigned"   = oi[2],                      # sky blue
  "chrC"               = oi[3],                      # bluish green
  "chrC_debris"        = lighten(oi[3], amount = 0.5),
  "chrM"               = oi[6],                      # vermillion
  "chrM_debris"        = lighten(oi[6], amount = 0.5),
  "organelle_complete" = oi[3],                      # bluish green (same as chrC)
  "organelle_debris"   = oi[4],                      # yellow
  "rDNA"               = oi[1],                      # orange
  "contaminant"        = oi[7],                      # reddish purple
  "debris"             = oi[9],                      # gray
  "unclassified"       = oi[8]                       # dark gray (black)
)

# Prepare data for plotting
# Separate chrC and chrM organelles into their own categories
# Also split organelle_debris into chrC_debris and chrM_debris based on assigned_ref_id
df_depth <- df %>%
  filter(!is.na(depth_mean)) %>%
  mutate(
    # Create display classification that separates chrC and chrM
    # For organelle_debris, use assigned_ref_id to determine chrC_debris vs chrM_debris
    display_class = case_when(
      contig == "chrC" ~ "chrC",
      contig == "chrM" ~ "chrM",
      classification == "organelle_debris" & assigned_ref_id == "chrC" ~ "chrC_debris",
      classification == "organelle_debris" & assigned_ref_id == "chrM" ~ "chrM_debris",
      TRUE ~ classification
    ),
    display_class = factor(display_class, levels = classification_order),
    length_mb = as.numeric(length) / 1e6,
    # For chrom_assigned, create chrom_order for sorting
    chrom_num = as.integer(str_extract(assigned_chrom_id, "\\d+")),
    chrom_num = if_else(is.na(chrom_num), 999L, chrom_num)
  ) %>%
  arrange(display_class, chrom_num, assigned_chrom_id, desc(length_mb))

# ----------------------------
# Top: Mean depth by classification (box + violin)
# Organelles (chrC and chrM) are plotted separately
# ----------------------------
p_class_depth <- ggplot(df_depth, aes(x = display_class, y = depth_mean, fill = display_class)) +
  geom_violin(alpha = 0.5, scale = "width", width = 0.7) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA, width = 0.3) +
  geom_jitter(aes(color = display_class), width = 0.12, alpha = 0.5, size = 1.2, show.legend = FALSE) +
  scale_fill_manual(values = classification_colors, drop = FALSE, guide = "none") +
  scale_color_manual(values = classification_colors, drop = FALSE) +
  scale_x_discrete(
    labels = function(x) str_replace_all(x, "_", "\n")
  ) +
  theme_classic(base_family = base_family, base_size = base_font_pt) +
  axis_theme +
  theme(
    plot.title = element_text(hjust = 0.5, size = base_font_pt, family = base_family),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = base_font_pt - 1),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.3)
  ) +
  labs(
    title = paste0("Read depth by classification (", plot_title_suffix, ")"),
    x = "Classification",
    y = "Mean read depth (x)"
  )

# ----------------------------
# Middle: Depth for chromosome-assigned contigs, ordered by chromosome
# Includes chrC and chrM organelles if present
# ----------------------------
df_chrom <- df_depth %>%
  filter(
    classification == "chrom_assigned" |
    contig == "chrC" |
    contig == "chrM"
  ) %>%
  mutate(
    # Use contig name for organelles, assigned_chrom_id for others
    chrom_id = case_when(
      contig == "chrC" ~ "chrC",
      contig == "chrM" ~ "chrM",
      TRUE ~ assigned_chrom_id
    ),
    # Extract numeric part for sorting (organelles get special values)
    chrom_sort = case_when(
      contig == "chrC" ~ 9998L,
      contig == "chrM" ~ 9999L,
      TRUE ~ chrom_num
    ),
    # Label for display
    # Note: assigned_subgenome may be NA (R's missing) or "NA" (string)
    # Check for both cases to avoid labels like "chr1NA"
    chrom_label = if_else(
      classification == "chrom_assigned" &
        !is.na(assigned_subgenome) &
        assigned_subgenome != "" &
        assigned_subgenome != "NA",
      paste0(chrom_id, assigned_subgenome),
      chrom_id
    )
  ) %>%
  arrange(chrom_sort, chrom_label, desc(length_mb)) %>%
  mutate(plot_order = row_number())

if (nrow(df_chrom) > 0) {
  # Create x-axis breaks at chromosome boundaries
  chrom_breaks <- df_chrom %>%
    group_by(chrom_label) %>%
    summarise(mid_pos = mean(plot_order), .groups = "drop")

  # Color by subgenome for regular chromosomes, distinct colors for organelles
  df_chrom <- df_chrom %>%
    mutate(
      fill_category = case_when(
        contig == "chrC" ~ "chrC",
        contig == "chrM" ~ "chrM",
        !is.na(assigned_subgenome) ~ as.character(assigned_subgenome),
        TRUE ~ "NA"
      )
    )

  p_chrom_depth <- ggplot(df_chrom, aes(x = plot_order, y = depth_mean)) +
    geom_col(aes(fill = fill_category), width = 0.8, alpha = 0.85) +
    scale_fill_manual(
      values = c("A" = "#1F77B4", "B" = "#FF7F0E", "NA" = "#7F7F7F",
                 "At" = "#1F77B4", "Dt" = "#FF7F0E", "G" = "#2CA02C",
                 "chrC" = "#2CA02C", "chrM" = "#D62728"),
      name = "Subgenome"
    ) +
    scale_x_continuous(
      breaks = chrom_breaks$mid_pos,
      labels = str_replace(chrom_breaks$chrom_label, "^chr", ""),
      expand = expansion(mult = c(0.01, 0.01))
    ) +
    theme_classic(base_family = base_family, base_size = base_font_pt) +
    axis_theme +
    theme(
      plot.title = element_text(hjust = 0.5, size = base_font_pt, family = base_family),
      axis.text.x = element_text(size = base_font_pt - 1),
      panel.grid.major.y = element_line(color = "grey90", linewidth = 0.3),
      legend.position = "top"
    ) +
    labs(
      title = "Chromosome contig depth",
      x = "Chromosome",
      y = "Mean depth (x)"
    )

  p_chrom_depth <- apply_legend_theme(p_chrom_depth, text_pt = 6, key_pt = 6, tight = TRUE) +
    theme(legend.margin = margin(1, 1, 1, 1))
} else {
  p_chrom_depth <- ggplot() +
    theme_void(base_family = base_family, base_size = base_font_pt) +
    labs(title = "Chromosome depth (no data)")
}

# ----------------------------
# Bottom: Breadth coverage (1x and 10x) by classification
# Organelles (chrC and chrM) are plotted separately
# ----------------------------
df_breadth <- df_depth %>%
  select(contig, classification, display_class, depth_breadth_1x, depth_breadth_10x) %>%
  tidyr::pivot_longer(
    cols = c(depth_breadth_1x, depth_breadth_10x),
    names_to = "metric",
    values_to = "breadth"
  ) %>%
  mutate(
    metric = case_when(
      metric == "depth_breadth_1x" ~ ">=1x",
      metric == "depth_breadth_10x" ~ ">=10x",
      TRUE ~ metric
    ),
    metric = factor(metric, levels = c(">=1x", ">=10x"))
  )

p_breadth <- ggplot(df_breadth, aes(x = display_class, y = breadth, fill = metric)) +
  geom_boxplot(alpha = 0.7, position = position_dodge(width = 0.75), width = 0.6, outlier.size = 0.8) +
  scale_fill_manual(values = c(">=1x" = "#66C2A5", ">=10x" = "#1F77B4"), name = "Coverage") +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2), labels = scales::percent) +
  scale_x_discrete(
    labels = function(x) str_replace_all(x, "_", "\n")
  ) +
  theme_classic(base_family = base_family, base_size = base_font_pt) +
  axis_theme +
  theme(
    plot.title = element_text(hjust = 0.5, size = base_font_pt, family = base_family),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = base_font_pt - 1),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.3),
    legend.position = "top"
  ) +
  labs(
    title = "Coverage breadth by classification",
    x = "Classification",
    y = "Fraction of bases covered"
  )

p_breadth <- apply_legend_theme(p_breadth, text_pt = 6, key_pt = 6, tight = TRUE) +
  theme(legend.margin = margin(1, 1, 1, 1))

# ----------------------------
# Combine plots
# ----------------------------
full_plot <- p_class_depth / p_chrom_depth / p_breadth +
  plot_layout(heights = c(1, 1, 1))

ggsave(out_pdf, plot = full_plot, width = 8, height = 10, units = "in", dpi = 300)
message("Depth overview plot saved to: ", out_pdf)

# ----------------------------
# Interactive HTML version
# ----------------------------
if (plot_html) {
  # Top plot with tooltips
  df_depth_tooltip <- df_depth %>%
    mutate(
      tooltip = paste0(
        "Contig: ", contig,
        "\nClassification: ", display_class,
        "\nLength: ", sprintf("%.2f", length_mb), " Mbp",
        "\nMean depth: ", sprintf("%.1f", depth_mean), "x",
        "\nMedian depth: ", sprintf("%.1f", depth_median), "x",
        "\nBreadth 1x: ", sprintf("%.1f%%", depth_breadth_1x * 100),
        "\nBreadth 10x: ", sprintf("%.1f%%", depth_breadth_10x * 100)
      )
    )

  p_class_depth_html <- ggplot(df_depth_tooltip, aes(x = display_class, y = depth_mean, fill = display_class)) +
    geom_violin(alpha = 0.5, scale = "width", width = 0.7) +
    geom_boxplot(alpha = 0.7, outlier.shape = NA, width = 0.3) +
    ggiraph::geom_point_interactive(
      aes(color = display_class, tooltip = tooltip),
      position = position_jitter(width = 0.12),
      alpha = 0.5, size = 1.2, show.legend = FALSE
    ) +
    scale_fill_manual(values = classification_colors, drop = FALSE, guide = "none") +
    scale_color_manual(values = classification_colors, drop = FALSE) +
    scale_x_discrete(
      labels = function(x) str_replace_all(x, "_", "\n")
    ) +
    theme_classic(base_family = base_family, base_size = base_font_pt) +
    axis_theme +
    theme(
      plot.title = element_text(hjust = 0.5, size = base_font_pt, family = base_family),
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = base_font_pt - 1),
      panel.grid.major.y = element_line(color = "grey90", linewidth = 0.3)
    ) +
    labs(
      title = paste0("Read depth by classification (", plot_title_suffix, ")"),
      x = "Classification",
      y = "Mean read depth (x)"
    )

  # Chromosome depth with tooltips
  if (nrow(df_chrom) > 0) {
    df_chrom_tooltip <- df_chrom %>%
      mutate(
        tooltip = paste0(
          "Contig: ", contig,
          "\nChromosome: ", chrom_label,
          "\nLength: ", sprintf("%.2f", length_mb), " Mbp",
          "\nMean depth: ", sprintf("%.1f", depth_mean), "x"
        )
      )

    p_chrom_depth_html <- ggplot(df_chrom_tooltip, aes(x = plot_order, y = depth_mean)) +
      ggiraph::geom_col_interactive(
        aes(fill = fill_category, tooltip = tooltip),
        width = 0.8, alpha = 0.85
      ) +
      scale_fill_manual(
        values = c("A" = "#1F77B4", "B" = "#FF7F0E", "NA" = "#7F7F7F",
                   "At" = "#1F77B4", "Dt" = "#FF7F0E", "G" = "#2CA02C",
                   "chrC" = "#2CA02C", "chrM" = "#D62728"),
        name = "Subgenome"
      ) +
      scale_x_continuous(
        breaks = chrom_breaks$mid_pos,
        labels = str_replace(chrom_breaks$chrom_label, "^chr", ""),
        expand = expansion(mult = c(0.01, 0.01))
      ) +
      theme_classic(base_family = base_family, base_size = base_font_pt) +
      axis_theme +
      theme(
        plot.title = element_text(hjust = 0.5, size = base_font_pt, family = base_family),
        axis.text.x = element_text(size = base_font_pt - 1),
        panel.grid.major.y = element_line(color = "grey90", linewidth = 0.3),
        legend.position = "top"
      ) +
      labs(
        title = "Chromosome contig depth",
        x = "Chromosome",
        y = "Mean depth (x)"
      )

    p_chrom_depth_html <- apply_legend_theme(p_chrom_depth_html, text_pt = 6, key_pt = 6, tight = TRUE) +
      theme(legend.margin = margin(1, 1, 1, 1))
  } else {
    p_chrom_depth_html <- p_chrom_depth
  }

  full_plot_html <- p_class_depth_html / p_chrom_depth_html / p_breadth +
    plot_layout(heights = c(1, 1, 1))

  girafe_obj <- ggiraph::girafe(ggobj = full_plot_html, width_svg = 8, height_svg = 10)
  girafe_obj <- ggiraph::girafe_options(
    girafe_obj,
    opts_tooltip(css = paste0(
      "font-family: ", base_family, "; ",
      "font-size: ", base_font_pt, "pt; ",
      "background: white; padding: 4px; border: 1px solid #666;"
    )),
    opts_hover(css = paste0("font-family: ", base_family, ";"))
  )
  htmlwidgets::saveWidget(girafe_obj, out_html, selfcontained = TRUE)
  message("Interactive depth plot saved to: ", out_html)
}
