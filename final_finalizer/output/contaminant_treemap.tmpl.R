#!/usr/bin/env Rscript
# Contamination treemap for final_finalizer
# Shows phylogenetic breakdown of detected contaminants as nested rectangles

if (!requireNamespace("pacman", quietly = TRUE)) {
  install.packages("pacman", repos = "https://cloud.r-project.org")
}
library(pacman)
pacman::p_load(
  readr, dplyr, stringr, ggplot2, tibble, tidyr,
  treemapify, scales, ggiraph, htmlwidgets, patchwork, showtext, sysfonts
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
contam_file <- "__CONTAMINANTS_TSV__"
out_pdf     <- "__OUTPDF__"
out_html    <- "__OUTHTML__"
plot_html   <- as.logical("__PLOTHTML__")
plot_suffix <- "__SUFFIX__"

base_family <- "Liberation Sans"
base_font_pt <- 8

# Read contaminant data
df <- read_tsv(contam_file, show_col_types = FALSE)

# Check if we have any data
if (nrow(df) == 0) {
  message("No contaminants detected. Skipping treemap generation.")
  quit(status = 0)
}

# Note: Coverage filtering already applied in Python before writing TSV
df_filtered <- df

# Check if we have actual taxonomy data or just genus from sci_name parsing
has_taxonomy <- any(!is.na(df_filtered$kingdom) & df_filtered$kingdom != "") ||
                any(!is.na(df_filtered$phylum) & df_filtered$phylum != "")
has_genus <- any(!is.na(df_filtered$genus) & df_filtered$genus != "")

# If genus column is empty, parse it from sci_name
if (!has_genus && "sci_name" %in% names(df_filtered)) {
  df_filtered <- df_filtered %>%
    mutate(
      genus = if_else(
        is.na(genus) | genus == "",
        str_extract(str_replace_all(sci_name, "_", " "), "^\\S+"),
        genus
      )
    )
  has_genus <- any(!is.na(df_filtered$genus) & df_filtered$genus != "")
}

# If kingdom is empty but phylum exists, infer kingdom from phylum
if (has_taxonomy && !any(!is.na(df_filtered$kingdom) & df_filtered$kingdom != "")) {
  df_filtered <- df_filtered %>%
    mutate(
      kingdom = case_when(
        !is.na(phylum) & phylum != "" & str_detect(phylum, "(?i)(archaeota|archaea)$") ~ "Archaea",
        !is.na(phylum) & phylum != "" & str_detect(phylum, "(?i)(ota|bacteria)$") ~ "Bacteria",
        !is.na(phylum) & phylum != "" ~ "Unknown",
        TRUE ~ kingdom
      )
    )
}

# Replace empty/NA values with "Unknown" for visualization
df_plot <- df_filtered %>%
  mutate(
    across(c(kingdom, family, genus, species), ~if_else(is.na(.) | . == "", "Unknown", .))
  )

# Aggregate by taxonomy path: Kingdom → Family → Genus → Species
df_agg <- df_plot %>%
  group_by(kingdom, family, genus, species) %>%
  summarise(
    total_bp = sum(length, na.rm = TRUE),
    n_contigs = n(),
    mean_depth = mean(depth_mean, na.rm = TRUE),
    mean_coverage = mean(coverage, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  # Create display label (species or genus if species is unknown)
  mutate(
    display_label = if_else(
      species == "Unknown" | species == genus,
      genus,
      paste0(genus, "\n", species)
    ),
    # Tooltip for interactive version
    tooltip = paste0(
      "Domain: ", kingdom,
      "\nFamily: ", family,
      "\nGenus: ", genus,
      "\nSpecies: ", species,
      "\n\nContigs: ", n_contigs,
      "\nTotal: ", label_bytes(units = "auto_binary")(total_bp),
      "\nMean coverage: ", sprintf("%.1f%%", mean_coverage * 100),
      if_else(!is.na(mean_depth), paste0("\nMean depth: ", sprintf("%.1fx", mean_depth)), "")
    )
  )

# Get unique families and create a distinct color for each
families <- unique(df_agg$family)
n_families <- length(families)

# Use a qualitative color palette - each family gets a distinct color
# This makes it easy to distinguish different contamination sources
if (n_families <= 8) {
  family_colors <- scales::brewer_pal(palette = "Set2")(n_families)
} else if (n_families <= 12) {
  family_colors <- scales::brewer_pal(palette = "Set3")(n_families)
} else {
  # For many families, use hue palette
  family_colors <- scales::hue_pal()(n_families)
}
names(family_colors) <- families

# Get domain summary for subtitle (kingdom column contains domain-level taxonomy)
domain_summary <- df_agg %>%
  group_by(kingdom) %>%
  summarise(bp = sum(total_bp), .groups = "drop") %>%
  arrange(desc(bp)) %>%
  mutate(pct = sprintf("%.0f%%", bp / sum(bp) * 100)) %>%
  mutate(label = paste0(kingdom, " (", pct, ")")) %>%
  pull(label) %>%
  paste(collapse = ", ")

# Calculate summary statistics
total_contigs <- sum(df_agg$n_contigs)
total_bp <- sum(df_agg$total_bp)

# Decide which plot to create based on available taxonomy data
if (has_taxonomy) {
  message("Creating treemap with full taxonomy")

  # Treemap with subgroups: Kingdom > Family > Genus
  p <- ggplot(df_agg, aes(
    area = total_bp,
    fill = family,
    label = display_label,
    subgroup = kingdom,
    subgroup2 = family
  )) +
    treemapify::geom_treemap(color = "white", size = 1) +
    treemapify::geom_treemap_subgroup_border(colour = "grey20", size = 2) +
    treemapify::geom_treemap_subgroup2_border(colour = "grey40", size = 1) +
    treemapify::geom_treemap_subgroup_text(
      place = "topleft",
      grow = FALSE,
      alpha = 0.7,
      colour = "black",
      fontface = "bold",
      size = 9,
      family = base_family,
      min.size = 2
    ) +
    treemapify::geom_treemap_text(
      colour = "white",
      place = "centre",
      size = 7,
      family = base_family,
      min.size = 2,
      grow = FALSE
    ) +
    scale_fill_manual(values = family_colors, name = "Family") +
    theme_minimal(base_family = base_family, base_size = base_font_pt) +
    theme(
      plot.title = element_text(hjust = 0.5, size = base_font_pt + 2, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, size = base_font_pt, color = "grey40"),
      legend.position = "bottom",
      legend.title = element_text(size = base_font_pt),
      legend.text = element_text(size = base_font_pt - 1)
    ) +
    guides(fill = guide_legend(nrow = 2, keywidth = unit(0.4, "cm"), keyheight = unit(0.4, "cm"))) +
    labs(
      title = paste0("Contaminant contigs"),
      subtitle = paste0(
        domain_summary, "\n",
        total_contigs, " contigs, ",
        label_bytes(units = "auto_binary")(total_bp)
      )
    )

} else {
  # Fallback: simple treemap by genus
  message("No full taxonomy available; creating treemap by genus")
  message("Tip: Install taxonkit and NCBI taxonomy database for full phylogenetic breakdown")

  # Aggregate by genus
  df_genus <- df_agg %>%
    group_by(genus) %>%
    summarise(
      total_bp = sum(total_bp),
      n_contigs = sum(n_contigs),
      mean_depth = weighted.mean(mean_depth, total_bp, na.rm = TRUE),
      mean_coverage = weighted.mean(mean_coverage, total_bp, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      tooltip = paste0(
        genus,
        "\nContigs: ", n_contigs,
        "\nTotal: ", label_bytes(units = "auto_binary")(total_bp),
        "\nMean coverage: ", sprintf("%.1f%%", mean_coverage * 100),
        if_else(!is.na(mean_depth), paste0("\nMean depth: ", sprintf("%.1fx", mean_depth)), "")
      )
    )

  n_genera <- nrow(df_genus)
  genus_colors <- scales::hue_pal()(min(n_genera, 20))
  names(genus_colors) <- df_genus$genus[order(-df_genus$total_bp)][1:min(n_genera, 20)]

  p <- ggplot(df_genus, aes(
    area = total_bp,
    fill = genus,
    label = paste0(genus, "\n", label_bytes(units = "auto_binary")(total_bp))
  )) +
    treemapify::geom_treemap(color = "white", size = 1) +
    treemapify::geom_treemap_text(
      colour = "white",
      place = "centre",
      size = 10,
      family = base_family,
      min.size = 3,
      grow = FALSE
    ) +
    scale_fill_manual(values = genus_colors, name = "Genus") +
    theme_minimal(base_family = base_family, base_size = base_font_pt) +
    theme(
      plot.title = element_text(hjust = 0.5, size = base_font_pt + 2, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, size = base_font_pt, color = "grey40"),
      plot.caption = element_text(size = base_font_pt - 1, color = "grey50", hjust = 0),
      legend.position = "none"  # Labels are on tiles
    ) +
    labs(
      title = paste0("Contamination by genus (", plot_suffix, ")"),
      subtitle = paste0(
        total_contigs, " contigs, ",
        label_bytes(units = "auto_binary")(total_bp)
      ),
      caption = "Note: Install taxonkit with NCBI taxonomy for hierarchical breakdown"
    )
}

# Create top 10 by depth table panel
# Get top 10 species by mean depth, with binomial names
top_by_depth <- df_agg %>%
  filter(!is.na(mean_depth) & mean_depth > 0) %>%
  arrange(desc(mean_depth)) %>%
  head(10) %>%
  mutate(
    # Clean species name - remove genus prefix if duplicated
    species_clean = str_replace_all(species, "_", " "),
    species_clean = if_else(
      str_starts(species_clean, paste0(genus, " ")),
      str_replace(species_clean, paste0("^", genus, " "), ""),
      species_clean
    ),
    # Create proper binomial name
    binomial = if_else(
      species_clean == "Unknown" | species_clean == genus | species_clean == "",
      paste0(genus, " sp."),
      paste0(genus, " ", species_clean)
    ),
    # Clean up any remaining underscores
    binomial = str_replace_all(binomial, "_", " "),
    # Rank for y-axis positioning
    rank = row_number()
  )

# Create the table panel
if (nrow(top_by_depth) > 0) {
  p_table <- ggplot(top_by_depth, aes(y = reorder(binomial, -rank))) +
    # Family color indicator (thin vertical line, height matches legend keys)
    geom_tile(aes(x = 0.02, fill = family), width = 0.015, height = 0.35) +
    # Binomial name (italic)
    geom_text(
      aes(x = 0.06, label = binomial),
      hjust = 0, size = (base_font_pt - 1) / .pt, fontface = "italic", family = base_family
    ) +
    # Depth value (right-aligned)
    geom_text(
      aes(x = 0.95, label = sprintf("%.1fx", mean_depth)),
      hjust = 1, size = (base_font_pt - 1) / .pt, family = base_family
    ) +
    scale_fill_manual(values = family_colors, guide = "none") +
    scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
    labs(title = "Top 10 by depth") +
    theme_void(base_family = base_family) +
    theme(
      plot.title = element_text(size = base_font_pt, face = "bold", hjust = 0.5, margin = margin(b = 5)),
      plot.margin = margin(t = 10, r = 5, b = 10, l = 5)
    )

  # Combine treemap and table using patchwork
  # Keep legend at bottom of treemap
  p_combined <- p + p_table +
    plot_layout(widths = c(2, 1))

  # Save combined PDF (max 7 inches wide)
  ggsave(out_pdf, p_combined, width = 7, height = 5, units = "in", dpi = 300)
} else {
  # No depth data - save treemap only
  ggsave(out_pdf, p, width = 7, height = 5, units = "in", dpi = 300)
}

message("Contamination treemap saved to: ", out_pdf)

# Note: Interactive HTML not supported for treemaps (treemapify lacks ggiraph integration)
# The static PDF provides the full visualization
if (plot_html) {
  message("Note: Interactive HTML not available for treemap visualization (PDF only)")
}
