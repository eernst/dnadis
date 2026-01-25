#!/usr/bin/env Rscript
# Contamination summary table using gt
# Shows top contaminants ranked by abundance (depth × length)
# Uses CSS gradient bars with text overlay for Total Mb and Depth columns

if (!requireNamespace("pacman", quietly = TRUE)) {
  install.packages("pacman", repos = "https://cloud.r-project.org")
}
library(pacman)
pacman::p_load(readr, dplyr, stringr, tidyr, gt, scales)

# Placeholders - replaced by Python
contam_file <- "__CONTAMINANTS_TSV__"
out_html    <- "__OUTHTML__"
plot_suffix <- "__SUFFIX__"
top_n       <- as.integer("__TOP_N__")

# Note: PDF output is not supported for this table because it uses HTML/CSS
# for the inline bar visualizations. HTML output works in any browser and
# can be printed to PDF if needed.

# Read contaminant data
df <- read_tsv(contam_file, show_col_types = FALSE)

# Check if we have any data
if (nrow(df) == 0) {
  message("No contaminants detected. Skipping table generation.")
  quit(status = 0)
}

# Check for depth data
has_depth <- "depth_mean" %in% names(df) && !all(is.na(df$depth_mean))

if (has_depth) {
  # Filter to contigs with non-NA depth (include depth=0)
  df <- df %>%
    filter(!is.na(depth_mean))
}

if (nrow(df) == 0) {
  message("No contaminants with depth data. Skipping table.")
  quit(status = 0)
}

# Determine contig topology from name suffix and create binomial names
df <- df %>%
  mutate(
    is_circular = str_detect(contig, "c$"),
    # Clean species name - remove genus prefix if duplicated
    species_clean = str_replace_all(species, "_", " "),
    species_clean = if_else(
      !is.na(species_clean) & str_starts(species_clean, paste0(genus, " ")),
      str_replace(species_clean, paste0("^", genus, " "), ""),
      species_clean
    ),
    # Create binomial name
    binomial = case_when(
      is.na(species_clean) | species_clean == "" | species_clean == genus | species_clean == "Unknown" ~
        paste0(genus, " sp."),
      TRUE ~ paste0(genus, " ", species_clean)
    ),
    binomial = str_replace_all(binomial, "_", " "),
    # Abundance metric: depth × length if available, otherwise just length
    abundance = if (has_depth) depth_mean * length else length
  )

# Aggregate by binomial with min/max for spread displays
df_agg <- df %>%
  group_by(binomial, kingdom, family) %>%
  summarise(
    n_contigs = n(),
    total_bp = sum(length),
    abundance = sum(abundance),
    depth_mean_agg = if (has_depth) mean(depth_mean) else NA_real_,
    depth_min = if (has_depth) min(depth_mean) else NA_real_,
    depth_max = if (has_depth) max(depth_mean) else NA_real_,
    length_min = min(length),
    length_max = max(length),
    n_circular = sum(is_circular),
    mean_coverage = mean(coverage),
    .groups = "drop"
  ) %>%
  arrange(desc(abundance)) %>%
  head(top_n) %>%
  mutate(rank = row_number())

if (nrow(df_agg) == 0) {
  message("No contaminants to display. Skipping table.")
  quit(status = 0)
}

# Domain abbreviations and colors
domain_info <- tibble(
  kingdom = c("Bacteria", "Archaea", "Viruses", "Animalia", "Plantae", "Fungi", "Chromista", "Protozoa"),
  domain_abbr = c("Bac", "Arc", "Vir", "Ani", "Pla", "Fun", "Chr", "Pro"),
  domain_color = c("#E41A1C", "#984EA3", "#4DAF4A", "#377EB8", "#4DAF4A", "#FF7F00", "#A65628", "#999999")
)

df_agg <- df_agg %>%
  left_join(domain_info, by = "kingdom") %>%
  mutate(
    domain_abbr = if_else(is.na(domain_abbr), "Unk", domain_abbr),
    domain_color = if_else(is.na(domain_color), "#999999", domain_color)
  )

# Family color palette (consistent with other plots)
families <- unique(df_agg$family)
n_families <- length(families)
if (n_families <= 8) {
  family_colors <- scales::brewer_pal(palette = "Set2")(max(n_families, 3))[1:n_families]
} else if (n_families <= 12) {
  family_colors <- scales::brewer_pal(palette = "Set3")(max(n_families, 3))[1:n_families]
} else {
  family_colors <- scales::hue_pal()(n_families)
}
names(family_colors) <- families

df_agg <- df_agg %>%
  mutate(
    family_color = family_colors[family],
    # Circular indicator
    circ_indicator = case_when(
      n_circular == 0 ~ "\u25CB",
      n_circular == n_contigs ~ "\u25CF",
      TRUE ~ "\u25D0"  # Mixed
    ),
    total_mb = total_bp / 1e6
  )

# Calculate bar percentages (relative to max in each column)
max_total <- max(df_agg$total_mb)

df_agg <- df_agg %>%
  mutate(
    total_pct = total_mb / max_total * 100,
    # Text displays with spread for multi-contig entries
    # Total: shows aggregate Mb plus individual contig size range
    total_text = if_else(
      n_contigs == 1,
      sprintf("%.2f Mb", total_mb),
      sprintf("%.2f Mb (%.2f\u2013%.2f)", total_mb, length_min / 1e6, length_max / 1e6)
    )
  )

# Add depth columns only if depth data available
if (has_depth) {
  max_depth <- max(df_agg$depth_mean_agg)
  df_agg <- df_agg %>%
    mutate(
      depth_pct = depth_mean_agg / max_depth * 100,
      depth_text = if_else(
        n_contigs == 1,
        sprintf("%.0f\u00D7", depth_mean_agg),
        sprintf("%.0f\u00D7 (%.0f\u2013%.0f)", depth_mean_agg, depth_min, depth_max)
      )
    )
}

# Create table data - conditionally include depth columns
if (has_depth) {
  tbl_data <- df_agg %>%
    select(
      rank, domain_abbr, domain_color, family, family_color, binomial,
      n_contigs, total_text, total_pct, depth_text, depth_pct,
      mean_coverage, circ_indicator
    )
} else {
  tbl_data <- df_agg %>%
    select(
      rank, domain_abbr, domain_color, family, family_color, binomial,
      n_contigs, total_text, total_pct,
      mean_coverage, circ_indicator
    )
}

# Build gt table - base structure
gt_tbl <- tbl_data %>%
  gt() %>%
  # Format coverage as percentage
  fmt_percent(columns = mean_coverage, decimals = 0) %>%
  # Style domain with background color badge
  text_transform(
    locations = cells_body(columns = domain_abbr),
    fn = function(x) {
      colors <- tbl_data$domain_color
      paste0(
        "<span style='background-color:", colors,
        "; color: white; padding: 2px 4px; border-radius: 3px; font-size: 10px;'>",
        x, "</span>"
      )
    }
  ) %>%
  # Style family with colored left border
  text_transform(
    locations = cells_body(columns = family),
    fn = function(x) {
      colors <- tbl_data$family_color
      paste0("<span style='border-left: 4px solid ", colors, "; padding-left: 6px;'>", x, "</span>")
    }
  ) %>%
  # Style binomial as italic
  text_transform(
    locations = cells_body(columns = binomial),
    fn = function(x) paste0("<em>", x, "</em>")
  ) %>%
  # Total Mb with CSS gradient bar background
  text_transform(
    locations = cells_body(columns = total_text),
    fn = function(x) {
      pcts <- tbl_data$total_pct
      paste0(
        "<div style='background: linear-gradient(to right, rgba(70,130,180,0.3) ", pcts,
        "%, transparent ", pcts, "%); padding: 2px 4px; border-radius: 2px; white-space: nowrap;'>",
        x, "</div>"
      )
    }
  ) %>%
  # Column alignment
  cols_align(align = "center", columns = c(rank, n_contigs, circ_indicator, mean_coverage)) %>%
  cols_align(align = "left", columns = total_text) %>%
  # Table styling
  tab_style(style = cell_text(size = px(11)), locations = cells_body()) %>%
  tab_style(style = cell_text(size = px(11), weight = "bold"), locations = cells_column_labels()) %>%
  tab_options(
    data_row.padding = px(4),
    column_labels.padding = px(6),
    table.font.size = px(11)
  ) %>%
  # Footnotes for indicators
  tab_footnote(
    footnote = "\u25CF = all circularized, \u25D0 = mixed, \u25CB = all linear",
    locations = cells_column_labels(columns = circ_indicator)
  ) %>%
  # Source note with family legend
  tab_source_note(
    source_note = md(paste0(
      "**Families:** ",
      paste(sapply(names(family_colors), function(f) {
        paste0("<span style='border-left: 4px solid ", family_colors[f],
               "; padding-left: 4px; margin-right: 12px;'>", f, "</span>")
      }), collapse = " ")
    ))
  )

# Add depth-specific formatting if depth data available
if (has_depth) {
  gt_tbl <- gt_tbl %>%
    cols_label(
      rank = "", domain_abbr = "Domain", family = "Family", binomial = "Species",
      n_contigs = "n", total_text = "Total Mb", depth_text = "Depth",
      mean_coverage = "Cov%", circ_indicator = "Circ"
    ) %>%
    cols_hide(c(domain_color, family_color, total_pct, depth_pct)) %>%
    text_transform(
      locations = cells_body(columns = depth_text),
      fn = function(x) {
        pcts <- tbl_data$depth_pct
        paste0(
          "<div style='background: linear-gradient(to right, rgba(128,128,128,0.3) ", pcts,
          "%, transparent ", pcts, "%); padding: 2px 4px; border-radius: 2px; white-space: nowrap;'>",
          x, "</div>"
        )
      }
    ) %>%
    cols_align(align = "left", columns = depth_text) %>%
    cols_width(
      rank ~ px(30), domain_abbr ~ px(50), family ~ px(120), binomial ~ px(180),
      n_contigs ~ px(30), total_text ~ px(150), depth_text ~ px(130),
      mean_coverage ~ px(45), circ_indicator ~ px(35)
    ) %>%
    tab_header(
      title = md(paste0("**Top ", nrow(tbl_data), " Contaminants**")),
      subtitle = md(paste0("Ranked by abundance (depth \u00D7 length) \u2014 ", plot_suffix))
    )
} else {
  gt_tbl <- gt_tbl %>%
    cols_label(
      rank = "", domain_abbr = "Domain", family = "Family", binomial = "Species",
      n_contigs = "n", total_text = "Total Mb",
      mean_coverage = "Cov%", circ_indicator = "Circ"
    ) %>%
    cols_hide(c(domain_color, family_color, total_pct)) %>%
    cols_width(
      rank ~ px(30), domain_abbr ~ px(50), family ~ px(120), binomial ~ px(180),
      n_contigs ~ px(30), total_text ~ px(150),
      mean_coverage ~ px(45), circ_indicator ~ px(35)
    ) %>%
    tab_header(
      title = md(paste0("**Top ", nrow(tbl_data), " Contaminants**")),
      subtitle = md(paste0("Ranked by total length \u2014 ", plot_suffix))
    )
}

# Save as HTML
gtsave(gt_tbl, out_html)
message("Contaminant table saved to: ", out_html)
