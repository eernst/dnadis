#!/usr/bin/env Rscript
# Contamination summary table using gt and gtExtras
# Shows top contaminants ranked by abundance (depth × length)

if (!requireNamespace("pacman", quietly = TRUE)) {
  install.packages("pacman", repos = "https://cloud.r-project.org")
}
library(pacman)
pacman::p_load(
  readr, dplyr, stringr, tidyr, gt, gtExtras, scales
)

# Placeholders - replaced by Python
contam_file <- "__CONTAMINANTS_TSV__"
out_html    <- "__OUTHTML__"
plot_suffix <- "__SUFFIX__"
top_n       <- as.integer("__TOP_N__")

# Note: PDF output is not supported for this table because gtExtras sparklines
# are SVG-based and require a browser engine (Chrome/Chromium) to render.
# HTML output works in any browser and can be printed to PDF if needed.

# Read contaminant data
df <- read_tsv(contam_file, show_col_types = FALSE)

# Check if we have any data
if (nrow(df) == 0) {
  message("No contaminants detected. Skipping table generation.")
  quit(status = 0)
}

# Check for depth data
if (!"depth_mean" %in% names(df) || all(is.na(df$depth_mean))) {
  message("No depth data available. Skipping contaminant table (requires --reads).")
  quit(status = 0)
}

# Filter to contigs with non-NA depth (include depth=0)
df <- df %>%
  filter(!is.na(depth_mean))

if (nrow(df) == 0) {
  message("No contaminants with depth data. Skipping table.")
  quit(status = 0)
}

# Determine contig topology from name suffix
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
    # Abundance metric: depth × length (proxy for sequencing effort)
    abundance = depth_mean * length
  )

# Aggregate by binomial
df_agg <- df %>%
  group_by(binomial, kingdom, family) %>%
  summarise(
    n_contigs = n(),
    total_bp = sum(length),
    abundance = sum(abundance),
    depths = list(depth_mean),
    lengths = list(length / 1e6),  # Convert to Mb for plotting
    coverages = list(coverage),
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

# Add family color to data
df_agg <- df_agg %>%
  mutate(family_color = family_colors[family])

# Circular indicator
df_agg <- df_agg %>%
  mutate(
    circ_indicator = case_when(
      n_circular == 0 ~ "\u25CB",
      n_circular == n_contigs ~ "\u25CF",
      TRUE ~ "\u25D0"  # Mixed
    )
  )

# Format total_bp for display (always Mb)
df_agg <- df_agg %>%
  mutate(
    total_mb = sprintf("%.2f", total_bp / 1e6)
  )

# Create the base table data
tbl_data <- df_agg %>%
  select(
    rank,
    domain_abbr,
    domain_color,
    family,
    family_color,
    binomial,
    n_contigs,
    total_mb,
    depths,
    lengths,
    mean_coverage,
    circ_indicator
  )

# Build gt table
gt_tbl <- tbl_data %>%
  gt() %>%
  # Column labels
  cols_label(
    rank = "",
    domain_abbr = "Domain",
    family = "Family",
    binomial = "Species",
    n_contigs = "n",
    total_mb = "Mb",
    depths = "Depth",
    lengths = "Mb",
    mean_coverage = "Cov%",
    circ_indicator = "Circ"
  ) %>%
  # Hide helper columns
  cols_hide(c(domain_color, family_color)) %>%
  # Format coverage as percentage
  fmt_percent(
    columns = mean_coverage,
    decimals = 0
  ) %>%
  # Style domain with background color
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
      paste0(
        "<span style='border-left: 4px solid ", colors,
        "; padding-left: 6px;'>", x, "</span>"
      )
    }
  ) %>%
  # Style binomial as italic
  text_transform(
    locations = cells_body(columns = binomial),
    fn = function(x) {
      paste0("<em>", x, "</em>")
    }
  ) %>%
  # Add sparklines for depth distribution
  gt_plt_dist(
    column = depths,
    type = "boxplot",
    line_color = "darkgray",
    fill_color = "lightgray",
    same_limit = TRUE
  ) %>%
  # Add sparklines for length distribution (in Mb)
  gt_plt_dist(
    column = lengths,
    type = "boxplot",
    line_color = "steelblue",
    fill_color = "lightblue",
    same_limit = TRUE
  ) %>%
  # Center columns
  cols_align(
    align = "center",
    columns = c(rank, n_contigs, circ_indicator, mean_coverage, total_mb)
  ) %>%
  # Right-align numeric columns
  cols_align(
    align = "right",
    columns = c(total_mb)
  ) %>%
  # Column widths
  cols_width(
    rank ~ px(30),
    domain_abbr ~ px(50),
    family ~ px(120),
    binomial ~ px(180),
    n_contigs ~ px(30),
    total_mb ~ px(50),
    depths ~ px(90),
    lengths ~ px(90),
    mean_coverage ~ px(45),
    circ_indicator ~ px(35)
  ) %>%
  # Table styling
  tab_style(
    style = cell_text(size = px(11)),
    locations = cells_body()
  ) %>%
  tab_style(
    style = cell_text(size = px(11), weight = "bold"),
    locations = cells_column_labels()
  ) %>%
  # Reduce row padding for compact display
  tab_options(
    data_row.padding = px(4),
    column_labels.padding = px(6),
    table.font.size = px(11)
  ) %>%
  # Title and subtitle
tab_header(
    title = md(paste0("**Top ", nrow(tbl_data), " Contaminants**")),
    subtitle = md(paste0("Ranked by abundance (depth \u00D7 length) \u2014 ", plot_suffix))
  ) %>%
  # Footnotes for indicators
  tab_footnote(
    footnote = "\u25CF = all circular, \u25D0 = mixed, \u25CB = all linear",
    locations = cells_column_labels(columns = circ_indicator)
  ) %>%
  # Spanner for distribution columns
  tab_spanner(
    label = "Distribution",
    columns = c(depths, lengths)
  ) %>%
  # Source note with family legend
  tab_source_note(
    source_note = md(paste0(
      "**Families:** ",
      paste(
        sapply(names(family_colors), function(f) {
          paste0("<span style='border-left: 4px solid ", family_colors[f],
                 "; padding-left: 4px; margin-right: 12px;'>", f, "</span>")
        }),
        collapse = " "
      )
    ))
  )

# Save as HTML (PDF not supported due to SVG sparklines requiring browser engine)
gtsave(gt_tbl, out_html)
message("Contaminant table saved to: ", out_html)
