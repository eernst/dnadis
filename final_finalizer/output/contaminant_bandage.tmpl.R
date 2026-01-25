#!/usr/bin/env Rscript
# Bandage-like contamination visualization
# Shows individual contigs as circles (circular) or pills (linear)

if (!requireNamespace("pacman", quietly = TRUE)) {
 install.packages("pacman", repos = "https://cloud.r-project.org")
}
library(pacman)
pacman::p_load(
 readr, dplyr, stringr, ggplot2, ggforce,
 scales, patchwork
)

# Placeholders - replaced by Python
contam_file <- "__CONTAMINANTS_TSV__"
out_pdf     <- "__OUTPDF__"
plot_suffix <- "__SUFFIX__"

base_family <- "Helvetica"
base_font_pt <- 8

# Custom byte formatting function (avoids locale issues with scales::label_bytes)
format_bytes <- function(x) {
 sapply(x, function(bytes) {
   if (is.na(bytes) || bytes < 0) return(NA_character_)
   if (bytes < 1024) return(paste0(bytes, " B"))
   if (bytes < 1024^2) return(sprintf("%.1f KB", bytes / 1024))
   if (bytes < 1024^3) return(sprintf("%.1f MB", bytes / 1024^2))
   return(sprintf("%.1f GB", bytes / 1024^3))
 })
}

# Read and filter data
df <- read_tsv(contam_file, show_col_types = FALSE)

# Check if we have any data
if (nrow(df) == 0) {
 message("No contaminants detected. Skipping bandage plot generation.")
 quit(status = 0)
}

# Filter to contigs with valid depth data (required for depth-ordered layout)
# Note: Coverage filtering already applied in Python before writing TSV
df <- df %>%
 filter(!is.na(depth_mean) & depth_mean > 0)

if (nrow(df) == 0) {
 message("No contaminants with valid depth data. Skipping bandage plot.")
 quit(status = 0)
}

# Determine contig topology from name suffix
# Use power scaling (length^0.4) for better size differentiation than log
df <- df %>%
 mutate(
   topology = if_else(str_detect(contig, "c$"), "circular", "linear"),
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
   # Power scaling for size differentiation (length^0.65 for strong contrast)
   # Normalize to reasonable visual range
   scaled_size = (length / max(length))^0.65,
   radius = scaled_size * 1.5 + 0.15,  # Scale to visual range [0.15, 1.65]
   # Format label text
   size_label = format_bytes(length)
 ) %>%
 arrange(desc(depth_mean))  # Order by depth (highest first)

# Fixed thickness for rings and pills (not proportional to size)
fixed_thickness <- 0.12

# Family color palette (same as treemap)
families <- unique(df$family)
n_families <- length(families)
if (n_families <= 8) {
 family_colors <- scales::brewer_pal(palette = "Set2")(max(n_families, 3))[1:n_families]
} else if (n_families <= 12) {
 family_colors <- scales::brewer_pal(palette = "Set3")(max(n_families, 3))[1:n_families]
} else {
 family_colors <- scales::hue_pal()(n_families)
}
names(family_colors) <- families

# Function to create grid layout (left-to-right, top-to-bottom by input order)
# Bottom-aligns contigs within each row, returns row baseline for label alignment
grid_layout <- function(radii, ncol = 4, padding = 0.15) {
 n <- length(radii)
 if (n == 0) return(data.frame(x = numeric(0), y = numeric(0), row_baseline = numeric(0)))

 # Calculate cell dimensions
 nrow <- ceiling(n / ncol)

 # Assign each item to a cell
 col_idx <- ((seq_len(n) - 1) %% ncol) + 1
 row_idx <- ((seq_len(n) - 1) %/% ncol) + 1

 # Calculate column widths (max radius * 2 + padding for each column)
 col_widths <- sapply(1:ncol, function(c) {
   items_in_col <- which(col_idx == c)
   if (length(items_in_col) == 0) return(0)
   max(radii[items_in_col]) * 2 + padding
 })

 # Calculate row heights based on max radius in row
 # Label space is fixed (text height), shape height varies
 label_space <- 0.55  # Space for labels below shapes
 row_max_radii <- sapply(1:nrow, function(r) {
   items_in_row <- which(row_idx == r)
   if (length(items_in_row) == 0) return(0)
   max(radii[items_in_row])
 })
 row_heights <- row_max_radii * 2 + padding + label_space

 # Calculate cumulative positions
 col_starts <- c(0, cumsum(col_widths[-ncol]))
 row_starts <- c(0, cumsum(row_heights[-nrow]))

 # Calculate center positions for each item
 x <- sapply(seq_len(n), function(i) {
   col_starts[col_idx[i]] + col_widths[col_idx[i]] / 2
 })

 # Y positions: bottom-aligned within each row
 # All contigs in a row have their bottom (y - radius) at the same level
 # The row bottom is at -(row_start + row_height - label_space)
 y <- sapply(seq_len(n), function(i) {
   row <- row_idx[i]
   row_bottom <- row_starts[row] + row_max_radii[row] * 2 + padding / 2
   # Position center so bottom of contig aligns with row bottom
   -(row_bottom - radii[i])
 })

 # Row baseline: the y-position where all labels in a row should start
 # This is the bottom of the tallest circle in the row (for consistent label alignment)
 row_baseline <- sapply(seq_len(n), function(i) {
   row <- row_idx[i]
   row_bottom <- row_starts[row] + row_max_radii[row] * 2 + padding / 2
   # Baseline is below the tallest shape (using circle bottom, not pill)
   -(row_bottom)
 })

 data.frame(x = x, y = y, row_baseline = row_baseline)
}

# Function to create plot for one domain
plot_domain <- function(df_domain, domain_name) {
 if (nrow(df_domain) == 0) return(NULL)

 # Grid layout: left-to-right, top-to-bottom by depth order
 ncol <- min(4, nrow(df_domain))  # Max 4 columns, fewer if less data
 layout <- grid_layout(df_domain$radius, ncol = ncol, padding = 0.2)
 df_domain <- df_domain %>%
   mutate(
     x = layout$x,
     y = layout$y,
     row_baseline = layout$row_baseline,  # Consistent label y-position per row
     pill_thickness = fixed_thickness  # Fixed thickness for all pills
   )

 # Domain-specific stats
 n_contigs <- nrow(df_domain)
 total_mb <- sum(df_domain$length) / 1e6

 # Separate circular and linear contigs
 df_circular <- df_domain %>% filter(topology == "circular")
 df_linear <- df_domain %>% filter(topology == "linear")

 # Calculate pill dimensions for drawing
 df_linear <- df_linear %>%
   mutate(
     pill_length = radius * 2  # Horizontal extent matches circle diameter
   )

 # Build pill polygons for linear contigs
 if (nrow(df_linear) > 0) {
   pill_polys <- do.call(rbind, lapply(1:nrow(df_linear), function(i) {
     row <- df_linear[i, ]
     half_len <- row$pill_length / 2
     half_thick <- row$pill_thickness / 2
     r <- half_thick

     # Number of points for semicircles
     n_arc <- 20

     # Right semicircle
     theta_right <- seq(-pi/2, pi/2, length.out = n_arc)
     right_x <- row$x + half_len - r + r * cos(theta_right)
     right_y <- row$y + r * sin(theta_right)

     # Left semicircle
     theta_left <- seq(pi/2, 3*pi/2, length.out = n_arc)
     left_x <- row$x - half_len + r + r * cos(theta_left)
     left_y <- row$y + r * sin(theta_left)

     data.frame(
       x = c(right_x, left_x),
       y = c(right_y, left_y),
       id = i,
       family = row$family,
       contig = row$contig
     )
   }))
   # Join back to get all columns for the legend
   pill_polys <- pill_polys %>%
     left_join(df_linear %>% select(contig, family), by = c("contig", "family"))
 } else {
   pill_polys <- NULL
 }

 # Label offset below shapes
 label_offset <- 0.15

 # Create plot
 p <- ggplot()

 # Draw circular contigs as thin rings (fixed thickness)
 if (nrow(df_circular) > 0) {
   p <- p +
     geom_circle(
       data = df_circular,
       aes(x0 = x, y0 = y, r = radius, fill = family),
       color = "white", linewidth = 0.3
     ) +
     geom_circle(
       data = df_circular,
       aes(x0 = x, y0 = y, r = pmax(radius - fixed_thickness, 0.02)),
       fill = "white", color = NA
     )
 }

 # Draw linear contigs as pills (stadium shapes)
 if (!is.null(pill_polys) && nrow(pill_polys) > 0) {
   p <- p +
     geom_polygon(
       data = pill_polys,
       aes(x = x, y = y, group = id, fill = family),
       color = "white", linewidth = 0.3
     )
 }

 # Labels beneath contigs (use row_baseline for consistent alignment across row)
 p <- p +
   geom_text(
     data = df_domain,
     aes(x = x, y = row_baseline - label_offset,
         label = paste0(binomial, "\n",
                       sprintf("%.0fx", depth_mean), ", ", size_label)),
     size = (base_font_pt - 2) / .pt,
     fontface = "italic",
     family = base_family,
     lineheight = 0.9,
     vjust = 1
   ) +
   scale_fill_manual(values = family_colors, name = "Family") +
   # Set explicit coordinate limits based on data extent to prevent excessive scaling
   coord_fixed(
     xlim = c(min(df_domain$x - df_domain$radius) - 0.3,
              max(df_domain$x + df_domain$radius) + 0.3),
     ylim = c(min(df_domain$row_baseline) - 0.8,  # Space for labels
              max(df_domain$y + df_domain$radius) + 0.2),
     expand = FALSE
   ) +
   theme_void(base_family = base_family) +
   labs(
     title = domain_name,
     subtitle = paste0(n_contigs, " contigs, ", sprintf("%.2f", total_mb), " Mb")
   ) +
   theme(
     plot.title = element_text(hjust = 0.5, size = base_font_pt + 1, face = "bold"),
     plot.subtitle = element_text(hjust = 0.5, size = base_font_pt - 1, color = "grey40")
   )

 return(p)
}

# Domain summary for main title
domain_summary <- df %>%
 group_by(kingdom) %>%
 summarise(bp = sum(length), .groups = "drop") %>%
 arrange(desc(bp)) %>%
 mutate(pct = sprintf("%.0f%%", bp / sum(bp) * 100)) %>%
 mutate(label = paste0(kingdom, " (", pct, ")")) %>%
 pull(label) %>%
 paste(collapse = ", ")

# Gating information (coverage filtering applied upstream in Python)
gating_info <- "depth > 0"

# Create faceted plot by domain
domains <- unique(df$kingdom)
plots <- lapply(domains, function(d) {
 plot_domain(df %>% filter(kingdom == d), d)
})
plots <- plots[!sapply(plots, is.null)]

if (length(plots) == 0) {
 message("No valid domains to plot. Skipping bandage plot.")
 quit(status = 0)
}

# Calculate appropriate height based on data
# Estimate rows needed per domain
rows_per_domain <- sapply(domains, function(d) {
 n <- sum(df$kingdom == d)
 max(ceiling(n / 4), 1)  # 4 columns, minimum 1 row
})
total_rows <- sum(rows_per_domain)
plot_height <- 1.5 + total_rows * 2.0

# Combine plots with heights proportional to rows needed
combined <- wrap_plots(plots, ncol = 1, heights = rows_per_domain) +
 plot_annotation(
   title = paste0("Contaminants (", plot_suffix, ")"),
   subtitle = paste0(domain_summary, "\n", gating_info),
   theme = theme(
     plot.title = element_text(hjust = 0.5, size = base_font_pt + 2, face = "bold", family = base_family,
                               margin = margin(t = 0, b = 2)),
     plot.subtitle = element_text(hjust = 0.5, size = base_font_pt, color = "grey40", family = base_family,
                                  margin = margin(t = 0, b = 5)),
     plot.margin = margin(t = 15, r = 2, b = 15, l = 2)
   )
 ) &
 theme(
   legend.position = "bottom",
   legend.margin = margin(t = 15, b = 0),
   legend.title = element_text(size = base_font_pt - 1, margin = margin(r = 8)),
   plot.margin = margin(t = 0, r = 2, b = 0, l = 2)
 ) &
 guides(fill = guide_legend(nrow = 2, keywidth = unit(0.4, "cm"), keyheight = unit(0.4, "cm")))

# Save (max 7" wide for uniformity with other plots)
ggsave(out_pdf, combined, width = 7, height = min(plot_height, 14), units = "in", dpi = 300)
message("Bandage-like contamination plot saved to: ", out_pdf)
