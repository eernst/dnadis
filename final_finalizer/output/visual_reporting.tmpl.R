#!/usr/bin/env Rscript

if (!requireNamespace("pacman", quietly = TRUE)) {
  install.packages("pacman", repos = "https://cloud.r-project.org")
}
library(pacman)
pacman::p_load(
  readr, dplyr, stringr, ggplot2, ggnewscale, tibble, tidyr, patchwork,
  grid, ggiraph, htmlwidgets
)

summary_file <- "__SUMMARY__"
ref_file     <- "__REF__"
seg_file     <- "__SEGMENTS__"
ev_file      <- "__EVIDENCE__"
macro_file   <- "__MACRO__"
out_pdf      <- "__OUTPDF__"
out_html     <- "__OUTHTML__"
plot_html    <- as.logical("__PLOTHTML__")
chr_like_minlen <- as.numeric("__CHRLIKE__")
plot_title_suffix <- "__SUFFIX__"
base_family <- "Helvetica"
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

df <- read_tsv(summary_file, show_col_types = FALSE)

ref_lengths <- read_tsv(ref_file, show_col_types = FALSE)

seg <- readr::read_tsv(
  seg_file,
  show_col_types = FALSE,
  col_types = cols(
    contig = col_character(),
    contig_len = col_double(),
    target_chrom_id = col_character(),
    target_subgenome = col_character(),
    strand = col_character(),
    chain_id = col_double(),
    qstart = col_double(),
    qend = col_double()
  )
)

macro <- readr::read_tsv(
  macro_file,
  show_col_types = FALSE,
  col_types = cols(
    contig = col_character(),
    contig_len = col_double(),
    ref_id = col_character(),
    chrom_id = col_character(),
    subgenome = col_character(),
    strand = col_character(),
    chain_id = col_double(),
    qstart = col_double(),
    qend = col_double(),
    qspan_bp = col_double(),
    union_bp = col_double(),
    matches = col_double(),
    aln_len = col_double(),
    identity = col_double(),
    score = col_double(),
    n_segments = col_double(),
    gene_count_chain = col_double()
  )
) %>%
  mutate(
    target_chrom_id = chrom_id,
    target_subgenome = subgenome
  )

ev <- readr::read_tsv(
  ev_file,
  show_col_types = FALSE,
  col_types = cols(
    contig = col_character(),
    ref_id = col_character(),
    chrom_id = col_character(),
    subgenome = col_character(),
    union_bp = col_double(),
    identity = col_double(),
    score_all = col_double()
  )
)

# Subgenome suffix labels present in the ref lengths file (single-letter A-Z)
sg_levels <- ref_lengths %>%
  distinct(subgenome) %>%
  pull(subgenome) %>%
  as.character()

sg_levels <- sg_levels[grepl("^[A-Z]$", sg_levels)]
sg_levels <- sort(sg_levels)

n_sets <- length(sg_levels)

if (n_sets > 4) {
  warning(paste0("Reference has ", n_sets, " chromosome sets; plotting supports up to 4. ",
                 "Will use first 4: ", paste(sg_levels[1:4], collapse=",")))
  sg_levels <- sg_levels[1:4]
  n_sets <- 4
}

has_subgenomes <- (n_sets >= 2)

# Custom palettes for up to 4 subgenomes
# Blue, Orange, Teal (harmonious triad), Green
pal_dark  <- c("#1F77B4", "#FF7F0E", "#17A589", "#2E7D32")
pal_light <- c("#9ECAE1", "#FDD0A2", "#A3E4D7", "#A5D6A7")

plot_levels <- if (has_subgenomes) sg_levels else c("G")

# Named color vectors
sg_cols_dark  <- setNames(pal_dark [seq_along(plot_levels)], plot_levels)
sg_cols_light <- setNames(pal_light[seq_along(plot_levels)], plot_levels)

col_dark  <- c(sg_cols_dark,  "NA"="grey45")
col_light <- c(sg_cols_light, "NA"="grey82")

if (has_subgenomes) {
  ref_sg <- ref_lengths %>%
    filter(subgenome %in% sg_levels) %>%
    mutate(
      chrom_label = chrom_id %>%
        str_replace("^chr", "") %>%
        str_replace("^Chr", ""),
      chrom_index = suppressWarnings(as.integer(chrom_label)),
      chrom_index = if_else(is.na(chrom_index), 9999L, chrom_index)
    )

  chrom_levels <- ref_sg %>%
    distinct(chrom_id, chrom_index) %>%
    arrange(chrom_index, chrom_id) %>%
    pull(chrom_id)

} else {
  # Non-hybrid reference: use chrom_id directly (chr1, chr2, ...)
  ref_sg <- ref_lengths %>%
    filter(!(ref_id %in% c("chrC","chrM"))) %>%
    mutate(
      chrom_label = chrom_id %>%
        str_replace("^chr", "") %>%
        str_replace("^Chr", ""),
      chrom_index = suppressWarnings(as.integer(chrom_label)),
      chrom_index = if_else(is.na(chrom_index), 9999L, chrom_index),
      subgenome = "G"
    )

  chrom_levels <- ref_sg %>%
    distinct(chrom_id, chrom_index) %>%
    arrange(chrom_index, chrom_id) %>%
    pull(chrom_id)

  seg <- seg %>%
    mutate(target_subgenome = "G")

  macro <- macro %>%
    mutate(target_subgenome = "G")
}

# Spacing factor between chromosome bins (1.0 = no extra space, 1.2 = 20% more)
chrom_spacing <- 1.15

x_levels <- c(chrom_levels, "Un")
x_tbl <- tibble(chrom_id = factor(x_levels, levels = x_levels),
                x_index = seq_along(x_levels) * chrom_spacing)

ref_lines_base <- ref_sg %>%
  mutate(chrom_id = factor(chrom_id, levels = x_levels)) %>%
  left_join(x_tbl, by = "chrom_id") %>%
  filter(!is.na(x_index)) %>%
  mutate(ref_len_mb = ref_len / 1e6)

df_plot <- df %>%
  mutate(
    contig_len = if_else(is.na(length), 0, as.numeric(length)),
    assigned_subgenome = case_when(
      !is.na(assigned_subgenome) & assigned_subgenome %in% plot_levels ~ as.character(assigned_subgenome),
      TRUE ~ "NA"
    ),
    assigned_chrom_id  = if_else(assigned_chrom_id %in% chrom_levels, assigned_chrom_id, "Un"),
    chrom_id = factor(assigned_chrom_id, levels = x_levels),
    best_identity = as.numeric(best_identity),
    # Orientation: TRUE if contig should be shown reversed to match reference
    is_reversed = if ("reversed" %in% names(.)) reversed == "yes" else FALSE,
    # Telomere columns (may not exist in older output files)
    # Raw values before orientation adjustment
    raw_5p_telo = if ("has_5p_telomere" %in% names(.)) has_5p_telomere == "yes" else FALSE,
    raw_3p_telo = if ("has_3p_telomere" %in% names(.)) has_3p_telomere == "yes" else FALSE,
    # Flip telomere positions for reversed contigs so they appear at correct visual ends
    has_5p_telo = if_else(is_reversed, raw_3p_telo, raw_5p_telo),
    has_3p_telo = if_else(is_reversed, raw_5p_telo, raw_3p_telo)
  ) %>%
  filter(contig_len >= chr_like_minlen)

# Evidence (for radar)
ev <- ev %>%
  mutate(
    subgenome = case_when(
      !is.na(subgenome) & subgenome %in% sg_levels ~ as.character(subgenome),
      TRUE ~ "NA"
    )
  ) %>%
  filter(subgenome %in% sg_levels)

# ----------------------------
# Top: Contig composition from segments
# ----------------------------
df_slots <- df_plot %>%
  mutate(contig_len_mb = contig_len / 1e6) %>%
  filter(assigned_chrom_id != "Un" | classification == "chrom_unassigned") %>%
  left_join(x_tbl, by = "chrom_id") %>%
  filter(!is.na(x_index)) %>%
  group_by(chrom_id) %>%
  arrange(assigned_subgenome, desc(contig_len), contig, .by_group = TRUE) %>%
  mutate(
    n_in_bin = n(),
    slot = row_number(),
    total_width = 0.80,
    dx = if_else(n_in_bin > 1, total_width / n_in_bin, 0),
    x_offset = if_else(n_in_bin > 1, (slot - (n_in_bin + 1) / 2) * dx, 0),
    x_plot = x_index + x_offset,
    # Use assigned_subgenome from df_plot (already extracted from summary file)
    # Must be a factor with levels matching col_dark for scale_color_manual to work
    ref_subgenome = factor(assigned_subgenome, levels = names(col_dark))
  ) %>%
  ungroup()

n_per_chrom <- df_slots %>% count(chrom_id, name = "n_contigs")
max_n <- max(n_per_chrom$n_contigs, na.rm = TRUE)
if (!is.finite(max_n) || max_n <= 0) max_n <- 1
band_xw <- max(0.02, min(0.08, 0.80 / max_n * 0.80))

# Reference lines: vertical lines to left of leftmost pill, one per subgenome with assigned contigs
# Basic data prep (x_ref calculation deferred until after pill_xw is defined)
ref_lines_data <- df_slots %>%
  group_by(chrom_id) %>%
  summarise(
    x_left = min(x_plot),  # leftmost pill position in this bin
    subgenomes_present = list(unique(assigned_subgenome)),
    .groups = "drop"
  ) %>%
  tidyr::unnest(subgenomes_present) %>%
  rename(assigned_subgenome = subgenomes_present) %>%
  left_join(
    ref_lines_base %>% select(chrom_id, subgenome, ref_len_mb),
    by = c("chrom_id" = "chrom_id", "assigned_subgenome" = "subgenome")
  ) %>%
  filter(!is.na(ref_len_mb)) %>%
  # Arrange subgenomes and count within each chrom_id
  group_by(chrom_id) %>%
  arrange(assigned_subgenome, .by_group = TRUE) %>%
  mutate(
    n_sg = n(),
    sg_idx = row_number()
  ) %>%
  ungroup()

pill_lwd <- dplyr::case_when(
  max_n <= 3  ~ 4,
  max_n <= 6  ~ 3,
  max_n <= 10 ~ 2.5,
  max_n <= 18 ~ 1.5,
  TRUE        ~ 1
) * 0.5

seg_plot <- seg %>%
  filter(contig %in% df_slots$original_name) %>%
  mutate(
    qstart_mb = qstart / 1e6,
    qend_mb   = qend   / 1e6,
    contig_len_mb = contig_len / 1e6
  ) %>%
  left_join(
    df_slots %>% select(original_name, chrom_id, x_plot, contig_len_mb_slot = contig_len_mb, is_reversed),
    by = c("contig" = "original_name")
  ) %>%
  filter(!is.na(x_plot)) %>%
  mutate(
    assigned_chrom_id = as.character(chrom_id),
    off_target = (target_chrom_id != assigned_chrom_id),

    max_len_mb = if_else(is.na(contig_len_mb_slot), contig_len_mb, contig_len_mb_slot),
    # Flip coordinates for reversed contigs to show in reference orientation
    qstart_flip = if_else(is_reversed, max_len_mb - qend_mb, qstart_mb),
    qend_flip   = if_else(is_reversed, max_len_mb - qstart_mb, qend_mb),
    qstart_mb = pmax(0, pmin(qstart_flip, max_len_mb)),
    qend_mb   = pmax(0, pmin(qend_flip,   max_len_mb))
  ) %>%
  filter(qend_mb > qstart_mb) %>%
  mutate(
    xmin = x_plot - band_xw/2,
    xmax = x_plot + band_xw/2,
    ymin = qstart_mb,
    ymax = qend_mb
  ) %>%
  mutate(
    tooltip = paste0(
      "contig: ", contig,
      "\nassigned: ", assigned_chrom_id,
      "\nsegment: ", target_chrom_id, " ", target_subgenome,
      "\nqstart: ", sprintf("%.3f", qstart_mb), " Mbp",
      "\nqend: ", sprintf("%.3f", qend_mb), " Mbp",
      "\nchain: ", chain_id
    )
  )

macro_plot <- macro %>%
  filter(contig %in% df_slots$original_name) %>%
  mutate(
    qstart_mb = qstart / 1e6,
    qend_mb   = qend   / 1e6,
    contig_len_mb = contig_len / 1e6
  ) %>%
  left_join(
    df_slots %>%
      transmute(
        original_name,
        assigned_chrom_id = as.character(chrom_id),
        x_plot,
        contig_len_mb_slot = contig_len_mb,
        is_reversed
      ),
    by = c("contig" = "original_name")
  ) %>%
  filter(!is.na(x_plot)) %>%
  mutate(
    off_target = (target_chrom_id != assigned_chrom_id),

    max_len_mb = if_else(is.na(contig_len_mb_slot), contig_len_mb, contig_len_mb_slot),
    # Flip coordinates for reversed contigs to show in reference orientation
    qstart_flip = if_else(is_reversed, max_len_mb - qend_mb, qstart_mb),
    qend_flip   = if_else(is_reversed, max_len_mb - qstart_mb, qend_mb),
    qstart_mb = pmax(0, pmin(qstart_flip, max_len_mb)),
    qend_mb   = pmax(0, pmin(qend_flip,   max_len_mb))
  ) %>%
  filter(qend_mb > qstart_mb) %>%
  mutate(
    xmin = x_plot - band_xw/2,
    xmax = x_plot + band_xw/2,
    ymin = qstart_mb,
    ymax = qend_mb
  ) %>%
  mutate(
    tooltip = paste0(
      "contig: ", contig,
      "\nassigned: ", assigned_chrom_id,
      "\nmacro: ", target_chrom_id, " ", target_subgenome,
      "\nqstart: ", sprintf("%.3f", qstart_mb), " Mbp",
      "\nqend: ", sprintf("%.3f", qend_mb), " Mbp",
      "\nidentity: ", sprintf("%.4f", identity),
      "\nscore: ", sprintf("%.3f", score),
      "\nchain: ", chain_id
    )
  )

# Width for pill backgrounds and alignment segments (same width for both)
pill_xw <- band_xw * 2.4

# Reference line parameters
ref_tick_offset <- pill_xw * 0.4  # gap between lines and pill
ref_line_width <- 0.4  # linewidth in mm
ref_line_spacing <- 0.04  # x-axis spacing between lines (data units)

# Calculate x positions for reference lines (now that pill_xw is defined)
ref_lines <- ref_lines_data %>%
  mutate(
    # Position from right to left (rightmost = closest to pill)
    x_ref = x_left - pill_xw/2 - ref_tick_offset - (n_sg - sg_idx) * ref_line_spacing * 2
  )

p_comp <- ggplot() +
  # Reference chromosome length lines (vertical, to left of leftmost pill per bin)
  geom_segment(
    data = ref_lines,
    aes(
      x = x_ref, xend = x_ref,
      y = 0, yend = ref_len_mb,
      color = assigned_subgenome
    ),
    linewidth = ref_line_width,
    alpha = 0.9,
    show.legend = FALSE
  ) +
  scale_color_manual(values = col_light, guide = "none", drop = FALSE) +
  ggnewscale::new_scale_color() +
  # Telomere indicators - slightly darker circles behind pill background
  # 5' telomere (bottom of contig)
  geom_point(
    data = df_slots %>% filter(has_5p_telo),
    aes(x = x_plot, y = 0),
    color = "grey75",
    size = pill_xw * 5,
    show.legend = FALSE
  ) +
  # 3' telomere (top of contig)
  geom_point(
    data = df_slots %>% filter(has_3p_telo),
    aes(x = x_plot, y = contig_len_mb),
    color = "grey75",
    size = pill_xw * 5,
    show.legend = FALSE
  ) +
  # Pill backgrounds - all light gray (same width as alignment blocks)
  geom_rect(
    data = df_slots,
    aes(xmin = x_plot - pill_xw/2, xmax = x_plot + pill_xw/2,
        ymin = 0, ymax = contig_len_mb),
    fill = "grey88",
    color = NA,
    show.legend = FALSE
  ) +
  ggnewscale::new_scale_fill() +
  # Alignment segments (full pill width)
  geom_rect(
    data = seg_plot %>% filter(!off_target),
    aes(xmin = x_plot - pill_xw/2, xmax = x_plot + pill_xw/2,
        ymin = ymin, ymax = ymax, fill = target_subgenome),
    color = NA,
    alpha = 1.0
  ) +
  geom_rect(
    data = seg_plot %>% filter(off_target),
    aes(xmin = x_plot - pill_xw/2, xmax = x_plot + pill_xw/2,
        ymin = ymin, ymax = ymax),
    fill = "red",
    color = NA,
    alpha = 1.0
  ) +
  scale_fill_manual(values = col_light, guide = "none", drop = FALSE) +
  ggnewscale::new_scale_fill() +
  # Macro blocks (full pill width)
  geom_rect(
    data = macro_plot %>% filter(!off_target),
    aes(xmin = x_plot - pill_xw/2, xmax = x_plot + pill_xw/2,
        ymin = ymin, ymax = ymax, fill = target_subgenome),
    color = NA,
    alpha = 1.0,
    show.legend = FALSE
  ) +
  geom_rect(
    data = macro_plot %>% filter(off_target),
    aes(xmin = x_plot - pill_xw/2, xmax = x_plot + pill_xw/2,
        ymin = ymin, ymax = ymax),
    fill = "red",
    color = NA,
    alpha = 1.0,
    show.legend = FALSE
  ) +
  scale_fill_manual(values = col_dark, guide = "none", drop = FALSE) +
  guides(
    fill = guide_legend(
      title.position = "top",
      title.hjust = 0.5,
      nrow = 1,
      byrow = TRUE,
      override.aes = list(alpha = 1)
    )
  ) +
  scale_x_continuous(
    breaks = x_tbl$x_index,
    labels = str_replace(as.character(x_tbl$chrom_id), "^chr", ""),
    expand = expansion(mult = c(0.01, 0.01))
  ) +
  theme_classic(base_family = base_family, base_size = base_font_pt) +
  axis_theme +
  theme(
    axis.ticks = element_blank(),
    plot.title  = element_text(hjust = 0.5, size = base_font_pt, family = base_family),
    panel.grid = element_blank(),
    legend.position = c(0.985, 0.985),
    legend.justification = c(1, 1)
  ) +
  labs(
    x = "assigned chromosome",
    y = "Contig position (Mbp)",
    title = paste0("Contig composition (", plot_title_suffix, ")")
  )

p_comp <- apply_legend_theme(p_comp, text_pt = 6, key_pt = 6, tight = TRUE) +
  theme(legend.margin = margin(1, 1, 1, 1))

# Interactive version (tooltips on segment/macro blocks)
if (plot_html) {
  p_comp_html <- ggplot() +
    # Reference chromosome length lines (vertical, to left of leftmost pill per bin)
    geom_segment(
      data = ref_lines,
      aes(
        x = x_ref, xend = x_ref,
        y = 0, yend = ref_len_mb,
        color = assigned_subgenome
      ),
      linewidth = ref_line_width,
      alpha = 0.9,
      show.legend = FALSE
    ) +
    scale_color_manual(values = col_light, guide = "none", drop = FALSE) +
    ggnewscale::new_scale_color() +
    # Telomere indicators - slightly darker circles behind pill background
    # 5' telomere (bottom of contig)
    ggiraph::geom_point_interactive(
      data = df_slots %>% filter(has_5p_telo),
      aes(x = x_plot, y = 0,
          tooltip = paste0(original_name, "\n5' telomere detected")),
      color = "grey75",
      size = pill_xw * 5,
      show.legend = FALSE
    ) +
    # 3' telomere (top of contig)
    ggiraph::geom_point_interactive(
      data = df_slots %>% filter(has_3p_telo),
      aes(x = x_plot, y = contig_len_mb,
          tooltip = paste0(original_name, "\n3' telomere detected")),
      color = "grey75",
      size = pill_xw * 5,
      show.legend = FALSE
    ) +
    # Pill backgrounds - all light gray (same width as alignment blocks)
    geom_rect(
      data = df_slots,
      aes(xmin = x_plot - pill_xw/2, xmax = x_plot + pill_xw/2,
          ymin = 0, ymax = contig_len_mb),
      fill = "grey88",
      color = NA,
      show.legend = FALSE
    ) +
    ggnewscale::new_scale_fill() +
    # Alignment segments (full pill width)
    ggiraph::geom_rect_interactive(
      data = seg_plot %>% filter(!off_target),
      aes(xmin = x_plot - pill_xw/2, xmax = x_plot + pill_xw/2,
          ymin = ymin, ymax = ymax, fill = target_subgenome, tooltip = tooltip),
      color = NA,
      alpha = 1.0
    ) +
    ggiraph::geom_rect_interactive(
      data = seg_plot %>% filter(off_target),
      aes(xmin = x_plot - pill_xw/2, xmax = x_plot + pill_xw/2,
          ymin = ymin, ymax = ymax, tooltip = tooltip),
      fill = "red",
      color = NA,
      alpha = 1.0
    ) +
    scale_fill_manual(values = col_light, guide = "none", drop = FALSE) +
    ggnewscale::new_scale_fill() +
    # Macro blocks (full pill width)
    ggiraph::geom_rect_interactive(
      data = macro_plot %>% filter(!off_target),
      aes(xmin = x_plot - pill_xw/2, xmax = x_plot + pill_xw/2,
          ymin = ymin, ymax = ymax, fill = target_subgenome, tooltip = tooltip),
      color = NA,
      alpha = 1.0,
      show.legend = FALSE
    ) +
    ggiraph::geom_rect_interactive(
      data = macro_plot %>% filter(off_target),
      aes(xmin = x_plot - pill_xw/2, xmax = x_plot + pill_xw/2,
          ymin = ymin, ymax = ymax, tooltip = tooltip),
      fill = "red",
      color = NA,
      alpha = 1.0,
      show.legend = FALSE
    ) +
    scale_fill_manual(values = col_dark, guide = "none", drop = FALSE) +
    guides(
      fill = guide_legend(
        title.position = "top",
        title.hjust = 0.5,
        nrow = 1,
        byrow = TRUE,
        override.aes = list(alpha = 1)
      )
    ) +
    scale_x_continuous(
      breaks = x_tbl$x_index,
      labels = str_replace(as.character(x_tbl$chrom_id), "^chr", ""),
      expand = expansion(mult = c(0.01, 0.01))
    ) +
    theme_classic(base_family = base_family, base_size = base_font_pt) +
    axis_theme +
    theme(
      axis.ticks = element_blank(),
      plot.title  = element_text(hjust = 0.5, size = base_font_pt, family = base_family),
      panel.grid = element_blank(),
      legend.position = c(0.985, 0.985),
      legend.justification = c(1, 1)
    ) +
    labs(
      x = "assigned chromosome",
      y = "Contig position (Mbp)",
      title = paste0("Contig composition (", plot_title_suffix, ")")
    )

  p_comp_html <- apply_legend_theme(p_comp_html, text_pt = 6, key_pt = 6, tight = TRUE) +
    theme(legend.margin = margin(1, 1, 1, 1))
}

# ----------------------------
# Bottom-left: radar scatter (hybrid only)
# ----------------------------
if (has_subgenomes) {

  # Support metric: union_bp * identity
  radar_support <- ev %>%
    mutate(support = pmax(0, union_bp) * pmax(0, pmin(1, identity))) %>%
    group_by(contig, subgenome) %>%
    summarise(support = sum(support, na.rm = TRUE), .groups = "drop") %>%
    tidyr::complete(contig, subgenome = plot_levels, fill = list(support = 0)) %>%
    group_by(contig) %>%
    mutate(
      wsum = sum(support),
      w = if_else(wsum > 0, support / wsum, NA_real_)
    ) %>%
    ungroup()

  get_vertices <- function(plot_levels) {
    n <- length(plot_levels)
    if (n == 2) {
      tibble(subgenome = plot_levels, vx = c(-0.5, 0.5), vy = c(0, 0))
    } else if (n == 3) {
      tibble(subgenome = plot_levels,
             vx = c(-0.5, 0.5, 0.0),
             vy = c( 0.0, 0.0, sqrt(3)/2))
    } else if (n == 4) {
      tibble(subgenome = plot_levels,
             vx = c(-0.5,  0.5,  0.5, -0.5),
             vy = c(-0.5, -0.5,  0.5,  0.5))
    } else {
      stop("get_vertices expects 2..4 subgenomes")
    }
  }
  verts <- get_vertices(plot_levels)

  label_verts <- if (n_sets == 2) {
    verts %>%
      mutate(
        lx = if_else(vx < 0, -0.58, 0.58),
        ly = 0,
        hjust = if_else(vx < 0, 1, 0),
        vjust = 0.5
      )
  } else if (n_sets == 3) {
    # Triangle: offset labels outward from center
    verts %>%
      mutate(
        # Push labels away from centroid (0, ~0.29)
        centroid_y = sqrt(3)/6,
        dir_x = vx - 0,
        dir_y = vy - centroid_y,
        dir_len = sqrt(dir_x^2 + dir_y^2),
        lx = vx + 0.12 * dir_x / dir_len,
        ly = vy + 0.12 * dir_y / dir_len,
        hjust = 0.5,
        vjust = 0.5
      ) %>%
      select(-centroid_y, -dir_x, -dir_y, -dir_len)
  } else {
    # Square (4 subgenomes): offset labels to corners
    verts %>%
      mutate(
        lx = vx * 1.2,
        ly = vy * 1.2,
        hjust = 0.5,
        vjust = 0.5
      )
  }

  radar_xy <- radar_support %>%
    left_join(verts, by = "subgenome") %>%
    group_by(contig) %>%
    summarise(
      rx = sum(w * vx, na.rm = TRUE),
      ry = sum(w * vy, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    left_join(df_plot %>% select(original_name, assigned_subgenome, classification),
              by = c("contig" = "original_name")) %>%
    mutate(
      assigned_subgenome = if_else(is.na(assigned_subgenome), "NA", as.character(assigned_subgenome)),
      classification = if_else(is.na(classification), "NA", as.character(classification))
    ) %>%
    filter(assigned_subgenome %in% plot_levels | classification == "chrom_unassigned") %>%
    mutate(assigned_subgenome = factor(assigned_subgenome, levels = c(plot_levels, "NA")))

  frame <- if (n_sets == 2) {
    tibble(x = c(-0.5, 0.5), y = c(0, 0))
  } else {
    verts %>%
      transmute(x = vx, y = vy) %>%
      bind_rows(.[1,])
  }

  radar_xlim <- c(-0.6, 0.6)
  radar_ylim <- if (n_sets == 2) {
    c(-0.6, 0.6)
  } else {
    c(-0.15, max(verts$vy) + 0.06)
  }
  axis_pos <- (0 - radar_ylim[1]) / (radar_ylim[2] - radar_ylim[1])
  legend_pos <- c(0.5, pmax(0.02, axis_pos - 0.08))

  p_radar <- ggplot() +
    geom_path(data = frame, aes(x = x, y = y), linewidth = 0.4, color = "grey40") +
    geom_text(
      data = label_verts,
      aes(x = lx, y = ly, label = subgenome, hjust = hjust, vjust = vjust),
      color = "grey20",
      size = base_font_pt / .pt,
      family = base_family
    ) +
    geom_point(
      data = radar_xy %>% filter(!is.na(rx), !is.na(ry)),
      aes(x = rx, y = ry, color = assigned_subgenome),
      shape = 16, alpha = 0.80, size = 1.6, show.legend = TRUE
    ) +
    scale_color_manual(
      values = col_dark,
      breaks = "NA",
      labels = "Unassigned contigs",
      drop = FALSE
    ) +
    guides(color = guide_legend(title = NULL, override.aes = list(size = 2))) +
    coord_equal(clip = "off", xlim = radar_xlim, ylim = radar_ylim) +
    theme_classic(base_family = base_family, base_size = base_font_pt) +
    theme(
      plot.title = element_text(hjust = 0.5, size = base_font_pt, family = base_family),
      axis.title = element_blank(),
      axis.text  = element_blank(),
      axis.ticks = element_blank(),
      axis.line = element_blank(),
      panel.grid = element_blank(),
      plot.margin = margin(1.5, 1.5, 1.5, 1.5)
    ) +
    labs(title = "Chromosome-set support", x = NULL, y = NULL)

  p_radar <- apply_legend_theme(p_radar, text_pt = 6, key_pt = 6, tight = TRUE) +
    theme(legend.position = legend_pos, legend.justification = c(0.5, 1))

} else {
  p_radar <- ggplot() + theme_void(base_family = base_family, base_size = base_font_pt) + labs(title = "Subgenome support (n/a)")
}

# ----------------------------
# Bottom-right: identity vs assigned subgenome (or single-genome)
# ----------------------------
if (has_subgenomes) {

  df_id_assigned <- df_plot %>%
    filter(!is.na(best_identity), assigned_subgenome %in% plot_levels) %>%
    mutate(assigned_sg = factor(assigned_subgenome, levels = plot_levels))

  p_id <- ggplot(df_id_assigned, aes(x = assigned_sg, y = pmax(0, pmin(1, best_identity)), color = assigned_sg)) +
    geom_jitter(width = 0.18, height = 0, alpha = 0.55, size = 1.2, show.legend = FALSE) +
    stat_summary(fun = median, geom = "point", shape = 95, size = 7, show.legend = FALSE) +
    scale_color_manual(values = col_dark, drop = FALSE) +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
    theme_classic(base_family = base_family, base_size = base_font_pt) +
    axis_theme +
    theme(
      plot.title = element_text(hjust = 0.5, size = base_font_pt, family = base_family),
      axis.ticks = element_blank(),
      panel.grid = element_blank(),
      plot.margin = margin(5.5, 5.5, 5.5, 5.5)
    ) +
    labs(
      title = "Identity to assigned reference chromosome",
      x = "Assigned subgenome",
      y = "Alignment identity proxy"
    )

} else {

  df_id_assigned <- df_plot %>%
    filter(!is.na(best_identity), assigned_ref_id != "NA")

  df_id_assigned <- df_id_assigned %>%
    mutate(assigned_bin = factor("Genome", levels = "Genome"))

  p_id <- ggplot(df_id_assigned, aes(x = assigned_bin, y = pmax(0, pmin(1, best_identity)))) +
    geom_jitter(width = 0.18, height = 0, alpha = 0.55, size = 1.2, color = "grey35", show.legend = FALSE) +
    stat_summary(fun = median, geom = "point", shape = 95, size = 7, color = "grey15", show.legend = FALSE) +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
    theme_classic(base_family = base_family, base_size = base_font_pt) +
    axis_theme +
    theme(
      plot.title = element_text(hjust = 0.5, size = base_font_pt, family = base_family),
      axis.title.x = element_blank(),
      axis.ticks = element_blank(),
      panel.grid = element_blank(),
      plot.margin = margin(5.5, 5.5, 5.5, 5.5)
    ) +
    labs(
      title = "Identity to assigned reference chromosome",
      x = NULL,
      y = "Alignment identity proxy"
    )
}

if (has_subgenomes) {
  bottom_row <- (p_radar + p_id + plot_layout(widths = c(1, 1.35), guides = "keep")) &
    theme(plot.margin = margin(8.5, 8.5, 5.5, 8.5))
  full_plot <- p_comp / bottom_row + plot_layout(heights = c(1, 1))
} else {
  full_plot <- p_comp / p_id + plot_layout(heights = c(1, 1))
}

#message("macro rows: ", nrow(macro_plot), "  seg rows: ", nrow(seg_plot))
ggsave(out_pdf, plot = full_plot, width = 6, height = 6, units = "in", dpi = 300)

if (plot_html) {
  if (has_subgenomes) {
    bottom_row_html <- (p_radar + p_id + plot_layout(widths = c(1, 1.35), guides = "keep")) &
      theme(plot.margin = margin(8.5, 8.5, 5.5, 8.5))
    full_plot_html <- p_comp_html / bottom_row_html + plot_layout(heights = c(1, 1))
  } else {
    full_plot_html <- p_comp_html / p_id + plot_layout(heights = c(1, 1))
  }

  girafe_obj <- ggiraph::girafe(ggobj = full_plot_html)
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
}
