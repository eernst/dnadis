#!/usr/bin/env Rscript

if (!requireNamespace("pacman", quietly = TRUE)) {
  install.packages("pacman", repos = "https://cloud.r-project.org")
}
library(pacman)
pacman::p_load(
  readr, dplyr, stringr, ggplot2, ggnewscale, tibble, tidyr, patchwork,
  grid, ggiraph, htmlwidgets, ggrepel, ggokabeito, colorspace, showtext, sysfonts
)

# Font setup: try Helvetica Neue / Helvetica (macOS) then Liberation Sans (Linux)
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

summary_file <- "__SUMMARY__"
ref_file     <- "__REF__"
seg_file     <- "__SEGMENTS__"
ev_file      <- "__EVIDENCE__"
macro_file   <- "__MACRO__"
rdna_file    <- "__RDNA_ANNOTATIONS__"
out_pdf      <- "__OUTPDF__"
out_html     <- "__OUTHTML__"
plot_html    <- as.logical("__PLOTHTML__")
chr_like_minlen <- as.numeric("__CHRLIKE__")
plot_title_suffix <- "__SUFFIX__"
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

df <- read_tsv(summary_file, show_col_types = FALSE,
               col_types = cols(assigned_subgenome = col_character()))

ref_lengths <- read_tsv(ref_file, show_col_types = FALSE,
                        col_types = cols(subgenome = col_character()))

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

# Load rDNA annotations (optional)
rdna_annot <- if (nzchar(rdna_file) && file.exists(rdna_file)) {
  readr::read_tsv(
    rdna_file,
    show_col_types = FALSE,
    col_types = cols(
      contig = col_character(),
      start = col_double(),
      end = col_double(),
      strand = col_character(),
      identity = col_double(),
      consensus_coverage = col_double(),
      copy_type = col_character(),
      sub_features = col_character(),
      is_nor_candidate = col_character(),
      contig_classification = col_character()
    )
  ) %>%
    filter(copy_type == "full")  # Only full-copy rDNA
} else {
  tibble()  # Empty tibble if no rDNA annotations
}

# Subgenome suffix labels present in the ref lengths file (single-letter A-Z)
sg_levels <- ref_lengths %>%
  distinct(subgenome) %>%
  pull(subgenome) %>%
  as.character()

sg_levels <- sg_levels[grepl("^[A-Z]$", sg_levels)]
sg_levels <- sort(sg_levels)

n_sets <- length(sg_levels)

if (n_sets > 8) {
  warning(paste0("Reference has ", n_sets, " chromosome sets; plotting supports up to 8. ",
                 "Will use first 8: ", paste(sg_levels[1:8], collapse=",")))
  sg_levels <- sg_levels[1:8]
  n_sets <- 8
}

has_subgenomes <- (n_sets >= 2)

# Okabe-Ito colorblind-friendly palette for subgenomes (up to 8)
# Use first 6 colors for subgenomes 1-6, then light blue (oi[2]) for 7, black (oi[8]) for 8
oi <- palette_okabe_ito()
# Subgenome order: orange, bluish green, blue, vermillion, reddish purple, yellow, sky blue, black
# Reserve sky blue and black for subgenomes 7-8 as requested
pal_dark <- c(
  oi[1],  # 1: orange
  oi[3],  # 2: bluish green
  oi[5],  # 3: blue
  oi[6],  # 4: vermillion
  oi[7],  # 5: reddish purple
  oi[4],  # 6: yellow
  oi[2],  # 7: sky blue (reserved)
  oi[8]   # 8: black (reserved)
)
# Light variants for reference lines and segments
pal_light <- lighten(pal_dark, amount = 0.5)

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

# For vertical orientation: chromosomes on y-axis, chr1 at top (highest y)
y_levels <- rev(chrom_levels)  # chr20, chr19, ..., chr1
y_tbl <- tibble(chrom_id = factor(y_levels, levels = y_levels),
                y_index = seq_along(y_levels) * chrom_spacing)  # chr20=lowest, chr1=highest

ref_lines_base <- ref_sg %>%
  mutate(chrom_id = factor(chrom_id, levels = y_levels)) %>%
  left_join(y_tbl, by = "chrom_id") %>%
  filter(!is.na(y_index)) %>%
  mutate(ref_len_mb = ref_len / 1e6)

df_plot <- df %>%
  mutate(
    contig_len = if_else(is.na(length), 0, as.numeric(length)),
    assigned_subgenome = case_when(
      !is.na(assigned_subgenome) & assigned_subgenome %in% plot_levels ~ as.character(assigned_subgenome),
      TRUE ~ "NA"
    ),
    assigned_chrom_id  = if_else(assigned_chrom_id %in% chrom_levels, assigned_chrom_id, "Un"),
    chrom_id = factor(assigned_chrom_id, levels = y_levels),
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
# Left: Contig composition from segments (VERTICAL orientation)
# ----------------------------
df_slots <- df_plot %>%
  mutate(contig_len_mb = contig_len / 1e6) %>%
  filter(assigned_chrom_id != "Un") %>%
  left_join(y_tbl, by = "chrom_id") %>%
  filter(!is.na(y_index)) %>%
  group_by(chrom_id) %>%
  arrange(assigned_subgenome, desc(contig_len), contig, .by_group = TRUE) %>%
  mutate(
    n_in_bin = n(),
    slot = row_number(),
    total_width = 0.80,
    dy = if_else(n_in_bin > 1, total_width / n_in_bin, 0),
    y_offset = if_else(n_in_bin > 1, (slot - (n_in_bin + 1) / 2) * dy, 0),
    y_plot = y_index + y_offset,
    # Use assigned_subgenome from df_plot (already extracted from summary file)
    # Must be a factor with levels matching col_dark for scale_color_manual to work
    ref_subgenome = factor(assigned_subgenome, levels = names(col_dark))
  ) %>%
  ungroup()

n_per_chrom <- df_slots %>% count(chrom_id, name = "n_contigs")
max_n <- max(n_per_chrom$n_contigs, na.rm = TRUE)
if (!is.finite(max_n) || max_n <= 0) max_n <- 1
band_yw <- max(0.02, min(0.08, 0.80 / max_n * 0.80))

# Width for pill backgrounds and alignment segments (same width for both)
pill_yw <- band_yw * 2.4

# ----------------------------
# Rearrangement hypothesis labels
# ----------------------------
# Helper to extract short chromosome IDs (e.g., "chr2A,chr5B" -> "2A,5B")
extract_short_chrom_ids <- function(candidates) {
  if (is.na(candidates) || candidates == "") return(NA_character_)
  ids <- str_split(candidates, ",")[[1]]
  short_ids <- str_replace(ids, "^[Cc]hr", "")
  paste(short_ids, collapse = ",")
}

# Check if rearrangement_candidates column exists
has_rearr_col <- "rearrangement_candidates" %in% names(df_plot)

# Filter contigs with rearrangement candidates and create labels
if (has_rearr_col) {
  df_rearr <- df_plot %>%
    filter(!is.na(rearrangement_candidates) & rearrangement_candidates != "") %>%
    mutate(
      short_ids = sapply(rearrangement_candidates, extract_short_chrom_ids, USE.NAMES = FALSE),
      rearr_label = paste0("<-> ", short_ids),  # rearrangement indicator + chromosome IDs
      # Extract first candidate for coloring (most significant off-target)
      first_candidate = str_split_fixed(rearrangement_candidates, ",", 2)[, 1],
      # Extract subgenome from first candidate (last character if single letter A-Z)
      offtarget_subgenome = str_extract(first_candidate, "[A-Z]$"),
      # Extract chrom_id from first candidate (everything before the subgenome suffix)
      offtarget_chrom_id = str_replace(first_candidate, "[A-Z]$", ""),
      # Determine if homeologous exchange (same chrom_id, different subgenome)
      is_homeologous = (offtarget_chrom_id == assigned_chrom_id),
      # Set color: use off-target subgenome dark color if homeologous, red if not
      label_color = if_else(
        is_homeologous & !is.na(offtarget_subgenome) & offtarget_subgenome %in% names(sg_cols_dark),
        sg_cols_dark[offtarget_subgenome],
        "red"
      )
    )

  # Join rearrangement data with positioning info
  rearr_labels <- df_rearr %>%
    inner_join(
      df_slots %>% select(original_name, y_plot, contig_len_mb),
      by = "original_name"
    ) %>%
    mutate(label_x = contig_len_mb)

  # Calculate minimum x for labels (to right of longest contig)
  max_contig_length <- max(df_slots$contig_len_mb, na.rm = TRUE)
  rearr_label_xmin <- max_contig_length + 0.5  # Labels must be to right of all pills

  # Calculate x-axis expansion needed for labels (as fraction of max length)
  # Add ~15% extra space at right for labels
  x_expand_right <- 0.15
} else {
  rearr_labels <- tibble()
  rearr_label_xmin <- Inf
  x_expand_right <- 0.02  # Minimal expansion when no labels
}

# Reference lines: horizontal lines above topmost pill, one per subgenome with assigned contigs
# Basic data prep (y_ref calculation deferred until after pill_yw is defined)
ref_lines_data <- df_slots %>%
  group_by(chrom_id) %>%
  summarise(
    y_top = max(y_plot),  # topmost pill position in this bin
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
    df_slots %>% select(original_name, chrom_id, y_plot, contig_len_mb_slot = contig_len_mb, is_reversed),
    by = c("contig" = "original_name")
  ) %>%
  filter(!is.na(y_plot)) %>%
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
    ymin = y_plot - pill_yw/2,
    ymax = y_plot + pill_yw/2,
    xmin = qstart_mb,
    xmax = qend_mb
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
        y_plot,
        contig_len_mb_slot = contig_len_mb,
        is_reversed
      ),
    by = c("contig" = "original_name")
  ) %>%
  filter(!is.na(y_plot)) %>%
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
    ymin = y_plot - pill_yw/2,
    ymax = y_plot + pill_yw/2,
    xmin = qstart_mb,
    xmax = qend_mb
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

# Prepare rDNA annotations for plotting (overlay on pills)
rdna_plot <- if (nrow(rdna_annot) > 0) {
  rdna_annot %>%
    filter(contig %in% df_slots$original_name) %>%
    mutate(
      start_mb = start / 1e6,
      end_mb = end / 1e6
    ) %>%
    left_join(
      df_slots %>% select(original_name, y_plot, contig_len_mb, is_reversed),
      by = c("contig" = "original_name")
    ) %>%
    filter(!is.na(y_plot)) %>%
    mutate(
      # Flip coordinates for reversed contigs to show in reference orientation
      start_flip = if_else(is_reversed, contig_len_mb - end_mb, start_mb),
      end_flip = if_else(is_reversed, contig_len_mb - start_mb, end_mb),
      start_mb = pmax(0, pmin(start_flip, contig_len_mb)),
      end_mb = pmax(0, pmin(end_flip, contig_len_mb))
    ) %>%
    filter(end_mb > start_mb) %>%
    mutate(
      ymin = y_plot - pill_yw/2,
      ymax = y_plot + pill_yw/2,
      xmin = start_mb,
      xmax = end_mb,
      tooltip = paste0(
        "rDNA full copy\n",
        "contig: ", contig, "\n",
        "start: ", sprintf("%.3f", start_mb), " Mbp\n",
        "end: ", sprintf("%.3f", end_mb), " Mbp\n",
        "identity: ", sprintf("%.4f", identity), "\n",
        "coverage: ", sprintf("%.4f", consensus_coverage)
      )
    )
} else {
  tibble()  # Empty tibble
}

# Reference line parameters
ref_tick_offset <- pill_yw * 0.5  # gap between lines and pill
ref_line_width <- 0.4  # linewidth in mm
ref_line_spacing <- 0.04  # y-axis spacing between lines (data units)

# Calculate y positions for reference lines (now that pill_yw is defined)
ref_lines <- ref_lines_data %>%
  mutate(
    # Position from bottom to top (bottommost = closest to pill)
    y_ref = y_top + pill_yw/2 + ref_tick_offset + (sg_idx - 1) * ref_line_spacing * 2
  )

p_comp <- ggplot() +
  # Reference chromosome length lines (horizontal, above topmost pill per bin)
  geom_segment(
    data = ref_lines,
    aes(
      x = 0, xend = ref_len_mb,
      y = y_ref, yend = y_ref,
      color = assigned_subgenome
    ),
    linewidth = ref_line_width,
    alpha = 0.9,
    show.legend = FALSE
  ) +
  scale_color_manual(values = col_light, guide = "none", drop = FALSE) +
  ggnewscale::new_scale_color() +
  # Telomere indicators - slightly darker circles behind pill background
  # 5' telomere (left end of contig)
  geom_point(
    data = df_slots %>% filter(has_5p_telo),
    aes(x = 0, y = y_plot),
    color = "grey75",
    size = pill_yw * 7,
    show.legend = FALSE
  ) +
  # 3' telomere (right end of contig)
  geom_point(
    data = df_slots %>% filter(has_3p_telo),
    aes(x = contig_len_mb, y = y_plot),
    color = "grey75",
    size = pill_yw * 7,
    show.legend = FALSE
  ) +
  # Pill backgrounds - all light gray (same height as alignment blocks)
  geom_rect(
    data = df_slots,
    aes(xmin = 0, xmax = contig_len_mb,
        ymin = y_plot - pill_yw/2, ymax = y_plot + pill_yw/2),
    fill = "grey88",
    color = NA,
    show.legend = FALSE
  ) +
  ggnewscale::new_scale_fill() +
  # Alignment segments (full pill height)
  geom_rect(
    data = seg_plot %>% filter(!off_target),
    aes(xmin = xmin, xmax = xmax,
        ymin = ymin, ymax = ymax, fill = target_subgenome),
    color = NA,
    alpha = 1.0
  ) +
  geom_rect(
    data = seg_plot %>% filter(off_target),
    aes(xmin = xmin, xmax = xmax,
        ymin = ymin, ymax = ymax),
    fill = "red",
    color = NA,
    alpha = 1.0
  ) +
  scale_fill_manual(values = col_light, guide = "none", drop = FALSE) +
  ggnewscale::new_scale_fill() +
  # Macro blocks (full pill height)
  geom_rect(
    data = macro_plot %>% filter(!off_target),
    aes(xmin = xmin, xmax = xmax,
        ymin = ymin, ymax = ymax, fill = target_subgenome),
    color = NA,
    alpha = 1.0,
    show.legend = FALSE
  ) +
  geom_rect(
    data = macro_plot %>% filter(off_target),
    aes(xmin = xmin, xmax = xmax,
        ymin = ymin, ymax = ymax),
    fill = "red",
    color = NA,
    alpha = 1.0,
    show.legend = FALSE
  ) +
  # rDNA full-copy annotations (bold purple, on top of everything)
  {
    if (nrow(rdna_plot) > 0) {
      geom_rect(
        data = rdna_plot,
        aes(xmin = xmin, xmax = xmax,
            ymin = ymin, ymax = ymax),
        fill = "#7B2D8E",
        color = "#7B2D8E",
        linewidth = 0.05,
        alpha = 1.0,
        show.legend = FALSE
      )
    }
  } +
  # Rearrangement hypothesis labels (colored by type: subgenome color for homeologous, red for non-homeologous)
  {
    if (nrow(rearr_labels) > 0) {
      # Create list of layers, one per unique color (ggrepel doesn't support per-point segment colors)
      label_layers <- lapply(unique(rearr_labels$label_color), function(col) {
        ggrepel::geom_text_repel(
          data = rearr_labels %>% filter(label_color == col),
          aes(x = label_x, y = y_plot, label = rearr_label),
          color = col,
          size = (base_font_pt + 2) / .pt * 0.5,
          family = base_family,
          direction = "x",
          xlim = c(rearr_label_xmin, NA),  # Force labels to right of all pills
          box.padding = 0.3,
          point.padding = 0.5,
          segment.color = col,
          segment.size = 0.2,
          max.overlaps = 20,
          min.segment.length = 0,
          show.legend = FALSE
        )
      })
      label_layers
    }
  } +
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
    expand = expansion(mult = c(0.02, x_expand_right)),
    position = "top"
  ) +
  scale_y_continuous(
    breaks = y_tbl$y_index,
    labels = str_replace(as.character(y_tbl$chrom_id), "^chr", ""),
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
    x = "Mbp",
    y = "Assigned chromosome",
    title = paste0("Contig composition (", plot_title_suffix, ")")
  )

p_comp <- apply_legend_theme(p_comp, text_pt = 6, key_pt = 6, tight = TRUE) +
  theme(legend.margin = margin(1, 1, 1, 1))

# Interactive version (tooltips on segment/macro blocks)
if (plot_html) {
  p_comp_html <- ggplot() +
    # Reference chromosome length lines (horizontal, above topmost pill per bin)
    geom_segment(
      data = ref_lines,
      aes(
        x = 0, xend = ref_len_mb,
        y = y_ref, yend = y_ref,
        color = assigned_subgenome
      ),
      linewidth = ref_line_width,
      alpha = 0.9,
      show.legend = FALSE
    ) +
    scale_color_manual(values = col_light, guide = "none", drop = FALSE) +
    ggnewscale::new_scale_color() +
    # Telomere indicators - slightly darker circles behind pill background
    # 5' telomere (left end of contig)
    ggiraph::geom_point_interactive(
      data = df_slots %>% filter(has_5p_telo),
      aes(x = 0, y = y_plot,
          tooltip = paste0(original_name, "\n5' telomere detected")),
      color = "grey75",
      size = pill_yw * 5,
      show.legend = FALSE
    ) +
    # 3' telomere (right end of contig)
    ggiraph::geom_point_interactive(
      data = df_slots %>% filter(has_3p_telo),
      aes(x = contig_len_mb, y = y_plot,
          tooltip = paste0(original_name, "\n3' telomere detected")),
      color = "grey75",
      size = pill_yw * 5,
      show.legend = FALSE
    ) +
    # Pill backgrounds - all light gray (same height as alignment blocks)
    geom_rect(
      data = df_slots,
      aes(xmin = 0, xmax = contig_len_mb,
          ymin = y_plot - pill_yw/2, ymax = y_plot + pill_yw/2),
      fill = "grey88",
      color = NA,
      show.legend = FALSE
    ) +
    ggnewscale::new_scale_fill() +
    # Alignment segments (full pill height)
    ggiraph::geom_rect_interactive(
      data = seg_plot %>% filter(!off_target),
      aes(xmin = xmin, xmax = xmax,
          ymin = ymin, ymax = ymax, fill = target_subgenome, tooltip = tooltip),
      color = NA,
      alpha = 1.0
    ) +
    ggiraph::geom_rect_interactive(
      data = seg_plot %>% filter(off_target),
      aes(xmin = xmin, xmax = xmax,
          ymin = ymin, ymax = ymax, tooltip = tooltip),
      fill = "red",
      color = NA,
      alpha = 1.0
    ) +
    scale_fill_manual(values = col_light, guide = "none", drop = FALSE) +
    ggnewscale::new_scale_fill() +
    # Macro blocks (full pill height)
    ggiraph::geom_rect_interactive(
      data = macro_plot %>% filter(!off_target),
      aes(xmin = xmin, xmax = xmax,
          ymin = ymin, ymax = ymax, fill = target_subgenome, tooltip = tooltip),
      color = NA,
      alpha = 1.0,
      show.legend = FALSE
    ) +
    ggiraph::geom_rect_interactive(
      data = macro_plot %>% filter(off_target),
      aes(xmin = xmin, xmax = xmax,
          ymin = ymin, ymax = ymax, tooltip = tooltip),
      fill = "red",
      color = NA,
      alpha = 1.0,
      show.legend = FALSE
    ) +
    # rDNA full-copy annotations (bold purple, on top of everything, with tooltips)
    {
      if (nrow(rdna_plot) > 0) {
        ggiraph::geom_rect_interactive(
          data = rdna_plot,
          aes(xmin = xmin, xmax = xmax,
              ymin = ymin, ymax = ymax, tooltip = tooltip),
          fill = "#7B2D8E",
          color = NA,
          alpha = 0.85,
          show.legend = FALSE
        )
      }
    } +
    # Rearrangement hypothesis labels (colored by type: subgenome color for homeologous, red for non-homeologous)
    {
      if (nrow(rearr_labels) > 0) {
        # Create list of layers, one per unique color (ggrepel doesn't support per-point segment colors)
        label_layers <- lapply(unique(rearr_labels$label_color), function(col) {
          ggrepel::geom_text_repel(
            data = rearr_labels %>% filter(label_color == col),
            aes(x = label_x, y = y_plot, label = rearr_label),
            color = col,
            size = (base_font_pt + 2) / .pt * 0.5,
            family = base_family,
            direction = "x",
            xlim = c(rearr_label_xmin, NA),  # Force labels to right of all pills
            box.padding = 0.3,
            segment.color = col,
            segment.size = 0.2,
            max.overlaps = 20,
            min.segment.length = 0,
            show.legend = FALSE
          )
        })
        label_layers
      }
    } +
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
      expand = expansion(mult = c(0.02, x_expand_right)),
      position = "top"
    ) +
    scale_y_continuous(
      breaks = y_tbl$y_index,
      labels = str_replace(as.character(y_tbl$chrom_id), "^chr", ""),
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
      x = "Contig position (Mbp)",
      y = "Assigned chromosome",
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
    filter(assigned_subgenome %in% plot_levels) %>%
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
      shape = 16, alpha = 0.80, size = 1.6, show.legend = FALSE
    ) +
    scale_color_manual(values = col_dark, drop = FALSE) +
    coord_fixed(ratio = 1, xlim = radar_xlim, ylim = radar_ylim, clip = "off") +
    theme_classic(base_family = base_family, base_size = base_font_pt) +
    theme(
      axis.title.x = element_text(hjust = 0.5, size = base_font_pt, family = base_family, margin = margin(t = -2)),
      axis.title.y = element_blank(),
      axis.text  = element_blank(),
      axis.ticks = element_blank(),
      axis.line = element_blank(),
      panel.grid = element_blank(),
      plot.margin = margin(2, 2, 2, 2)
    ) +
    labs(x = "Chromosome-set support", y = NULL)

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

  # Calculate mean and SD per subgenome for error bars
  id_stats <- df_id_assigned %>%
    group_by(assigned_sg) %>%
    summarise(
      mean_id = mean(best_identity, na.rm = TRUE),
      sd_id = sd(best_identity, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      ymin = pmax(0, mean_id - sd_id),
      ymax = pmin(1, mean_id + sd_id)
    )

  p_id <- ggplot(df_id_assigned, aes(x = assigned_sg, y = pmax(0, pmin(1, best_identity)), color = assigned_sg)) +
    geom_jitter(width = 0.18, height = 0, alpha = 0.55, size = 1.2, show.legend = FALSE) +
    # Show mean ± 1 SD as error bars
    geom_errorbar(
      data = id_stats,
      aes(x = assigned_sg, y = mean_id, ymin = ymin, ymax = ymax, color = assigned_sg),
      width = 0.1, linewidth = 0.4, show.legend = FALSE,
      position = position_nudge(x = 0.35)
    ) +
    geom_point(
      data = id_stats,
      aes(x = assigned_sg, y = mean_id, color = assigned_sg),
      shape = 16, size = 1.5, show.legend = FALSE,
      position = position_nudge(x = 0.35)
    ) +
    scale_color_manual(values = col_dark, drop = FALSE) +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
    theme_classic(base_family = base_family, base_size = base_font_pt) +
    axis_theme +
    theme(
      axis.ticks = element_blank(),
      panel.grid = element_blank(),
      plot.margin = margin(5.5, 5.5, 5.5, 5.5)
    ) +
    labs(
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
      axis.ticks = element_blank(),
      panel.grid = element_blank(),
      plot.margin = margin(5.5, 5.5, 5.5, 5.5)
    ) +
    labs(
      x = NULL,
      y = "Alignment identity proxy"
    )
}

# Placeholder for future third panel on right side
p_placeholder <- ggplot() +
  theme_void(base_family = base_family, base_size = base_font_pt) +
  theme(
    plot.background = element_rect(fill = "grey95", color = NA),
    plot.margin = margin(5.5, 5.5, 5.5, 5.5)
  ) +
  labs(title = NULL)

if (has_subgenomes) {
  right_col <- p_id / p_radar / p_placeholder + plot_layout(heights = c(1, 1, 1))
  full_plot <- (p_comp | right_col) + plot_layout(widths = c(7, 3))
} else {
  full_plot <- (p_comp | p_id / p_placeholder) + plot_layout(widths = c(7, 3), heights = c(1, 1))
}

#message("macro rows: ", nrow(macro_plot), "  seg rows: ", nrow(seg_plot))
ggsave(out_pdf, plot = full_plot, width = 7.2, height = 7.2, units = "in", dpi = 300, device = cairo_pdf)

if (plot_html) {
  if (has_subgenomes) {
    right_col_html <- p_id / p_radar / p_placeholder + plot_layout(heights = c(1, 1, 1))
    full_plot_html <- (p_comp_html | right_col_html) + plot_layout(widths = c(7, 3))
  } else {
    full_plot_html <- (p_comp_html | p_id / p_placeholder) + plot_layout(widths = c(7, 3), heights = c(1, 1))
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
