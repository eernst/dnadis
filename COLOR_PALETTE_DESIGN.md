# Color Palette Design for final_finalizer R Visualizations

This document describes the standardized color palette strategy across all R plot templates.

## Design Principles

All visualization templates now use **ggokabeito** (Okabe-Ito colorblind-friendly palette) for consistency and accessibility.

## Palette Overview

The Okabe-Ito palette provides 8 colors + gray (9 total):
1. Orange
2. Sky blue
3. Bluish green
4. Yellow
5. Blue
6. Vermillion (red-orange)
7. Reddish purple
8. Black
9. Gray

## Usage by Template

### 1. classification_summary_bar.tmpl.R

**Purpose**: Classification categories (9 categories)

**Mapping**:
- `chrom_assigned`: oi[5] (blue)
- `chrom_unassigned`: oi[7] (reddish purple)
- `organelle_complete`: oi[3] (bluish green)
- `organelle_debris`: oi[4] (yellow)
- `rdna`: oi[1] (orange)
- `contaminant`: oi[6] (vermillion)
- `chrom_debris`: oi[2] (sky blue)
- `debris`: oi[9] (gray)
- `unclassified`: oi[8] (black)

**Overflow strategy**: If categories exceed 9, use `colorspace::lighten()` to create lighter variants of base colors.

### 2. chromosome_overview.tmpl.R

**Purpose**: Reference subgenomes (up to 8)

**Mapping** (reserves slots 7-8 as requested):
- Subgenome 1: oi[1] (orange)
- Subgenome 2: oi[3] (bluish green)
- Subgenome 3: oi[5] (blue)
- Subgenome 4: oi[6] (vermillion)
- Subgenome 5: oi[7] (reddish purple)
- Subgenome 6: oi[4] (yellow)
- Subgenome 7: oi[2] (sky blue) **[RESERVED]**
- Subgenome 8: oi[8] (black) **[RESERVED]**

**Light variants**: Reference lines use `colorspace::lighten(pal_dark, amount = 0.5)` for lighter shades.

**Single genome fallback**: Uses oi[5] (blue) when no subgenomes detected.

### 3. depth_overview.tmpl.R

**Purpose**: 
- Classification colors (for depth by classification plot)
- Subgenome colors (for chromosome depth plot)

**Classification mapping**: Same as classification_summary_bar.tmpl.R (lines 104-119), with debris variants using `lighten()`:
- `chrom_assigned`: oi[5] (blue)
- `chrom_debris`: `lighten(oi[5], 0.5)` (light blue)
- `chrC`: oi[3] (bluish green)
- `chrC_debris`: `lighten(oi[3], 0.5)` (light bluish green)
- `chrM`: oi[6] (vermillion)
- `chrM_debris`: `lighten(oi[6], 0.5)` (light vermillion)
- etc.

**Subgenome mapping**: Same as chromosome_overview.tmpl.R

**Breadth coverage colors**: Fixed 2-color palette (not from Okabe-Ito):
- `>=1x`: #66C2A5 (teal)
- `>=10x`: #1F77B4 (matplotlib blue)

## Design Rationale

### Why Okabe-Ito?
- **Accessibility**: Designed to be distinguishable by people with color vision deficiencies
- **Print-friendly**: Works well in both color and grayscale
- **Scientific standard**: Widely used in academic and scientific publications

### Why reserve slots 7-8 for subgenomes?
- Light blue (oi[2]) and black (oi[8]) are visually distinct from the first 6 colors
- Provides clear visual separation for higher ploidy levels (hexaploid, octoploid)
- Black provides strong contrast for the highest subgenome level

### Why lighten() for debris categories?
- Creates visual hierarchy: parent category (dark) vs. debris (light)
- Maintains color family relationship while providing distinction
- Avoids exhausting the 9-color palette

## Color Selection Guidelines

When adding new visualization categories:

1. **First choice**: Use unused Okabe-Ito colors (oi[1] through oi[9])
2. **Second choice**: Use `colorspace::lighten()` or `colorspace::darken()` variants
3. **Third choice**: Use `colorspace::desaturate()` for muted variants
4. **Avoid**: Creating entirely new color palettes that break consistency

## Dependencies

All templates now require:
- `ggokabeito`: Provides `palette_okabe_ito()` function
- `colorspace`: Provides `lighten()`, `darken()`, `desaturate()` functions

These are added to the `pacman::p_load()` calls in each template.

## Testing Recommendations

When modifying color palettes:
1. View plots on screen (color)
2. Print plots in grayscale to verify distinction
3. Use color blindness simulators (e.g., https://www.color-blindness.com/coblis-color-blindness-simulator/)
4. Verify legend clarity and labels

## Future Considerations

If classification categories or subgenome counts exceed palette capacity:
- **Option 1**: Use pattern fills (requires ggpattern package)
- **Option 2**: Split into multiple sub-plots
- **Option 3**: Use interactive filtering (already supported via ggiraph HTML output)
