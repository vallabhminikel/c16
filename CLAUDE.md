# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

This is a biomedical research manuscript project analyzing C16 siRNA therapeutics targeting PRNP (prion protein) in animal models. The goal is to generate publication-ready figures, supplementary tables, and in-text statistics.

## Running the Analysis

There is no build system. Run R scripts directly:

```r
# Step 1: Prepare data (reads from ../a/ parent project, writes to data/)
Rscript src/prep_c16_data.R

# Step 2: Generate all figures, supplement, and stats_for_text.txt
Rscript src/c16_manuscript.R
```

Both scripts use `setwd()` internally and assume the working directory is the repo root, or they set it explicitly. Run from the repo root or from within R/RStudio.

## Dependencies

R packages required: `tidyverse`, `survival`, `janitor`, `magick`, `smoother`, `openxlsx` (implied by supplement generation). Install via `install.packages(c("tidyverse", "survival", "janitor", "magick", "smoother", "openxlsx"))`.

## Architecture

**Two-script pipeline:**

1. **`src/prep_c16_data.R`** — Reads raw experimental data from a sibling project directory (`../a/`) and writes 10 analysis-ready TSV files to `data/`. Run this when upstream raw data changes.

2. **`src/c16_manuscript.R`** — Reads from `data/`, performs all analyses, and writes to `display_items/`:
   - `figure-*.png` — Publication figures
   - `supplement.xlsx` — Supplementary data workbook
   - `stats_for_text.txt` — Statistics for embedding in manuscript narrative
   - `data-s*.tsv` — TSV exports of supplement tables

**Data flow:** `../a/` raw data → `prep_c16_data.R` → `data/*.tsv` → `c16_manuscript.R` → `display_items/`

## Key Concepts in the Code

**Cohort system:** Animals are grouped into cohorts defined in `data/meta.tsv`, which controls display names (`disp`), plot colors (`color` as hex), line types (`lty`), and whether to exclude from analysis. Every plot joins against this metadata.

**Five experiments:** A, YA, AF, AR, QS — each comparing siRNA-treated animals to aCSF controls, varying treatment timing, dose, or genetic background.

**Custom functions in c16_manuscript.R:**
- `spearman_karber()` — Calculates SD50 titers from dilution series (RT-QuIC infectivity assay)
- `make_sf()` / `plot_sf()` — Kaplan-Meier survival curve wrappers with cohort-consistent styling
- `norm_wts()` / `plot_wts()` — Normalize and plot longitudinal weight data (relative to baseline)
- `norm_nests()` — Aggregate cage-level nesting behavior scores
- `alpha()` — Adds transparency to hex color strings (used for CI shading in plots)

**Output management:** `c16_manuscript.R` opens a sink to `stats_for_text.txt` and an `openxlsx` workbook at startup; both are finalized at the end of the script. Figures are written with `png()`/`dev.off()` blocks.
