# CLAUDE.md — thyroid-exposomics-2025

Project conventions and architecture. Match these when reading or editing code.

## R code style

**Hierarchical comment structure** (use these markers, with dotted numbering):
- `#*` major section — `#* 9: Validation Plots Adjustment`
- `#+` subsection — `#+ 9.2: Plots that failed after review`
- `#-` detail / sub-subsection — `#- 9.2.1: Methylparaben (CP2252)`
- `#_` individual item — `#_ Compound Name (ID)`
- `#!` important note
- Plain `#` for regular inline comments inside code blocks.

**Density:** code is compact — NO blank lines between sections/subsections, and code
follows immediately after a comment header (no blank line after it). Comments are
brief (compound names, short phrases, action descriptions).

```r
#* 9: Major Section
if (!isTRUE(config$skip_something)) {
#+ 9.2: Subsection
#- 9.2.1: Compound Name (ID)
variable <- function_call(args)
# regular inline comment
#- 9.2.2: Next Compound (ID)
more_code <- function_call(args)
}
```

**Functions:** roxygen2 docs on utility functions; prefer tidyverse; respect YAML
config flags (skip logic); separate concerns (e.g. plot generation vs PDF compilation).

## Architecture & pipeline

- **Config:** `All_Run/config_dynamic.yaml` loaded via
  `load_dynamic_config(computer = "auto", config_path = "All_Run/config_dynamic.yaml")`.
  `computer = "auto"` only picks the machine profile *inside* the yaml; it does not
  locate the file (that's `config_path`, resolved via `here::here()` from the repo root).
- **Run order:** `All_Run/run.R` sources numbered scripts sequentially (`00a`, `00b`,
  `00c`, `01` … `18`).
- **Utilities auto-load:** `R/Scripts/00b_setup.R` sources every `R/Utilities/**/*.R`
  via `purrr::walk(list.files(...), source)`. A new utility function is available to
  later scripts only after `00b` runs (re-source it in a warm session).
- **Packages:** renv (auto-restore). Run `renv::snapshot()` after adding packages.
- **Palettes:** `R/Utilities/Visualization/themes.R` — `carcinogen_colors`,
  `IARC_colors`, `variant_colors`. Reference these; don't hardcode hex.
- **Outputs:** `Outputs/Figures/`, `Outputs/Tables/`, `Outputs/Validation/`.
  OneDrive for final storage; local temp dirs for run-specific I/O.

## Supplement build — SOURCE vs OUTPUT (do not edit generated files)

`Supplementary/Components/` is assembled by `R/Scripts/18_construct_supplementary.R`,
which **reads the source sections, concatenates them into a generated file, then
compiles the PDF**:

- **SOURCES (edit these):** `Sections/cover_page.Rmd` (YAML header + TOC),
  `Sections/note1.tex` (Supplementary Note 1), `Sections/figures.Rmd`,
  `Sections/tables.tex`.
- **GENERATED (never edit — overwritten each build):** `supplementary_material.Rmd`
  (by script 18); `Tables/ST1.tex`, `Tables/ST2.tex`, `Sections/abbreviations.tex`
  (by script 17 — fix the `build_ST*()` / build functions instead); the PDF; `Build_Logs/*`.
- **Supplementary figures:** validation figures are cached grobs compiled 2-per-page by
  `compile_sf_sub_pdf()` in script 15; their numbering comes from the `figure_order`
  sheet of `validation.xlsx` (OneDrive). Figure S1 is the carcinogenicity flowchart,
  rendered by `render_carcinogen_flowchart()`.
- **Page numbers** in the TOC, figure-caption `pp. X-Y` ranges, and the
  `\ifthenelse{\value{page}>N}` "(Continued)" thresholds are **hardcoded**. Do NOT
  convert to `\pageref` auto-numbering (it breaks with `\includepdf` spans). When
  front-matter shifts pages, finalize the numbers from the **actual rendered PDF**.

## Conventions

- **Fonts:** supplement body text = Times New Roman; all figures = Arial.
- **Units:** `ppm` / `ppb` are **lowercase everywhere** (displayed text, captions,
  axis labels). Code identifiers (data-frame column names like `mean_tumor_PPB`) are
  exempt — they are not displayed.
- **"variant" vs "type":** the code uses **"variant"** for the three DTC types (PTC,
  FTC, IEFVPTC) in all identifiers/filenames; the **manuscript and all displayed text
  use "type"** (per 2022 WHO classification — "variant" is reserved for genetic
  variants). Keep code identifiers as "variant" (renaming risks breaking the pipeline);
  use "type" in any displayed captions/prose. Disclaimer lives at the top of
  `00c_FTs.R` (first usage) and in the README.
