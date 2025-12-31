# ST4 (In Progress) - Integration Instructions

## Overview
This file contains all archived ST4 (Supplementary Table 4) information for future reintegration. ST4 was temporarily removed from the supplementary material construction but can be easily reactivated in a few months when ready.

## ST4 Details

### Table Title
**SUPPLEMENTARY TABLE 4**

### Table Caption
**Quantitative estimates for all chemicals detected in both tumor and non-neoplastic thyroid tissue.** This table lists the estimated PPB values for all chemicals annotated in both tumor samples and non-neoplastic thyroid tissue. The estimated PPB values are calculated by taking the mean area for each chemical and dividing it by the mean area from reference standards with known concentrations (0·47 ng/mL). Chemicals are sorted alphabetically.

---

## How to Reintegrate ST4

### Step 1: Prepare ST4 Table Data
When ready, generate the ST4 table using R script 16_supplementary_tables.R (section 11.4 is marked as "Pending"). Create a `build_ST4()` function following the same pattern as `build_ST1()`, `build_ST2()`, and `build_ST3()`. Export the final LaTeX to `Supplementary/Components/Tables/ST4.tex`.

### Step 2: Add ST4 to tables.tex
In `/Users/jdp2019/Desktop/thyroid-exposomics-2025/Supplementary/Components/Sections/tables.tex`, replace the ST4 placeholder section with the proper input statement:

**Location:** Lines 45-56 (currently contains `[INSERT ST4 HERE]`)

**Replace this:**
```tex
\markright{\textit{Supplementary Table 4}}

\begin{landscape}

\noindent{\fontsize{12}{14.4}\selectfont\textbf{SUPPLEMENTARY TABLE 4}}

\noindent{\fontsize{10}{12}\selectfont\textbf{Quantitative estimates for all chemicals detected in both tumor and non-neoplastic thyroid tissue.} This table lists the estimated PPB values for all chemicals annotated in both tumor samples and non-neoplastic thyroid tissue. The estimated PPB values are calculated by taking the mean area for each chemical and dividing it by the mean area from reference standards with known concentrations (0·47 ng/mL). Chemicals are sorted alphabetically.}

[INSERT ST4 HERE]

\end{landscape}
```

**With this:**
```tex
\markright{\textit{Supplementary Table 4}}

\begin{landscape}

\noindent{\fontsize{12}{14.4}\selectfont\textbf{SUPPLEMENTARY TABLE 4}}

\noindent{\fontsize{10}{12}\selectfont\textbf{Quantitative estimates for all chemicals detected in both tumor and non-neoplastic thyroid tissue.} This table lists the estimated PPB values for all chemicals annotated in both tumor samples and non-neoplastic thyroid tissue. The estimated PPB values are calculated by taking the mean area for each chemical and dividing it by the mean area from reference standards with known concentrations (0·47 ng/mL). Chemicals are sorted alphabetically.}

\input{../Tables/ST4.tex}

\end{landscape}
```

### Step 3: Update supplementary_material.Rmd
In `/Users/jdp2019/Desktop/thyroid-exposomics-2025/Supplementary/Components/supplementary_material.Rmd`, replace the ST4 placeholder section (lines ~470-480) with the same content from Step 2.

### Step 4: Update Table of Contents
The table of contents in supplementary_material.Rmd should automatically update to include ST4 once the sections above are restored (the TOC is typically generated from LaTeX section headers).

### Step 5: Ensure Dynamic Paths
Verify that script 17_construct_supplementary.R includes path substitution for ST4.tex. Add the following line if not present:

```r
tables_content <- gsub('../Tables/ST4.tex', file.path(tables_dir, "ST4.tex"), tables_content, fixed = TRUE)
```

### Step 6: Test Rendering
Run the supplementary material rendering pipeline to ensure ST4 appears correctly:
- Execute `R/Scripts/17_construct_supplementary.R`
- Verify the output PDF includes ST4 with proper formatting and table inclusion

---

## Notes
- ST4 requires significant development work on the underlying chemical detection and quantification algorithms
- The table structure follows the same gt → LaTeX conversion pipeline as ST1-ST3
- All table formatting specs (font sizes, column widths, borders) should match ST1-ST3 conventions
- Citation handling and special characters should follow existing patterns
