# Environmental Chemical Burden in Differentiated Thyroid Cancer

## 📖 Citation

This code is associated with the analysis presented in the following manuscript:
> Preston et al. (2026). Environmental Chemical Burden in Differentiated Thyroid Cancer.
> *Environmental International* (submitted).

## 📝 Terminology Note (code vs. manuscript)

Throughout this codebase, the three differentiated thyroid cancer **types** analyzed — papillary (PTC), follicular (FTC), and the invasive encapsulated follicular variant of papillary thyroid carcinoma (IEFVPTC) — are referred to internally as **"variant"** (in column names, object names, function names, script names, and file paths). This is an intentional, isolated inconsistency with the manuscript.

Per the 2022 WHO Classification of Thyroid Tumours (Jung et al., 2022; WHO Classification of Tumours Editorial Board, 2022), the term "variant" is now reserved for genetic variants, and IEFVPTC is recognized as a distinct entity rather than a subtype of PTC. The manuscript accordingly uses **"type" / "tumor type"** throughout. We deliberately retain "variant" as an internal identifier to preserve the integrity and reproducibility of a validated, working pipeline — renaming would risk introducing errors into analysis code that is otherwise verified and stable. **In every case, "variant" in the code is equivalent to "type" in the manuscript and current pathology literature**; the discrepancy is purely nominal and confined to code-level naming.

## 🚀 Quick Start for Reproduction

**⚠️ Data Availability Notice**: 
- **No raw data files** are included in this repository
- **All instructions below assume you have obtained data files or are using your own data**
- **To reproduce this analysis**: Contact the first author (Joshua D. Preston, joshua.preston@emory.edu) or senior author (M. Ryan Smith, matthew.ryan.smith@emory.edu) to obtain the data files—this is the easiest and recommended approach
- **Public data access**: Open-format spectra (`.mzML`) are deposited at Metabolomics Workbench — link TBD, added on publication. Vendor `.raw` files are **not** deposited (see below)
- **To run analyses with your own data or provided data files**: Update file paths in `All_Run/config_dynamic.yaml` to match your system

## 📤 Metabolomics Workbench deposition

Open-format spectra (`.mzML`) for all 191 acquisitions across both batches are deposited at the NIH
Common Fund's National Metabolomics Data Repository, **Metabolomics Workbench**.

**Deposit link: TBD** — Study ID, Project ID and DOI will be added here on publication.

Three things about that deposit are worth knowing before you try to reproduce anything from it.

**Deposited filenames are batch-prefixed.** Every file is named `GC080_*` (DTC tumor batch,
Dec 2022) or `GC097_*` (non-cancer cadaver batch, Aug 2024). The prefix exists because 26
QC/blank/standard filenames are byte-identical across the two runs — `wash-1_1`, `BP1_1`,
`Qstd-1_1` and so on — so a single flat archive would otherwise collide on 52 files. Prefixing every
file, rather than only the colliding ones, keeps the archive self-describing: batch is readable from
any filename.

**The original acquisition name is published for every file** as `acquisition_file` in the
`SUBJECT_SAMPLE_FACTORS` block, alongside `study` (`tumor`/`cadaver`) and `analytical_batch`. That is
the name this pipeline uses internally, and it is what the spectral-validation panels in Supplementary
Figures S2-S4 print as `Sample:`. A published panel therefore maps directly onto a deposited file.

**The validation figures were not generated from the deposited mzML.** `convert_raw_to_mzml.R`
produces a separate conversion into `mzML_validation/` from the vendor `.raw` files, and the figure
code reads only from there. The `.raw` files are not deposited — they carry biorepository specimen
accessions in their `Sample name` field, which cannot be removed without corrupting the format, and
those identifiers are not publicly shareable. Reproducing the validation figures from the deposit
means substituting the deposited mzML, which derive from the same acquisitions but a different
conversion run: equivalent inputs, not byte-identical ones.

### Installation

**Prerequisites**: 
- R >= 4.5.1
- Git (to clone repository)

**Note**: This project uses `renv` for package management to ensure reproducibility. The `renv.lock` file contains exact versions of all packages used in the manuscript.

```r
# 1. Clone the repository
# (from terminal)
git clone https://github.com/jdpreston30/thyroid-exposomics-2025.git
cd thyroid-exposomics-2025

# 2. Start R in the project directory
# (renv automatically activates via .Rprofile)

# 3. Restore all packages at exact versions (first time only, ~10-20 minutes)
renv::restore()

# 4. Update configuration paths
# Edit All_Run/config_dynamic.yaml to set paths for your system:
#   - computers: Define your computer's user_home and onedrive_path
#   - paths.base_data_path: Path to GC-MS raw data parent directory
#   - All other paths use dynamic templates that auto-populate from these base settings

# 5. Run the complete analysis pipeline
source("All_Run/run.R")
```

**What happens during `renv::restore()`**:
- Installs all R packages at exact versions from `renv.lock`
- Installs CRAN packages (e.g., ggplot2, dplyr, broom, tidyr)
- Installs Bioconductor packages (e.g., mzR for mass spectrometry data)
- Creates isolated project library (doesn't affect your system R packages)
- Only needed once per computer; subsequent runs use installed packages
- Packages are automatically loaded from `DESCRIPTION` file during pipeline execution

## 📁 Project Structure

```
├── DESCRIPTION                 # R package dependencies
├── renv.lock                   # Exact package versions for reproducibility
├── SESSION_INFO.txt            # Session record from the manuscript run
├── All_Run/                    # Pipeline execution
│   ├── config_dynamic.yaml     # Analysis configuration (update paths for your system)
│   └── run.R                   # Main pipeline execution script
├── R/                          # Analysis code
│   ├── Scripts/                # Analysis workflow scripts (00a-18)
│   │   ├── 00a_environment_setup.R
│   │   ├── 00b_setup.R
│   │   ├── 00c_FTs.R
│   │   ├── 00d_peakwalk_compile.R
│   │   ├── 01_clinical_data.R
│   │   ├── 02_detection.R
│   │   ├── 03_classes.R
│   │   ├── 04_variant_stats.R
│   │   ├── 05_variant_vis_prep.R
│   │   ├── 06_tumor_cadaver.R
│   │   ├── 07_validation_prep.R
│   │   ├── 08_validation_run.R
│   │   ├── 09_validation_plots_create.R
│   │   ├── 10_post_validation_clean.R
│   │   ├── 11_variant_vis.R
│   │   ├── 12_IARC_vis.R
│   │   ├── 13_confounding_analysis.R
│   │   ├── 14_render_figures.R
│   │   ├── 15_render_supplementary_figures.R
│   │   ├── 16_tables.R
│   │   ├── 17_supplementary_tables.R
│   │   └── 18_construct_supplementary.R
│   └── Utilities/              # Custom analysis functions
│       ├── Analysis/           # Statistical and carcinogen classification
│       ├── Clinical/           # AJCC 8th ed. staging and T-category assignment
│       ├── Helpers/            # Helper functions (config, validation, tables)
│       ├── Tabulation/         # Table generation (demographics, supplementary)
│       ├── Terminal/           # Terminal helper functions
│       ├── Validation/         # Spectral validation and fragment processing
│       └── Visualization/      # Plotting functions (balloons, heatmaps, donuts)
├── Outputs/                    # Generated results
│   ├── Figures/                # Publication figures (PNG, PDF)
│   ├── Tables/                 # Manuscript tables (Excel format)
│   └── Validation/             # Spectral validation plots and PDFs
│       ├── failed/             # Compounds that failed validation
│       ├── initial_compile/    # Initial validation compilation
│       ├── revised/            # Revised validation plots
│       └── top_fragments/      # Top fragment validations
└── Supplementary/              # Materials for compiled supplementary PDF
    ├── Components/             # R Markdown components
    └── Build_Logs/             # LaTeX build logs
```

Chemical metadata and reference libraries live outside the repository (OneDrive); paths are set in
`All_Run/config_dynamic.yaml`.

## 🔬 Analysis Workflow

The complete pipeline executes in sequence:

1. **00a-00d**: Environment setup, feature tables, peakwalk compilation
2. **01**: Clinical data, demographics, and Table 1
3. **02**: Detection frequency analysis
4. **03**: Chemical class distribution
5. **04**: Variant-specific statistical comparisons
6. **05**: Variant visualization data preparation
7. **06**: Tumor vs cadaver control comparisons
8. **07-10**: Spectral validation workflow (preparation, execution, plotting, cleanup)
9. **11-12**: Variant and IARC carcinogen visualizations
10. **13**: Confounding analysis (age, sex, collection timing)
11. **14-15**: Render main and supplementary figures
12. **16-17**: Generate manuscript and supplementary tables
13. **18**: Construct supplementary materials document

## 💻 System Requirements

### Computational Requirements
- **R**: Version 4.5.1 or higher
- **Platform**: Developed on macOS (M1/Apple Silicon) but cross-platform compatible
- **Memory**: Minimum 8 GB RAM recommended for large GC-MS datasets
- **Storage**: ~100 GB for raw data + processed outputs

### System Dependencies
- **TinyTeX/LaTeX**: PDF generation (automatically installed via tinytex package)
- **Mono framework**: Required for ThermoRawFileParser (.raw file conversion to mzML)
- **ThermoRawFileParser**: Converts Thermo .raw files to open mzML format
  - Installation: `~/bin/ThermoRawFileParser/`
  - Download: https://github.com/compomics/ThermoRawFileParser


## 📦 Package Dependencies

All CRAN, Bioconductor and GitHub dependencies are declared in `DESCRIPTION` and pinned in
`renv.lock`. See those files for the complete, authoritative list.

## 🔄 Reproducibility Features

This project implements best practices for computational reproducibility:

- ✅ **Version Control**: Complete analysis code on GitHub
- ✅ **Package Management**: `renv` with `renv.lock` pinning all packages to exact versions
- ✅ **Dependency Declaration**: All dependencies specified in `DESCRIPTION` with automatic loading
- ✅ **Configuration-Driven**: All parameters in `config_dynamic.yaml` (computer-specific paths)
- ✅ **Dynamic Path Resolution**: Automatic detection of computer/user for path configuration
- ✅ **Documentation**: Comprehensive function documentation (roxygen2 style) and workflow comments
- ✅ **Hierarchical Code Organization**: Clear comment structure (#*, #+, #-, #_) for workflow navigation
- ✅ **Modular Design**: Utilities separated by function type (Analysis, Visualization, Validation, etc.)

## 📧 Contact

**First Author & Repository Maintainer**: Joshua D. Preston
- **Email**: joshua.preston@emory.edu  
- **ORCID**: [0000-0001-9834-3017](https://orcid.org/0000-0001-9834-3017)  
- **Institution**: Emory University School of Medicine

**Senior & Corresponding Author**: M. Ryan Smith
- **Email**: matthew.ryan.smith@emory.edu
- **ORCID**: [0000-0002-8889-3477](https://orcid.org/0000-0002-8889-3477)  
- **Institution**: Emory University School of Medicine

---

**Repository**: https://github.com/jdpreston30/thyroid-exposomics-2025  
**Metabolomics Workbench**: TBD
