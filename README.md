# Environmental Chemical Burden in Differentiated Thyroid Cancer

## 📖 Citation

This code is associated with the analysis presented in the following manuscript:
> Preston et al. (2025). Environmental Chemical Burden in Differentiated Thyroid Cancer. (submitted).

## 🚀 Quick Start for Reproduction

**⚠️ Data Availability Notice**: 
- **No raw data files** (~77 GB of GC-MS .raw files) are included in this repository
- **All instructions below assume you have obtained data files or are using your own data**
- **To reproduce this analysis**: Contact the first author (Joshua D. Preston, joshua.preston@emory.edu) or senior author (M. Ryan Smith, matthew.ryan.smith@emory.edu) to obtain the data files—this is the easiest and recommended approach
- **Public data access**: Raw GC-MS data will be made publicly available upon manuscript acceptance
- **To run analyses with your own data or provided data files**: Update file paths in `All_Run/config_dynamic.yaml` to match your system

### Option 1: Using Docker (Recommended for Exact Reproducibility)

**Status**: In Progress

### Option 2: Manual Installation (Without Docker)

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
- Installs ~30 R packages at exact versions from `renv.lock`
- Installs CRAN packages (e.g., ggplot2, dplyr, broom, tidyr)
- Installs Bioconductor packages (e.g., mzR for mass spectrometry data)
- Creates isolated project library (doesn't affect your system R packages)
- Only needed once per computer; subsequent runs use installed packages
- Packages are automatically loaded from `DESCRIPTION` file during pipeline execution

## 📁 Project Structure

```
├── DESCRIPTION                 # R package dependencies
├── renv.lock                   # Exact package versions for reproducibility
├── All_Run/                    # Pipeline execution
│   ├── config_dynamic.yaml     # Analysis configuration (update paths for your system)
│   ├── run.R                   # Main pipeline execution script
│   └── TO DO                   # Development notes
├── R/                          # Analysis code
│   ├── Scripts/                # Analysis workflow scripts (00a-17)
│   │   ├── 00a_environment_setup.R
│   │   ├── 00b_setup.R
│   │   ├── 00c_clinical_data.R
│   │   ├── 00d_FTs.R
│   │   ├── 00e_peakwalk_compile.R
│   │   ├── 01_demographics.R
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
│   │   ├── 13_render_figures.R
│   │   ├── 14_render_supplementary_figures.R
│   │   ├── 15_tables.R
│   │   ├── 16_supplementary_tables.R
│   │   └── 17_construct_supplementary.R
│   └── Utilities/              # Custom analysis functions
│       ├── Analysis/           # Statistical and carcinogen classification
│       ├── Helpers/            # Helper functions (config, validation, tables)
│       ├── Tabulation/         # Table generation (demographics, supplementary)
│       ├── Validation/         # Spectral validation and fragment processing
│       └── Visualization/      # Plotting functions (balloons, heatmaps, donuts)
├── Outputs/                    # Generated results
│   ├── Figures/                # Publication figures (PNG, PDF)
│   ├── Tables/                 # Manuscript tables
│   └── Validation/             # Spectral validation plots and PDFs
│       ├── failed/             # Compounds that failed validation
│       ├── initial_compile/    # Initial validation compilation
│       ├── revised/            # Revised validation plots
│       └── top_fragments/      # Top fragment validations
├── Supplementary/              # Supplementary materials
│   ├── Components/             # R Markdown components
│   ├── Build_Logs/             # LaTeX build logs
│   └── ST4_ip.md               # Supplementary Table 4 working document
├── metadata_files/             # Chemical metadata, libraries, tissue weights
└── utilities_hold/             # Archived unused functions (for testing)
```

## 🔬 Analysis Workflow

The complete pipeline executes in sequence:

1. **00a-00e**: Environment setup, clinical metadata, feature tables, peakwalk compilation
2. **01**: Demographics and clinical characteristics (Table 1)
3. **02**: Detection frequency analysis
4. **03**: Chemical class distribution
5. **04**: Variant-specific statistical comparisons
6. **05**: Variant visualization data preparation
7. **06**: Tumor vs cadaver control comparisons
8. **07-10**: Spectral validation workflow (preparation, execution, plotting, cleanup)
9. **11-12**: Variant and IARC carcinogen visualizations
10. **13-14**: Render main and supplementary figures
11. **15-16**: Generate manuscript and supplementary tables
12. **17**: Construct supplementary materials document

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

*Note: System dependencies will be automatically installed in the Docker container (In Progress). For manual installation, see above.*

## 📦 Package Dependencies

All R package dependencies are specified in `DESCRIPTION`. Key packages include:

### CRAN Packages
- **Data manipulation**: tidyverse (dplyr, tidyr, purrr, readr, stringr, tibble, forcats)
- **Visualization**: ggplot2, ggtext, cowplot, gridExtra, magick
- **Mass spectrometry**: mzR (Bioconductor)
- **Chemical informatics**: webchem
- **Statistical analysis**: broom, multcompView
- **Parallel processing**: doParallel, foreach
- **Document generation**: rmarkdown, tinytex, gt, openxlsx
- **Configuration**: yaml, here, jsonlite

### Bioconductor Packages
- **mzR**: Mass spectrometry data import and processing

*See `DESCRIPTION` file for complete list of all dependencies.*

## 🔄 Reproducibility Features

This project implements best practices for computational reproducibility:

- ✅ **Version Control**: Complete analysis code on GitHub
- ✅ **Package Management**: `renv` with `renv.lock` pinning all packages to exact versions
- ✅ **Dependency Declaration**: All dependencies specified in `DESCRIPTION` with automatic loading
- ✅ **Containerization**: In Progress
- ✅ **Docker Hub Distribution**: In Progress
- ✅ **Configuration-Driven**: All parameters in `config_dynamic.yaml` (computer-specific paths)
- ✅ **Dynamic Path Resolution**: Automatic detection of computer/user for path configuration
- ✅ **Documentation**: Comprehensive function documentation (roxygen2 style) and workflow comments
- ✅ **Hierarchical Code Organization**: Clear comment structure (#*, #+, #-, #_) for workflow navigation
- ✅ **Modular Design**: Utilities separated by function type (Analysis, Visualization, Validation, etc.)

## 📧 Contact

**First Author & Repository Maintainer**: Joshua D. Preston
- **Email**: joshua.preston@emory.edu  
- **ORCID**: [0000-0001-9834-3017](https://orcid.org/0000-0001-9834-3017)  
- **Institution**: Department of Surgery, Emory University School of Medicine

**Senior & Corresponding Author**: TBD
- **Email**: TBD
- **Institution**: Department of Surgery, Emory University School of Medicine

---

**Repository**: https://github.com/jdpreston30/thyroid-exposomics-2025  
**Docker Hub**: In Progress  
**Zenodo Archive**: In Progress
