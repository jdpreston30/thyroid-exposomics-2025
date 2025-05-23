# Thyroid Exposomics Analysis
This repository contains code and data used for the analysis presented in the submitted manuscript 'Environmental Chemical Burden in Differentiated Thyroid Cancer' by Preston et al. 2025.

## Requirements

- R version >= 4.3.1
- Packages: broom, dplyr, forcats, ggplot2, multcompView, purrr, readxl, stringr, tibble, tidyr, writexl

## Files

- `source_code.R`: Main analysis script. All data required to run the code, excluding data required for the blocks under `#+ Demographics (Table 1)`, can be imported from the "Data and Metadata Files" folder or loaded into the global environment from the compiled R serialized object (`processed_data.rds`). All required packages can be found in the header `#* Dependencies`.

- `processed_data.rds`: R serialized object containing all objects and the full dataset, excluding demographic data.

### Data and Metadata Files

Folder containing all primary data and metadata files. These data can be directly imported into the analysis script.

- `primary_data.xlsx`: Contains the primary data for this project.
  - `lib.subject.summary`: Feature table for all annotated features.
  - `tumors_sequence`: Sequence file and associated sample metadata.
  - `lib.subject.qusummary`: Quantitative feature table for tumors.
  - `lib.subject.qsummary.cadaver`: Quantitative feature table for cadaver control thyroids.
  - `library`: Information on the full library of standards used.
  - `tissue_weights`: Weights of each tissue sample.

- `chemical_metadata.xlsx`: Metadata for all detected chemicals.
  - `feature_metadata`: Detailed metadata for all unique identifications.
  - `Endogenous Excluded Features`: Features excluded from analysis due to being primarily endogenous.

- `Outputs`: Folder with all output data files (`.csv`, `.xlsx`, etc.).

- `Figures.prism`: GraphPad Prism file used to create all manuscript figures.