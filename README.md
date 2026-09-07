# Environmental Chemical Burden in Differentiated Thyroid Cancer

Analysis code for a GC–MS environmental chemical screening study of 60 differentiated thyroid
cancer (DTC) specimens — papillary, follicular, and invasive encapsulated follicular variant of
papillary thyroid carcinoma (n = 20 each) — compared with 8 non-cancer cadaver thyroids.

## Citation

> Preston JD, Liang Y, Szabo Yamashita T, Teeny S, Weinberg J, Crandall WJ, Jarrell ZR, Hu X,
> Safley SA, Robertson JM, Tran V, Jackson AS, Patel SG, Glosser LD, Weber CJ, Sharma J,
> Saunders ND, Go Y, Jones DP, Smith MR. Environmental Chemical Burden in Differentiated Thyroid
> Cancer. *Environment International*. Under revision, 2026.

## Terminology: "variant" in the code, "type" in the manuscript

Throughout this codebase, the three differentiated thyroid cancer **types** analyzed — papillary
(PTC), follicular (FTC), and the invasive encapsulated follicular variant of papillary thyroid
carcinoma (IEFVPTC) — are referred to internally as **"variant"**, in column names, object names,
function names, script names and file paths. This is an intentional, isolated inconsistency with the
manuscript.

Per the 2022 WHO Classification of Thyroid Tumours (Jung et al., 2022; WHO Classification of Tumours
Editorial Board, 2022), "variant" is now reserved for genetic variants, and IEFVPTC is recognized as
a distinct entity rather than a subtype of PTC. The manuscript accordingly uses **"type" / "tumor
type"** throughout. We retain "variant" as an internal identifier to preserve the integrity of a
validated, working pipeline — renaming would risk introducing errors into analysis code that is
otherwise verified and stable. **In every case, "variant" in the code is equivalent to "type" in the
manuscript and the current pathology literature**; the discrepancy is purely nominal and confined to
code-level naming.

## Data availability

**No data files are included in this repository.** Demographic and clinical data on study
participants are withheld.

Open-format spectra (`.mzML`) and the processed feature tables are deposited at Metabolomics
Workbench, [Study ST005162](https://doi.org/10.21228/M87P2Q). Vendor `.raw` files are **not**
deposited — see below.

To reproduce the analysis, the simplest route is to contact the first author (Joshua D. Preston,
joshua.preston@emory.edu) or the senior author (M. Ryan Smith, matthew.ryan.smith@emory.edu) for the
data files. To run the pipeline against those files or your own, edit the paths in
`All_Run/config_dynamic.yaml` to match your system.

## Metabolomics Workbench deposition

Open-format spectra (`.mzML`) for all 191 acquisitions across both batches, together with the
processed feature tables, are deposited at the NIH Common Fund's National Metabolomics Data
Repository, **Metabolomics Workbench**.

| | |
|---|---|
| Project ID | PR003330 |
| Study ID | ST005162 |
| Project DOI | https://doi.org/10.21228/M87P2Q |

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

## Requirements

- **R ≥ 4.5.1**
- **Memory**: 8 GB RAM minimum for the GC–MS feature tables
- **Storage**: ~100 GB for raw data plus processed outputs
- **Platform**: developed on macOS (Apple Silicon); cross-platform compatible

System tools:

- **TinyTeX / LaTeX** — supplementary PDF generation, installed automatically via `tinytex`
- **Mono** — required by ThermoRawFileParser
- **[ThermoRawFileParser](https://github.com/compomics/ThermoRawFileParser)** — converts Thermo
  `.raw` files to open `.mzML`; expected at `~/bin/ThermoRawFileParser/`

Package versions are pinned with **renv** (`renv.lock`); CRAN, Bioconductor and GitHub dependencies
are declared in `DESCRIPTION`. Those two files are the authoritative list.

## Running the analysis

```bash
git clone https://github.com/jdpreston30/thyroid-exposomics-2025.git
cd thyroid-exposomics-2025
```

```r
# renv activates automatically via .Rprofile
renv::restore()                 # first time only, ~10-20 min

# edit All_Run/config_dynamic.yaml so the paths match your system:
#   computers          — your machine's user_home and onedrive_path
#   paths.base_data_path — parent directory of the GC-MS raw data
# every other path is a template that populates from those two.

source("All_Run/run.R")
```

`renv::restore()` builds an isolated project library and leaves your system R untouched. It is
needed once per machine; the pipeline loads packages from `DESCRIPTION` on each run.

Chemical metadata and reference libraries live outside the repository, on OneDrive; their paths are
also set in `All_Run/config_dynamic.yaml`.

## Pipeline

`run.R` executes `R/Scripts/` in order, `00a` through `20`.

| Script | Purpose |
|---|---|
| `00a_environment_setup.R` | Conflict preferences and package loading from `DESCRIPTION` |
| `00b_setup.R` | Configuration; sources every utility under `R/Utilities/` |
| `00c_FTs.R` | Feature table import and preprocessing |
| `00d_peakwalk_compile.R` | PeakWalk compilation |
| `01_clinical_data.R` | Clinical data, demographics, Table 1 |
| `02_detection.R` | Detection frequency analysis |
| `03_classes.R` | Use classes and detection distribution |
| `04_variant_stats.R` | Statistical comparisons between tumor types |
| `05_variant_vis_prep.R` | Visualization data preparation |
| `06_tumor_cadaver.R` | Tumor vs. cadaver control comparison |
| `07_validation_prep.R` | Prepares manual spectral validation QC |
| `08_validation_run.R` | Manual spectral validation |
| `09_validation_plots_create.R` | Validation plot adjustment and manual review |
| `10_post_validation_clean.R` | Post-validation cleaning |
| `11_variant_vis.R` | Tumor-type differences, post-validation |
| `12_IARC_vis.R` | Tumor vs. control IARC carcinogen plots |
| `13_confounding_analysis.R` | Confounding analysis (age, sex, collection timing) |
| `14_render_figures.R` | Renders the main figures |
| `15_render_supplementary_figures.R` | Renders the supplementary figures |
| `16_tables.R` | Manuscript tables |
| `17_supplementary_tables.R` | Supplementary tables |
| `18_construct_supplementary.R` | Compiles the supplementary PDF |
| `19_results_validate.R` | Validates numerical claims against pipeline outputs |
| `20_session_info.R` | Writes `SESSION_INFO.txt` |

## Repository layout

```
├── DESCRIPTION                 # dependencies (CRAN, Bioconductor, GitHub)
├── renv.lock                   # pinned package versions
├── SESSION_INFO.txt            # session record from the manuscript run
├── All_Run/
│   ├── config_dynamic.yaml     # paths and analysis options — edit before running
│   └── run.R                   # pipeline entry point
├── R/
│   ├── Scripts/                # analysis workflow, 00a-20
│   └── Utilities/              # Analysis, Clinical, Helpers, Tabulation,
│                               #   Terminal, Validation, Visualization
├── Outputs/
│   ├── Figures/                # publication figures (PNG, TIFF, PDF)
│   ├── Tables/                 # manuscript tables
│   ├── Revisions/              # figures produced for the revision responses
│   └── Validation/             # spectral validation plots and PDFs
│       ├── failed/             #   compounds that failed validation
│       ├── initial_compile/    #   first-pass validation compilation
│       ├── revised/            #   revised validation plots (feed the supplement)
│       └── top_fragments/      #   top-fragment validations
└── Supplementary/
    ├── Components/             # Sections, Tables, Figures, References,
    │                           #   abbreviations.tsv
    ├── Build_Logs/             # LaTeX build logs
    └── Supplementary Data.pdf  # compiled supplement
```

## Contact

**Joshua D. Preston** — joshua.preston@emory.edu ·
[ORCID 0000-0001-9834-3017](https://orcid.org/0000-0001-9834-3017) ·
Emory University School of Medicine

**Corresponding author: M. Ryan Smith** — matthew.ryan.smith@emory.edu ·
[ORCID 0000-0002-8889-3477](https://orcid.org/0000-0002-8889-3477) ·
Emory University School of Medicine
