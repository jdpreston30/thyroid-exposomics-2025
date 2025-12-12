# GitHub Copilot Instructions

**R Code Style & Project Structure Preferences:**

## Comments
- Use hierarchical comment structure with symbols: `#*` (major section), `#+` (subsection), `#-` (detail), `#!` (important note)
- Number sections sequentially (e.g., `#* 8: Manual Spectral Validation`, `#+ 8.1: Convert files`, `#- 8.1.1: Setup`)
- Keep comments brief and descriptive
- Never use extra blank lines between sections unless absolutely necessary - code should be compact

## Project Architecture
- Use YAML configuration files for dynamic, computer-specific paths (e.g., `config_dynamic.yaml`)
- Store utility functions in organized subdirectories: `R/Utilities/Helpers/`, `R/Utilities/Validation/`, etc.
- Use renv for package management with automatic restore on first run
- Pipeline structure: Main run script (`All_Run/run.R`) sources numbered scripts sequentially
- Load all utility functions at setup with: `purrr::walk(list.files(...), source)`
- Store configuration in `.GlobalEnv$config` or similar for global access
- Use `load_dynamic_config()` pattern for automatic computer detection and path resolution

## Functions
- Use roxygen2-style documentation for utility functions
- Prefer tidyverse functions when appropriate
- Functions should respect YAML configuration flags (e.g., skip logic based on config parameters)
- Separate concerns: e.g., plot generation vs PDF compilation in separate functions

## File Organization
- Scripts: Numbered for pipeline order (00a, 00b, 01, 02, etc.)
- Outputs: Organized in `Outputs/Figures/`, `Outputs/Tables/`, `Outputs/Validation/`
- Metadata: Static reference files in `metadata_files/`
- Use OneDrive for final output storage, local temp dirs for run-specific I/O
