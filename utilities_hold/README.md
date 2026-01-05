# Utilities Hold - Quarantine Folder

This folder contains R utility functions that were identified as **completely unused** in the pipeline through static code analysis. These functions are quarantined here to allow pipeline testing without them before permanent deletion.

## Summary

- **Total functions moved**: 34 functions
- **Total files moved/extracted**: 17 files
- **Purpose**: Enable pipeline validation testing to confirm these functions are truly unnecessary
- **Next step**: Run full pipeline (`All_Run/run.R`) to verify it executes without errors

## Files in This Folder

### Entire Files Moved (No Used Functions)
These files contained only unused functions and were moved completely:

1. **add_file_references.R** - `add_file_references()` function
2. **classify_carcinogen.R** - `classify_carcinogenicity()` function  
3. **format_ppb_value.R** - `format_ppb_value()` function
4. **get_rt_range.R** - `get_rt_range()` function
5. **individual_fragments.R** - `individual_fragments()` function
6. **load_dynamic_config.R** - `load_dynamic_config()` function
7. **load_validation_plots.R** - `load_validation_plots()` function
8. **move_superscript.R** - `move_superscript()` function
9. **preview_table.R** - `preview_table()` function
10. **print_to_svg.R** - `print_to_svg()` function
11. **save_table_latex.R** - `save_table_latex()` function
12. **vp_top.R** - `vp_top()` function
13. **global_VP_adjust.R** - Duplicate copy (one copy remains in `R/Utilities/Validation/VP/`)

### Extracted Functions from Mixed-Content Files
These files contained both used and unused functions. Only unused functions were extracted:

1. **themes_unused_color_functions.R** (16 functions extracted from `R/Utilities/Visualization/themes.R`)
   - `get_EDC_color`, `get_IARC_color`, `get_carcinogen_color`, `get_tumor_noncancer_color`, `get_variant_color`
   - `scale_EDC_color`, `scale_EDC_fill`, `scale_IARC_color`, `scale_IARC_fill`
   - `scale_carcinogen_color`, `scale_carcinogen_fill`, `scale_tumor_noncancer_color`, `scale_tumor_noncancer_fill`
   - `scale_variant_color`, `scale_variant_fill`, `theme_thyroid`
   - **Remaining in original file**: Color palette definitions (e.g., `tep_palette`, `iarc_palette`)



3. **process_single_compound_unused.R** (extracted from `R/Utilities/Validation/rtx.R`)
   - Large worker function (739 lines) used internally by validation pipeline
   - **Remaining in original file**: Main `rtx()` function (confirmed used in pipeline)

## Functions Kept (Used in Pipeline)

These functions remain in their original locations as they are called by the pipeline:

**Helpers:**
- `adjust_plot()` - Modifies plot appearance (used in 8 scripts)
- `build_validation_table()` - Constructs validation tables
- `convert_superscript()` - Text formatting
- `fix_ST2_latex()` - LaTeX table fixes
- `fix_ST3_latex()` - LaTeX table fixes
- `fix_latex_header_fill()` - LaTeX formatting
- `restore_renv()` - Environment setup
- `superscript()` - Text formatting

**Visualization:**
- `plot_*()` functions - 10 plotting functions
- `print_to_png()` / `print_to_png_pdf()` - Output functions
- `figure_labels()` - Figure panel labels (vis_tools.R)
- Color palettes in themes.R (e.g., `tep_palette`)

**Validation & Analysis:**
- `rtx()` - Main validation processing function (rtx.R)
- Cross-referenced functions in VP/ subdirectory (used by `global_VP_adjust()`)

## Testing the Pipeline

To verify the pipeline works without the quarantined functions:

```r
# Run from R console in the project directory
source("All_Run/run.R")
```

**Expected behavior:**
- All scripts execute without errors related to missing functions
- Output files generated in `Outputs/` directories
- No warnings about undefined functions from utilities_hold

## If Pipeline Fails

If the pipeline fails due to a missing function:
1. Check which function is missing
2. Return the file from utilities_hold to its original location
3. Document the missing function (it was incorrectly classified as unused)

## Cleanup

Once pipeline testing confirms these functions are truly unused, permanently delete this folder:

```bash
rm -rf utilities_hold
```

---

**Audit Method**: Static code analysis using grep-based function extraction and pipeline call detection  
**Analysis Date**: 2025-01-05  
**Auditor**: Automated analysis tool  
**Cross-reference functions**: 9 functions identified as used internally by other utilities
