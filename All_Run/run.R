{
source("R/Utilities/Helpers/restore_renv.R")
source("R/Utilities/Helpers/load_dynamic_config.R")
config <- load_dynamic_config(computer = "auto", config_path = "All_Run/config_dynamic.yaml")
source("R/Scripts/00a_environment_setup.R")
source("R/Scripts/00b_setup.R")
source("R/Scripts/00c_clinical_data.R")
source("R/Scripts/00d_FTs.R")
source("R/Scripts/00e_peakwalk_compile.R")
source("R/Scripts/01_demographics.R")
source("R/Scripts/02_detection.R")
source("R/Scripts/03_classes.R")
source("R/Scripts/04_variant_stats.R")
source("R/Scripts/05_variant_vis_prep.R")
source("R/Scripts/06_tumor_cadaver.R")
source("R/Scripts/07_validation_prep.R")
source("R/Scripts/08_validation_run.R")
source("R/Scripts/09_validation_plots_create.R")
source("R/Scripts/10_post_validation_clean.R")
source("R/Scripts/11_variant_vis.R")
source("R/Scripts/12_IARC_vis.R")
source("R/Scripts/13_render_figures.R")
source("R/Scripts/14_render_supplementary_figures.R")
source("R/Scripts/15_tables.R")
source("R/Scripts/16_supplementary_tables.R")
source("R/Scripts/17_construct_supplementary.R")
}


source("R/Scripts/16_supplementary_tables.R")
# Read the ST1 LaTeX and insert it into tables.tex
st1_latex <- readLines("Outputs/Tables/ST1.tex")
tables_tex <- readLines("Supplementary/Components/Sections/tables.tex")
# Replace the placeholder line with actual ST1 content
tables_tex <- gsub("\\[INSERT ST1 HERE - TO BE GENERATED PROGRAMMATICALLY\\]", 
                   paste(st1_latex, collapse = "\n"), 
                   tables_tex)
writeLines(tables_tex, "Supplementary/Components/Sections/tables.tex")

# Now render the supplementary
source("R/Scripts/17_construct_supplementary.R")
