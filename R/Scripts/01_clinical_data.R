#* 1: Load and Clean Clinical Data
#+ 1.1: Import Clinical Data
#- 1.1.1: Tumor pathology data
tumor_pathology_raw <- read_excel(config$paths$tumor_pathology, sheet = "Logan Update")
#- 1.1.2: Cadaver metadata
cadaver_metadata_raw <- read_excel(config$paths$cadaver_metadata, sheet = "Validated_Pared")
#+ 1.2: Clean Clinical Data
#- 1.2.1: Tumor pathology data
clinical_data <- tumor_pathology_raw |>
  select(-T) |>
  assign_T_cat(ld_col = "LD", ete_col = "ETE", units = "cm", out_col = "T_stage_comp") |>
    mutate(year_bin = cut(as.numeric(year),
    breaks = seq(2006, 2022, length.out = 5),
    labels = c("2006-2009", "2010-2013", "2014-2017", "2018-2021"),
    include.lowest = TRUE
  )) |>
  mutate(
      Variant = factor(case_when(
        str_detect(Patient_ID, "^FVPTC\\d+$") ~ "IEFVPTC",
        str_detect(Patient_ID, "^F\\d+$") ~ "FTC",
        str_detect(Patient_ID, "^P\\d+$") ~ "PTC",
        TRUE ~ "Unknown"
      ), levels = c("FTC", "IEFVPTC", "PTC")),
      T_stage_comp = factor(T_stage_comp, levels = c("T1", "T2", "T3", "T4"), ordered = TRUE),
      N = factor(case_when(
        N == 0 ~ "N0", 
        N == 1 ~ "N1", 
        TRUE ~ NA_character_
      )),
      M = factor(case_when(
        M == 0 ~ "M0", 
        M == 1 ~ "M1", 
        TRUE ~ NA_character_
      )),
      Sex = factor(
        dplyr::case_when(
          is.na(Sex) ~ NA_character_,
          Sex == 1 ~ "Female",
          Sex == 0 ~ "Male"
        ),
        levels = c("Female", "Male")
      ),
      Age = as.numeric(Age)
    ) |>
  assign_AJCC8_stage(
      t_col = "T_stage_comp",
      n_col = "N", 
      m_col = "M",
      age_col = "Age",
      out_col = "Stage"
    ) |>
  select(Variant, Sex, Age, Sample_Collection_Timing = year_bin, T_Category = T_stage_comp, N_Category = N, M_Category = M, Pathologic_Stage = Stage, year)
#- 1.2.2: Cadaver metadata
cadaver_metadata <- cadaver_metadata_raw |>
  filter(seq != "006") |> # Did not run this sample
  mutate(collection_year = as.integer(format(DOD, "%Y"))) |> # extract collection year from DOD
  select(seq, age, sex, collection_year)
  
