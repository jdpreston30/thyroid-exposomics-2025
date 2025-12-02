#+ Create subsetted table from master (Table 3)
    #- Bring in final info
      table_classes <- feature_metadata |>
        select(cas, Table_Subclass, Table_Qualifier)
    #- Create table columns and format
      table_3 <- MT_final |>
        #trouble populating OD-PABA so doing ident here manually
        mutate(
          annot_ident = if_else(cas == "21245-02-3", "Identification", annot_ident)
        ) |>
          mutate(
            Table_Class = case_when(
              Table_Class == "Insecticides and Pesticides" ~ "Insecticide/Pesticide",
              str_detect(Table_Class, "Dye intermediates") ~ "Dye intermediate",
              str_detect(Table_Class, "Chemical Synthesis Intermediates") ~ "Chemical Synthesis Intermediate",
              str_detect(Table_Class, "Carcinogenic Research Chemicals") ~ "Carcinogenic Research Chemical",
              str_detect(Table_Class, "Combustion Byproducts") ~ "Combustion Byproduct",
              str_detect(Table_Class, "Side-Reaction Byproducts") ~ "Side-Reaction Byproduct",
              str_detect(Table_Class, "Fungicides") ~ "Fungicide",
              str_detect(Table_Class, "Herbicides") ~ "Herbicide",
              str_detect(Table_Class, "Disinfectant Breakdown Products") ~ "Disinfectant Breakdown Product",
              str_detect(Table_Class, "Flavoring or Fragrance Agents") ~ "Flavoring or Fragrance Agent",
              str_detect(Table_Class, "Plasticizers and Plastic Additives") ~ "Plasticizer",
              str_detect(Table_Class, "Preservatives") ~ "Preservative",
              str_detect(Table_Class, "Flame Retardants") ~ "Flame Retardant",
              str_detect(Table_Class, "Plant Growth Regulators") ~ "Plant Growth Regulator",
              str_detect(Table_Class, "Humectants") ~ "Humectant",
              TRUE ~ Table_Class # Leave PFAS, Organic UV Filters, etc. unchanged
            )
          ) |>
        left_join(table_classes, by = "cas") |>
          mutate(
            Table_Class = ifelse(
              !is.na(Table_Subclass) & Table_Subclass != "",
              paste0(Table_Class, " (", Table_Subclass, ")"),
              Table_Class
            )
          ) |>
        mutate(
          short_name = paste0(
            short_name,
            if_else(annot_ident == "Annotation", "*", ""),
            if_else(
              Carcinogenicity %in% c("Likely Carcinogen", "Possible Carcinogen", "Known Carcinogen"),
              "\u2020", ""
            ),
            ifelse(Potential_EDC == "Y", "\u2021", "") # <-- FIXED HERE
          ),
          FTC_let = coalesce(FTC_let, FTC),
          FV_PTC_let = coalesce(FV_PTC_let, FV_PTC),
          PTC_let = coalesce(PTC_let, PTC),
          p_value = round(p_value, 3)
        ) |>
        mutate(short_name = str_replace(short_name, "NA$", "")) |>
        select(name_sub_lib_id,short_name, cas,annot_ident, Table_Class, FTC_let:PTC_let, p_value) |>
        # error with OD-PABA (identification) populating, so doing it manually
        mutate(
          short_name = if_else(cas == "21245-02-3", "OD-PABA", short_name),
          Table_Class = if_else(cas == "21245-02-3", "Organic UV Filter", Table_Class))
    #- Export
      write_xlsx(table_3, "table_3_export.xlsx")
  
    #! At this point, just exported and graphed these in Prism
    #- Run a ttest on log transformed data to determine signfiicant differences
      IARC_ttests <- full_joiner |>
        select(sample_ID, any_of(IARC_tumors_ctrl_filtered)) |>
        mutate(across(
          .cols = -c(sample_ID),
          .fns = ~ {
            col_min <- min(.x, na.rm = TRUE)
            replace(.x, is.na(.x), 0.5 * col_min)
          }
        )) |>
        mutate(across(where(is.numeric), ~ log2(. + 1e-6))) |>
        mutate(group = if_else(str_starts(sample_ID, "T00"), "Control", "Tumor")) |>
        pivot_longer(cols = where(is.numeric) & !any_of("group"),
                    names_to = "chemical", values_to = "value") |>
        group_by(chemical) |>
        summarise(
          p_value = t.test(value ~ group)$p.value,
          mean_control = mean(value[group == "Control"]),
          mean_tumor = mean(value[group == "Tumor"]),
          .groups = "drop"
        ) |>
        arrange(p_value)
        #! Data not shown but listed in manuscript results
  #+ Manual Spectral Validation (Supplementary)
    #- Create vector of variant significant CAS which are identifications
      sig_id <- table_3 |>
      filter(annot_ident == "Identification") |>
      mutate(id = str_extract(name_sub_lib_id, "(?<=_)[^_]+$")) |>
      pull(id)
    #- Import the spectral validation data
      spectral_validation_i <- read_excel(config$paths$primary_data, sheet = "lib.subject.summary") |>
        filter(id %in% sig_id) |>
        mutate(
          mzMed = paste0(
            round(mzMed, 4)
          ),
          rtMed = paste0(
            round(rtMed, 2)
          )
        ) |>
        select(
          name_sub_lib_id, lib, cas, name, short_display_name, id, subid, tmz,
          mzMed, trt, rtMed, libFragCnt, F1:F20
        ) |>
        arrange(cas)
    #- Now, systematically determine the ideal sample to validate on
      #_Select candidate sample columns
        f1_index <- which(names(spectral_validation_i) == "F1")
        f20_index <- which(names(spectral_validation_i) == "F20")
        candidate_cols <- names(spectral_validation_i)[f1_index:f20_index]
      #_Summarise n and min for each candidate column per ID
        fragment_summary <- spectral_validation_i |>
          group_by(id) |>
          summarise(
            across(all_of(candidate_cols), list(
              n = ~ sum(!is.na(.)),
              min = ~ if (sum(!is.na(.)) >= 2) min(., na.rm = TRUE) else NA_real_
            ), .names = "{.col}__{.fn}"),
            .groups = "drop"
          )
      #_Reshape long and extract best columns by fragment count + min value
        long_summary <- fragment_summary |>
          pivot_longer(-id, names_to = "measure", values_to = "value") |>
          separate(measure, into = c("column", "stat"), sep = "__") |>
          pivot_wider(names_from = stat, values_from = value)
        ideal_queries <- long_summary |>
          filter(n >= 2 & !is.na(min)) |>
          group_by(id) |>
          arrange(desc(n), desc(min)) |>
          summarise(
            ideal_query = first(column),
            ideal_query_backup = paste(head(column[-1], 5), collapse = ","),
            .groups = "drop"
          )
      #_Join back and finalize
        spectral_validation <- spectral_validation_i |>
          left_join(ideal_queries, by = "id") |>
          select(ideal_query, ideal_query_backup, id, everything()) |>
          arrange(id)
    #- Now, pare this down and join with the feature library
      #_Pare down
        spectral_validation_pared <- spectral_validation |>
          select(short_display_name,cas,ideal_query,ideal_query_backup,lib,id,subid,tmz,trt,mzMed,rtMed,libFragCnt) |>
          arrange(id) |>
          distinct()
      #_Get a sum of libfrags
      libfrag_summary <- spectral_validation_pared |>
        group_by(id) |>
        summarise(libFragCnt = sum(libFragCnt, na.rm = TRUE), .groups = "drop")
      #_Rejoin with feature library, export
        spectral_validation_wide <- spectral_validation_pared |>
          rename(
            mz = mzMed,
            rt = rtMed
          ) |>
          select(id, subid, mz, rt, ideal_query, ideal_query_backup) |>
          pivot_wider(
            names_from = subid,
            values_from = c(mz, rt),
            names_glue = "Exp_{.value}{subid}"
          ) |>
          relocate(id) |>
          left_join(feature_lib, by = "id") |>
          mutate(
            short_display_name = if_else(is.na(short_display_name), name, short_display_name)
          ) |>
          left_join(libfrag_summary, by = "id") |>
          select(
            short_display_name, cas, id, trt,
            starts_with("mz"),
            starts_with("Exp_mz"),
            starts_with("Exp_rt"),
            ideal_query, ideal_query_backup, libFragCnt
          )
      write_xlsx(spectral_validation_wide, "spectral_validation_wide.xlsx")
      #! At this point, exported to excel and formatted there