
  
    #! At this point, just exported and graphed these in Prism

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