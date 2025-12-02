#* Analysis
  #+ IARC and EDC Classes (Figure 2B and C)
    IARC <- feature_metadata %>%
      select(IARC_Group) %>%
      count(IARC_Group) %>%
      arrange(desc(n))
    EDC <- feature_metadata %>%
      select(Potential_EDC) %>%
      count(Potential_EDC) %>%
      arrange(desc(n))
    #! Graphed in prism at this point
  #+ Chemical Superclasses and Classes (Supl Figure 1A and B)
    Superclass <- feature_metadata %>%
      select(Superclass) %>%
      count(Superclass) %>%
      arrange(desc(n))
    Class <- feature_metadata %>%
      select(Class) %>%
      count(Class) %>%
      arrange(desc(n))
    write.csv(Superclass,"superclass.csv")
    write.csv(Class, "Class.csv")
  #+ ANOVA and Fisher's (3A and 3B)
    #- ANOVA Stats
      #_Import the weights of the tissues
        weights <- read_excel("Data and Metadata Files/primary_data.xlsx", sheet = "tissue_weights") %>%
          filter(samples == "Tumor") %>%
          select(ID, weight_mg) %>%
          rename(patient_ID = ID)
      #_Normalize by weights
        tumors_quant_wt_i <- tumor_column %>%
          pivot_longer(-name_sub_lib_id, names_to = "ID", values_to = "Value") %>%
          pivot_wider(names_from = name_sub_lib_id, values_from = Value) %>%
          mutate(patient_ID = ID) %>%
          left_join(tumor_seq, by = "ID") %>%
          select(patient_ID,all_of(cols_quant)) %>%
          mutate(across(-c(variant,patient_ID), ~ as.numeric(.))) %>%
          mutate(across(-c(variant,patient_ID), ~ ifelse(is.na(.), 0.5 * min(., na.rm = TRUE), .))) %>%
          left_join(weights, by = "patient_ID") %>%  # join tissue weights in
          mutate(across(where(is.numeric) & !matches("weight_mg"), ~ . / weight_mg)) %>%
          select(-c(weight_mg,patient_ID)) %>%
          mutate(across(-variant,~log2(.)))
      #_Filter out endogenous features
        cas_key_endog <- tumor_raw %>%
          select(name_sub_lib_id, cas) %>%
          filter(cas %in% endog_cas) %>%
          pull(name_sub_lib_id)
        tumors_quant_wt <- tumors_quant_wt_i %>%
          select(-any_of(cas_key_endog))
      #_Run ANOVA, subset FT, save
        anova_results_sig <- tumors_quant_wt %>%
          pivot_longer(-variant, names_to = "name_sub_lib_id", values_to = "value") %>%
          group_by(name_sub_lib_id) %>%
          summarise(
            p_value = broom::tidy(aov(value ~ variant))$p.value[1] # Extract p-value from ANOVA
          ) %>%
          filter(p_value < 0.05) %>%
          arrange(p_value) %>%
          arrange(name_sub_lib_id) %>%
          pull(name_sub_lib_id) # Extract significant feature names
      #_Filter tumors_quant to keep only significant compounds and apply z-score
        tumors_quant_sig <- tumors_quant_wt %>%
          select(variant, all_of(anova_results_sig)) %>% 
          mutate(across(-variant, ~ scale(.)[,1]))
      #_Make a second cas key
        cas_key_2 <- tumor_raw %>%
          select(name_sub_lib_id, cas,short_display_name) %>%
          rename(Name = short_display_name)
      #_Create summary table with reran ANOVA
        summary_table_i <- tumors_quant_sig %>%
          pivot_longer(-variant, names_to = "name_sub_lib_id", values_to = "value") %>%
          group_by(name_sub_lib_id, variant) %>%
          mutate(variant = if_else(variant == "FV-PTC", "FV_PTC", variant)) %>%
          summarise(
            mean_sd = paste0(round(mean(value, na.rm = TRUE), 2), " ± ", round(sd(value, na.rm = TRUE), 2)),
            .groups = "drop"
          ) %>%
          pivot_wider(names_from = variant, values_from = mean_sd) %>%
          left_join(
            tumors_quant_sig %>%
              pivot_longer(-variant, names_to = "name_sub_lib_id", values_to = "value") %>%
              group_by(name_sub_lib_id) %>%
              summarise(
                p_value = round(summary(aov(value ~ variant))[[1]][["Pr(>F)"]][1], 4),
                .groups = "drop"
              ),
            by = "name_sub_lib_id"
          ) %>%
          select(name_sub_lib_id, Follicular, FV_PTC, Papillary, p_value) %>%
          arrange(p_value) %>%
          mutate(mode = "quantitative") %>%
          left_join(cas_key_2, by = "name_sub_lib_id") %>%
          select(Name, cas, name_sub_lib_id, mode, p_value, everything()) %>%
          arrange(p_value)
    #- Fishers
      #_ Run Fisher's
        fisher_results <- tumors_qual %>%
          pivot_longer(-variant, names_to = "name_sub_lib_id", values_to = "detected") %>%
          group_by(name_sub_lib_id) %>%
          summarise(
            p_value = tryCatch(
              fisher.test(table(variant, detected))$p.value,
              error = function(e) NA_real_
            ),
            Follicular = (sum(detected[variant == "Follicular"]) / sum(variant == "Follicular")) * 100,
            Papillary = (sum(detected[variant == "Papillary"]) / sum(variant == "Papillary")) * 100,
            FV_PTC = (sum(detected[variant == "FV-PTC"]) / sum(variant == "FV-PTC")) * 100
          ) %>%
          arrange(p_value) %>%
          filter(p_value < 0.05) %>%
          arrange(name_sub_lib_id) %>%
          arrange(p_value) %>%
          left_join(cas_key_2, by = "name_sub_lib_id") %>%
          mutate(mode = "qualitative") %>%
          select(Name,everything())
    #- Visually inspect for duplicates (with different fragments) between the ANOVA and Fisher's
      #_bind names and modes to check for duplicates
        summary_table_i_dup <- summary_table_i %>%
          select(Name, cas, name_sub_lib_id, mode)
        fisher_results_dup <- fisher_results %>%
          select(Name, cas, name_sub_lib_id, mode)
      #_Add duplicate check column
        dupl_check <- rbind(summary_table_i_dup, fisher_results_dup) %>%
          arrange(cas) %>%
            mutate(
              duplicate = case_when(
                Name %in% Name[base::duplicated(Name)] |
                  cas %in% cas[base::duplicated(cas)] ~ "*",
                TRUE ~ ""
              )
            ) %>%
          arrange(Name) %>%
          arrange(desc(duplicate))
      #_Inspect          
        print(dupl_check, n = Inf)
      #_Duplicates which have both quant and qual frags (quant chosen)
        dupl_check_qual_quant <- dupl_check %>%
          group_by(Name, cas) %>%
          filter(n() == 2, sum(mode == "qualitative") == 1, sum(mode == "quantitative") == 1) %>%
          ungroup()
      #! Benz(a)anthracene_1_BP2.GC2_CP2220 (quant) chosen over Benz(a)anthracene_0_BP3.GC2_CP3027
      #! 2-amino-4-nitrotoluene_1_BP3.GC2_CP3003 (quant) chosen over 2-amino-4-nitrotoluene_0_BP3.GC2_CP3003
      #_Duplicates which have both qual
        # Pull these features
        dupl_check_qual <- dupl_check %>%
          group_by(Name, cas) %>%
          filter(n() > 1, all(mode == "qualitative")) %>%
          ungroup()
        # Subset to only those, then systematically choose lowest detection fragment to remove
        qual_only_ids <- dupl_check_qual$name_sub_lib_id
        qual_dup_removed_frags <- fisher_results %>%
          filter(name_sub_lib_id %in% qual_only_ids) %>%
          mutate(sum_det = Follicular + Papillary + FV_PTC) %>%
          select(Name, cas, name_sub_lib_id, sum_det) %>%
          group_by(Name) %>%
          slice_min(sum_det, with_ties = FALSE) %>%
          ungroup()
        print(qual_dup_removed_frags, n = Inf)
      #! Aldicarb: aldicarb_0_BP2.GC2_CP2507 > detection than aldicarb_1_BP2.GC2_CP2507
      #! Fenvalerate: Fenvalerate_1_3_BP3.GC2_CP3164 > detection than Fenvalerate_2_3_BP3.GC2_CP3165
      #! Octyl-dimethyl-PABA: Octyl-dimethyl-PABA_1_BP2.GC2_CP2331 > detection than Octyl-dimethyl-PABA_2_BP2.GC2_CP2331
      #! tris(tribromoneopentyl): tris(tribromoneopentyl)_2_BP2.GC2_CP2302 > detection than tris(tribromoneopentyl)_3_BP2.GC2_CP2302
      #_Pull the fragment names to remove from both
        # qual/quant duplicates
        qual_quant_dupl_removed <- dupl_check_qual_quant %>%
          filter(mode == "qualitative") %>%
          pull(name_sub_lib_id)
        # qual/qual duplicates
        qual_qual_dupl_removed <- qual_dup_removed_frags %>%
          pull(name_sub_lib_id)
    #- Remove duplicate fragments with lower detection or qual when a quant is available:
      qual_single_frag <- fisher_results %>%
        filter(!name_sub_lib_id %in% qual_qual_dupl_removed) %>%
        filter(!name_sub_lib_id %in% qual_quant_dupl_removed) %>%
        select(cas,name_sub_lib_id,mode,everything()) %>%
        arrange(Name)
        print(qual_single_frag, n = Inf)
    #- Compiled ANOVA and Fisher's, add short name, export qual for Prism
      #_Add updated quant names to the excel
        short_name <- read_excel("Data and Metadata Files/chemical_metadata.xlsx", sheet = "feature_metadata") %>%
          select(cas, name, Potential_EDC, IARC_Group, GHS_var_diff_only,Short_display_name, Graph_Class, Superclass, Table_Class, Annotation_or_Identification) %>%
          rename(annot_ident = Annotation_or_Identification)
      #_Filter quant to be safe
        summary_table <- summary_table_i %>%
          filter(!name_sub_lib_id %in% qual_qual_dupl_removed) %>%
          filter(!name_sub_lib_id %in% qual_quant_dupl_removed)
      #_Bind and name
        quant_qual_results <- rbind(qual_single_frag, summary_table) %>%
          left_join(short_name, by = "cas") %>%
          rename(usage_class = Graph_Class) %>%
          select(cas, name, Short_display_name, name_sub_lib_id, mode, p_value, Follicular, FV_PTC, Papillary, everything()) %>%
          rename(FTC = Follicular, PTC = Papillary) %>%
          mutate(across(c(FTC, FV_PTC, PTC),
            ~ if_else(mode == "qualitative", paste0(.x, "%"), as.character(.x)),
            .names = "{.col}"
          )) %>%
          rename(short_name = Short_display_name)
      #_Filter to only the qualitative
        qual_prism_i <- quant_qual_results %>%
          filter(mode == "qualitative") %>%
          arrange(p_value) %>%
          mutate(across(c(FTC, FV_PTC, PTC), ~ as.numeric(gsub("%", "", .)) / 100)) %>% 
          mutate(short_name = if_else(
          annot_ident == "Annotation",
          paste0(short_name, "\u2020"),
          short_name)) %>%
          select(short_name,FTC:PTC)
      #_Make it a matrix
        qual_matrix <- qual_prism_i %>%
          select(-short_name) %>%
          as.matrix()
      #_Perform hierarchical clustering on rows
        row_clusters <- hclust(dist(qual_matrix))
      #_Reorder the rows based on the clustering
        qual_prism_i_reordered <- qual_prism_i %>%
          slice(row_clusters$order) %>% #multiply by 100
          mutate(across(c(FTC,FV_PTC,PTC), ~ round(.x * 100, 2))) %>%
          mutate(order = row_number()) %>%
          arrange(desc(order)) %>%
          select(-order)
      #_Export as an excel to preserve formatting
        write_xlsx(qual_prism_i_reordered, "qual_prism_clusters_ordered.xlsx")
    #- Get the top 5 quant for Prism
      #_Get the top 5 compounds
        top_5_quant <- summary_table %>%
          arrange(p_value) %>%
          slice_head(n = 5) %>%
          select(name_sub_lib_id, p_value) %>%
          pull(name_sub_lib_id)
      #_Extract the CAS numbers for the top 5 compounds, merge
        column_renames <- tibble(name_sub_lib_id = top_5_quant) %>%
          left_join(cas_key_2, by = c("name_sub_lib_id" = "name_sub_lib_id")) %>%
          left_join(short_name, by = "cas") %>%
          mutate(
            Short_display_name = ifelse(annot_ident == "Annotation", 
                                        paste0(Short_display_name, "\u2020"), 
                                        Short_display_name)
          ) %>%
          select(name_sub_lib_id, Short_display_name) %>%
          deframe()
      #_Now subset the z-scored quant data to the top 5,rename
        tumors_quant_sig_top5 <- tumors_quant_sig %>%
          select(variant, all_of(top_5_quant)) %>%
          mutate(variant = factor(variant, levels = c("Follicular", "FV-PTC", "Papillary"))) %>%
          arrange(variant) %>%
          rename_with(~ column_renames[.x], .cols = all_of(names(column_renames)))
      #_Export as excel to preserve formatting
        write_xlsx(tumors_quant_sig_top5, "top_5_quant.xlsx")
  #+ Create pre-master table and p-values for 3A
    #- Premaster posthoc
      #_Perform post-hoc test on quants and notate with compact letter display on values
        posthoc_quant <- tumors_quant_sig %>%
          mutate(variant = as.factor(variant)) %>%
          mutate(variant = recode(variant,
            "Papillary" = "PTC",
            "Follicular" = "FTC",
            "FV-PTC" = "FV_PTC"
          ))
      #_Helper to make superscripts for 'a', 'b', 'c'...
        superscript <- function(letter) {
          sapply(strsplit(letter, NULL)[[1]], function(char) {
            switch(char,
              "a" = "\u1d43",
              "b" = "\u1d47",
              "c" = "\u1d9c",
              "d" = "\u1d48",
              "e" = "\u1d49",
              "f" = "\u1da0",
              char # fallback if undefined
            )
          }) %>% paste0(collapse = "")
        }
      # _Get the names for the quant data
        posthoc_compound_names <- posthoc_quant %>%
          select(where(is.numeric)) %>%
          names()
      #_Main posthoc creation
        posthoc_table <- dplyr::bind_rows(
          lapply(posthoc_compound_names, function(compound) {
            formula <- as.formula(paste0("`", compound, "` ~ variant"))
            aov_model <- aov(formula, data = posthoc_quant)
            p_value <- summary(aov_model)[[1]][["Pr(>F)"]][1]
            # Get group means
            means <- posthoc_quant %>%
              group_by(variant) %>%
              summarise(mean = mean(.data[[compound]], na.rm = TRUE), .groups = "drop") %>%
              deframe()
            # Generate Tukey group letters
            letters <- multcompLetters(TukeyHSD(aov_model)[[1]][, "p adj"])$Letters
            # Create formatted output with superscript letters
            tibble(
              name_sub_lib_id = compound,
              FTC = paste0(sprintf("%.2f", means["FTC"]), superscript(letters["FTC"])),
              FV_PTC = paste0(sprintf("%.2f", means["FV_PTC"]), superscript(letters["FV_PTC"])),
              PTC = paste0(sprintf("%.2f", means["PTC"]), superscript(letters["PTC"])),
              p_value = round(p_value, 5)
            )
          }))%>% 
          left_join(cas_key_2, by = "name_sub_lib_id") %>%
          select(cas, FTC, FV_PTC, PTC) %>%
          rename(FTC_let = FTC, FV_PTC_let = FV_PTC, PTC_let = PTC)
      #_Write function to convert superscript letters to regular letters
        convert_superscript <- function(x) {
          str_replace_all(x, c(
            "ᵃ" = "a",
            "ᵇ" = "b",
            "ᶜ" = "c",
            "ᵈ" = "d",
            "ᵉ" = "e",
            "ᶠ" = "f"
          ))
        }
      #_Separate superscripts from means
        posthoc_decomp <- posthoc_table %>%
          # Extract means and letters
          mutate(
            FTC = as.numeric(str_extract(FTC_let, "-?[0-9.]+")),
            FV_PTC = as.numeric(str_extract(FV_PTC_let, "-?[0-9.]+")),
            PTC = as.numeric(str_extract(PTC_let, "-?[0-9.]+")),
            FTC_L = convert_superscript(str_extract(FTC_let, "[ᵃᵇᶜᵈᵉᶠ]+")),
            FV_PTC_L = convert_superscript(str_extract(FV_PTC_let, "[ᵃᵇᶜᵈᵉᶠ]+")),
            PTC_L = convert_superscript(str_extract(PTC_let, "[ᵃᵇᶜᵈᵉᶠ]+"))
          ) %>%
          select(cas, FTC:PTC_L) %>%
          rowwise()
      #_Systematically determine the ranking
        posthoc_ranked <- posthoc_decomp %>%
        rowwise() %>%
        mutate(
          ranking = {
            # Create a data frame for the current row
            groups <- tibble(
              group = c("FTC", "FV_PTC", "PTC"),
              mean = c(FTC, FV_PTC, PTC),
              letter = c(FTC_L, FV_PTC_L, PTC_L)
            )

            # Order by means (descending) and letters (ascending)
            ordered_groups <- groups %>%
              arrange(desc(mean), letter)

            # Construct the ranking string with correct symbols
            ranking_string <- paste(
              ordered_groups$group[1],
              if (ordered_groups$letter[1] == ordered_groups$letter[2]) "=" else ">",
              ordered_groups$group[2],
              if (ordered_groups$letter[2] == ordered_groups$letter[3]) "=" else ">",
              ordered_groups$group[3]
            )

            ranking_string
          }
        ) %>%
        ungroup()
      #_Extract the highest group
        highest_posthoc <- posthoc_ranked %>%
          rowwise() %>%
          mutate(
            highest = {
              # Extract the ranking from the previously computed 'ranking' column
              ranking_parts <- strsplit(ranking, " > ")[[1]]

              # Check if all groups are equal
              if (length(ranking_parts) == 1) {
                "Equal"
              } else {
                # Extract the first part before the first ">"
                first_part <- ranking_parts[1]

                # If the first part has "=", return both groups
                if (grepl("=", first_part)) {
                  strsplit(first_part, " = ")[[1]] %>% paste(collapse = ", ")
                } else {
                  first_part
                }
              }
            }
          ) %>%
          select(cas, highest,ranking) %>%
          ungroup()
      #_Join posthoc with original, and add simple ranking for qual
       MTii <- quant_qual_results %>%
         left_join(posthoc_table, by = "cas") %>%
         left_join(highest_posthoc %>% select(cas, highest), by = "cas") %>%
         select(short_name, FTC:PTC, highest, everything()) %>%
         # Update NA 'highest' for qualitative mode
         rowwise() %>%
         mutate(
           highest = ifelse(
             is.na(highest) & mode == "qualitative",
             {
               # Extract percentages as numeric
               FTC_num <- as.numeric(str_replace(FTC, "%", ""))
               FV_PTC_num <- as.numeric(str_replace(FV_PTC, "%", ""))
               PTC_num <- as.numeric(str_replace(PTC, "%", ""))

               # Create a ranking tibble
               groups <- tibble(
                 group = c("FTC", "FV_PTC", "PTC"),
                 percent = c(FTC_num, FV_PTC_num, PTC_num)
               ) %>%
                 arrange(desc(percent))

               # Handle ties for the highest percentage
               if (groups$percent[1] == groups$percent[2]) {
                 paste(groups$group[groups$percent == max(groups$percent)], collapse = ", ")
               } else {
                 groups$group[1]
               }
             },
             highest
           )
         ) %>%
         ungroup()
    #- Compute actual p-values for top 5 post-hocs which are displayed on 3A as stars
      #_Get values for top 5 ANOVA
        posthoc_table_pvalues <- dplyr::bind_rows(
          lapply(posthoc_compound_names, function(name_sub_lib_id) {
            formula <- as.formula(paste0("`", name_sub_lib_id, "` ~ variant"))
            aov_model <- aov(formula, data = posthoc_quant)
            p_value <- summary(aov_model)[[1]][["Pr(>F)"]][1]
            tibble(name_sub_lib_id = name_sub_lib_id, p_value = p_value)
          })
          ) %>%
          arrange(p_value) %>%
          rowwise() %>%
          mutate(
            pairwise = list({
              name_sub_lib_id <- name_sub_lib_id
              formula <- as.formula(paste0("`", name_sub_lib_id, "` ~ variant"))
              aov_model <- aov(formula, data = posthoc_quant)
              tukey <- TukeyHSD(aov_model)[["variant"]]
              pvals <- as.data.frame(tukey)[, "p adj"]
              names(pvals) <- rownames(tukey)
              signif_labels <- sapply(pvals, function(p) {
                if (p < 0.001) {
                  "***"
                } else if (p < 0.01) {
                  "**"
                } else if (p < 0.05) {
                  "*"
                } else {
                  "ns"
                }
              })
              tibble(
                comparison = names(signif_labels),
                star = signif_labels
              )
            })
          ) %>%
          unnest(pairwise) %>%
          select(name_sub_lib_id, p_value,comparison, star) %>%
          left_join(cas_key_2, by = "name_sub_lib_id") %>%
          select(-name_sub_lib_id) %>%
          left_join(MTii %>% select(cas,FTC_let:PTC_let), by = "cas") %>%
          slice_head(n = 15)
      #! At this point, just screenshotted and did mannually in BioRender
  #+ Balloon plot usage classes (3C)
    #- Prepare Data for Balloon Plot
      #_Organize data
        balloon_data <- MTii %>%
          filter(!is.na(highest)) %>%
          select(cas, usage_class, highest) %>%
          # Expand rows for double counting when multiple highest groups
          mutate(
            highest = strsplit(highest, ", ")
          ) %>%
          unnest(highest) %>%
          # Map the group names to the variant labels
          mutate(
            Variant = recode(highest,
              "FTC" = "Follicular",
              "FV_PTC" = "FV-PTC",
              "PTC" = "Papillary"
            )
          ) %>%
          # Remove 'Equal' and NA values
          filter(Variant != "Equal", !is.na(Variant)) %>%
          # Count occurrences for each variant within each usage class
          count(usage_class, Variant) %>%
          mutate(
            usage_class = str_replace_all(usage_class, "[†‡]", "")
          )
      #_Determine the correct order for usage_class
        usage_class_order <- balloon_data %>%
          pivot_wider(names_from = Variant, values_from = n, values_fill = 0) %>%
          arrange(Papillary, `FV-PTC`, Follicular) %>%
          pull(usage_class)
      #_Apply the correct order for y-axis
        balloon_data_graph <- balloon_data %>%
          mutate(
            Variant = factor(Variant, levels = c("Follicular", "FV-PTC", "Papillary")),
            usage_class = factor(usage_class, levels = usage_class_order)
          )
    #- Create and Save the Balloon Plot
      #_Set Variant Colors
        variant_colors <- c(
          "Follicular" = "#294b88",
          "FV-PTC" = "#23744d",
          "Papillary" = "#df8d09"
        )
      #_Plot
        balloon <- ggplot(balloon_data_graph, aes(x = Variant, y = usage_class, size = n, fill = Variant)) +
          geom_point(shape = 21, color = "black") + # Black outline for visibility
          scale_size_continuous(
            range = c(0.9, 4.5),
            limits = c(0, max(balloon_data$n, na.rm = TRUE)),
            breaks = c(1, 3, 5, 7),
            labels = as.character(c(1, 3, 5, 7))
          ) +
          scale_fill_manual(values = variant_colors, guide = "none") + # Apply custom colors & remove Variant legend
          labs(size = "Count") +
          theme_minimal() +
          theme(
            text = element_text(family = "Arial", face = "bold", size = 14, color = "black"),
            axis.text.x = element_blank(),
            # axis.text.x = element_text(angle = 45, hjust = 1, family = "Arial", face = "bold", size = 10.25, color = "black"),
            axis.text.y = element_text(family = "Arial", face = "bold", size = 8, color = "black"),
            legend.text = element_text(family = "Arial", face = "bold", size = 10, color = "black"),
            legend.title = element_text(size = 10),
            axis.title.x = element_blank(),
            axis.title.y = element_blank(),
            panel.grid.major = element_line(color = "grey80", linewidth = 0.3, linetype = "solid"),
            panel.grid.minor = element_line(color = "grey90", linewidth = 0.2, linetype = "solid")
          ) +
          guides(size = guide_legend(override.aes = list(fill = "white", color = "black", stroke = 0.5)))
      #_Save  
        ggsave("balloon3C.svg", plot = balloon, width = 4.4, height = 3, dpi = 1200)
  #+ Carcinogen classification based on GHS and IARC (3E)
    #- Define function to classify carcinogenicity
      classify_carcinogenicity <- function(ghs, iarc) {
        # Treat explicit "NA" as true NA for consistency
        if (is.na(ghs) || ghs == "NA" || ghs == "no ghs carcinogen statement") {
          ghs <- NA
        }
        if (is.na(iarc) || iarc == "NA" || iarc == "Not classified") {
          iarc <- NA
        }
        # Extract numeric values from H350, H350i, H351
        h350_perc <- as.numeric(str_extract(tolower(ghs), "(?<=h350 \\()\\d+\\.?\\d*"))
        h350i_perc <- as.numeric(str_extract(tolower(ghs), "(?<=h350i \\()\\d+\\.?\\d*"))
        h351_perc <- as.numeric(str_extract(tolower(ghs), "(?<=h351 \\()\\d+\\.?\\d*"))
        # Known Carcinogen (Highest Priority)
        if (!is.na(iarc) && iarc == "1") {
          return("Known Carcinogen")
        }
        # Likely Carcinogen (High Priority)
        if (!is.na(iarc) && iarc == "2A") {
          return("Likely Carcinogen")
        }
        if (!is.na(h350_perc) && h350_perc >= 50) {
          return("Likely Carcinogen")
        }
        if (!is.na(h350i_perc) && h350i_perc >= 50) {
          return("Likely Carcinogen")
        }
        # Possible Carcinogen (Moderate Priority)
        if (!is.na(iarc) && iarc == "2B") {
          return("Possible Carcinogen")
        }
        if (!is.na(h350_perc) && h350_perc > 0 && h350_perc < 50) {
          return("Possible Carcinogen")
        }
        if (!is.na(h350i_perc) && h350i_perc > 0 && h350i_perc < 50) {
          return("Possible Carcinogen")
        }
        if (!is.na(h351_perc) && h351_perc > 0) {
          return("Possible Carcinogen")
        }
        # Uncertain Risk (Low Priority)
        if (!is.na(iarc) && iarc == "3" && is.na(h350_perc) && is.na(h350i_perc) && is.na(h351_perc)) {
          return("Uncertain Risk")
        }
        # Unclassified (No IARC and No GHS)
        if (is.na(iarc) && is.na(ghs)) {
          return("Unclassified")
        }
        # Catch-all for unexpected cases
        return("Uncertain Risk")
      }
    #- Apply function to classify carcinogenicity
      MTi <- MTii %>%
        mutate(
          Carcinogenicity = mapply(
            classify_carcinogenicity,
            ghs = GHS_var_diff_only,
            iarc = IARC_Group
          )
        ) %>%
        select(short_name,highest, Carcinogenicity,IARC_Group,GHS_var_diff_only,everything()) %>%
        arrange(Carcinogenicity)
    #- Leftjoin and summarize for Prism
      carc_by_variant <- MTi %>%
        # Expand rows for double counting when multiple highest groups
        mutate(
          highest = strsplit(highest, ", ")
        ) %>%
        unnest(highest) %>%
        # Map the group names to the variant labels
        mutate(
          Variant = recode(highest,
            "FTC" = "Follicular",
            "FV_PTC" = "FV-PTC",
            "PTC" = "Papillary"
          )
        ) %>%
        # Count occurrences for each Variant and Carcinogenicity
        count(Variant, Carcinogenicity) %>%
        # Pivot to wide format
        pivot_wider(
          names_from = Carcinogenicity,
          values_from = n,
          values_fill = 0 # Fill missing counts with 0
        ) %>%
        # Reorder columns for readability
        select(
          Variant,
          "Known Carcinogen",
          "Likely Carcinogen",
          "Possible Carcinogen",
          "Uncertain Risk",
          "Unclassified",
          everything()
        ) %>%
        filter(Variant != "Equal")
      write.csv(carc_by_variant, "carc_by_variant.csv")
  #+ Chemical Library Display (ST1)
    #- Create feature library table which has subid pivoted columns and each id is treated individually
      feature_lib <- read_excel("Data and Metadata Files/primary_data.xlsx", sheet = "library") %>%
          filter(Disposition != "Endogenous") %>%
          mutate(subid_col = paste0("mz", subid)) %>%
          select(id, name, short_display_name, trt, monoisotopic, cas, formula, Disposition, subid_col, tmz) %>%
          distinct() %>%
          pivot_wider(
            names_from = subid_col,
            values_from = tmz) %>%
          arrange(cas)
      write_xlsx(feature_lib, "feature_library.xlsx", col_names = TRUE)
      #! In excel, then pared down and formatted, but reimporting here to double check the molecular formulas
  #+ Determine est PPM/PPB for detected samples and controls, join (ST3/4)
    #- Import the quantified data
      #_Pull the differing features by variant
        diff_by_var <- MTi %>%
          pull(name_sub_lib_id)
        diff_by_var_quant <- MTi %>%
          filter(mode == "quantitative") %>%
          pull(name_sub_lib_id)
        diff_by_var_qual <- MTi %>%
          filter(mode == "qualitative") %>%
          pull(name_sub_lib_id)
      #_Import raw concentration data
        conc_raw <- read_excel("Data and Metadata Files/primary_data.xlsx", sheet = "lib.subject.qsummary",col_type = "text") %>%
          select(name_sub_lib_id,F1:F20) %>%
          mutate(across(-name_sub_lib_id, as.numeric))
    #- Create a mean PPM/PPB in detected samples value for each chemical, add to master table
      #_Pivot long and normalize by weight
        ppm_ppb_long <- conc_raw %>%
          pivot_longer(cols = -name_sub_lib_id, names_to = "patient_ID", values_to = "Ce") %>%
          mutate(C = as.numeric(gsub(",", "", Ce))) %>%
          left_join(weights, by = "patient_ID") %>%
          mutate(
            PPM = (Ce * 10^2) / weight_mg,
            PPB = (Ce * 10^5) / weight_mg
          ) %>%
          select(name_sub_lib_id, patient_ID, PPM, PPB)
      #_Convert to PPM and PPB feature tables
        ppm_raw <- ppm_ppb_long %>%
          select(name_sub_lib_id, patient_ID, PPM) %>%
          pivot_wider(names_from = patient_ID, values_from = PPM)
        ppb_raw <- ppm_ppb_long %>%
          select(name_sub_lib_id, patient_ID, PPB) %>%
          pivot_wider(names_from = patient_ID, values_from = PPB)
      #_Now, use these to make estimated PPM and PPB per chemical in DETECTED samples
        ppm_ppb_summary_detected <- ppm_raw %>%
          rowwise() %>%
          mutate(
            pct_det = sum(!is.na(c_across(-c(name_sub_lib_id)))) / (ncol(.) - 1) * 100,
            PPM_mean = mean(c_across(-name_sub_lib_id), na.rm = TRUE),
            PPM_sd = sd(c_across(-name_sub_lib_id), na.rm = TRUE),
            PPM = paste0(format(PPM_mean, digits = 3), " ± ", format(PPM_sd, digits = 3))
          ) %>%
          select(name_sub_lib_id, PPM,pct_det) %>%
          left_join(
            ppb_raw %>%
              rowwise() %>%
              mutate(
                PPB_mean = mean(c_across(-name_sub_lib_id), na.rm = TRUE),
                PPB_sd = sd(c_across(-name_sub_lib_id), na.rm = TRUE),
                PPB = paste0(format(PPB_mean, digits = 3), " ± ", format(PPB_sd, digits = 3))
              ) %>%
              select(name_sub_lib_id, PPB),
            by = "name_sub_lib_id"
          ) %>%
          left_join(cas_key_2, by = "name_sub_lib_id") %>%
          select(cas, everything(),-c(Name))
    #- Determine PPM and PPB for controls
      #_Import control tissue weights
        cadaver_tissue_wts <- read_excel("Data and Metadata Files/primary_data.xlsx", sheet = "tissue_weights") %>%
          filter(samples == "Control") %>%
          select(ID, weight_mg) %>%
          rename(control_ID = ID)
      #_Import control conc data, normalize by tissue weight, compute PPM/PPB
        cadaver_qraw <- read_excel("Data and Metadata Files/primary_data.xlsx", sheet = "lib.subject.qsummary.cadaver") %>%
          select(name_sub_lib_id, T001:T009) %>%
          mutate(across(-name_sub_lib_id, as.numeric)) %>%
          pivot_longer(cols = -name_sub_lib_id, names_to = "control_ID", values_to = "Ce") %>%
          mutate(C = as.numeric(gsub(",", "", Ce))) %>%
          left_join(cadaver_tissue_wts, by = "control_ID") %>%
          mutate(
            PPM = (Ce * 10^2) / weight_mg,
            PPB = (Ce * 10^5) / weight_mg
          ) %>%
          select(name_sub_lib_id, control_ID, PPM, PPB)
      #_Convert to PPM and PPB feature tables
        ppm_raw_ctrl <- cadaver_qraw %>%
          select(name_sub_lib_id, control_ID, PPM) %>%
          pivot_wider(names_from = control_ID, values_from = PPM)
        ppb_raw_ctrl <- cadaver_qraw %>%
          select(name_sub_lib_id, control_ID, PPB) %>%
          pivot_wider(names_from = control_ID, values_from = PPB)
      #_Now, use these to make estimated PPM and PPB per chemical in DETECTED samples
        ppm_ppb_summary_detected_ctrl <- ppm_raw_ctrl %>%
          rowwise() %>%
          mutate(
            pct_det_ctrl = sum(!is.na(c_across(-c(name_sub_lib_id)))) / (ncol(.) - 1) * 100,
            PPM_mean = mean(c_across(-name_sub_lib_id), na.rm = TRUE),
            PPM_sd = sd(c_across(-name_sub_lib_id), na.rm = TRUE),
            PPM_ctrl = paste0(format(PPM_mean, digits = 3), " ± ", format(PPM_sd, digits = 3))
          ) %>%
          select(name_sub_lib_id, PPM_ctrl, pct_det_ctrl) %>%
          left_join(
            ppb_raw_ctrl %>%
              rowwise() %>%
              mutate(
                PPB_mean = mean(c_across(-name_sub_lib_id), na.rm = TRUE),
                PPB_sd = sd(c_across(-name_sub_lib_id), na.rm = TRUE),
                PPB_ctrl = paste0(format(PPB_mean, digits = 3), " ± ", format(PPB_sd, digits = 3))
              ) %>%
              select(name_sub_lib_id, PPB_ctrl),
            by = "name_sub_lib_id"
          ) %>%
          left_join(cas_key_2, by = "name_sub_lib_id") %>%
          select(cas, everything(), -c(Name))
    #- Merge the control and tumor PPM/PPB with the master table
      #_Create pared down filtered control and tumor for only variant diff
        joiner_ctrl <- ppm_ppb_summary_detected_ctrl %>%
          filter(name_sub_lib_id %in% diff_by_var) %>%
          arrange(name_sub_lib_id)
        joiner_tumor <- ppm_ppb_summary_detected %>%
          filter(name_sub_lib_id %in% diff_by_var) %>%
          select(-cas)
      #_Inner join these, use later below
        joiner <- joiner_ctrl %>%
          inner_join(joiner_tumor, by = "name_sub_lib_id") %>%
          select(-name_sub_lib_id)
  #+ All Inner join (only retain features in TUMORS that match controls) (ST3/4)
    #- Define the tumor variant columns for taking means
      ptc_cols <- names(ppm_raw)[grepl("^P\\d+$", names(ppm_raw))]
      ftc_cols <- names(ppm_raw)[grepl("^F\\d+$", names(ppm_raw))]
      fvptc_cols <- names(ppm_raw)[grepl("^FVPTC\\d+$", names(ppm_raw))]
    #- Bring in appropriate metadata to make prioritization decisions
      fragment_quality_info <- read_excel("Data and Metadata Files/primary_data.xlsx", sheet = "lib.subject.qsummary") %>%
        arrange(cas) %>%
        select(name_sub_lib_id, iMean) %>%
        rename(iMean_tumors = iMean)
    #- Fully join and gather metadata for detection and means
      SF4_ppm_fulltable <- ppm_raw_ctrl %>%
        inner_join(ppm_raw, by = "name_sub_lib_id") %>%
        left_join(cas_key_2, by = "name_sub_lib_id") %>%
        filter(!cas %in% endog_cas) %>%
          arrange(Name) %>%
          rowwise() %>%
          mutate(
            half_min = (min(c_across(where(is.numeric)), na.rm = TRUE)) / 2,
            max_value = max(c_across(where(is.numeric)), na.rm = TRUE)
          ) %>%
          mutate(
            pct_det_ctrl = sum(!is.na(c_across(T001:T009))) / length(c_across(T001:T009)),
            pct_det_tumor = sum(!is.na(c_across(F1:F20))) / length(c_across(F1:F20))
          ) %>% #Filter to ideal fragment based on pct_det_tumor, then ctrl, then iMean_tumors below
        left_join(fragment_quality_info, by = "name_sub_lib_id") %>%
        arrange(desc(pct_det_tumor), desc(pct_det_ctrl), desc(iMean_tumors)) %>%
        group_by(cas) %>%
        slice_head(n = 1) %>%
        ungroup() %>%
        select(-iMean_tumors) %>%
        mutate(across(
          T001:F20,
          ~ if_else(is.na(.x), half_min, .x)
        )) %>%
        rowwise() %>%
          mutate(
            mean_ctrl = mean(c_across(T001:T009), na.rm = TRUE),
            mean_tumor = mean(c_across(F1:F20), na.rm = TRUE),
            mean_PTC = mean(c_across(all_of(ptc_cols)), na.rm = TRUE),
            mean_FTC = mean(c_across(all_of(ftc_cols)), na.rm = TRUE),
            mean_FVPTC = mean(c_across(all_of(fvptc_cols)), na.rm = TRUE)
          ) %>%
        ungroup() %>%
        left_join(short_name %>% select(cas, annot_ident,IARC_Group), by = "cas") %>%
        select(Name, cas, annot_ident, IARC_Group, pct_det_ctrl, pct_det_tumor, half_min, max_value, name_sub_lib_id, mean_ctrl:mean_FVPTC,T001:F20) %>%
        arrange(cas)
    #- Format and export to excel
      format_ppb_value <- function(x) {
        if (is.na(x)) {
          return(NA_character_)
        }
        if (x < 1) {
          formatC(x, format = "e", digits = 1)
        } else if (x >= 1000) {
          formatC(x, format = "e", digits = 1)
        } else {
          round(x)
        }
      }
      PPB_version <- SF4_ppm_fulltable %>%
        select(Name:mean_tumor,mean_FTC,mean_FVPTC,mean_PTC,-c(T001:F20,name_sub_lib_id)) %>%
        mutate(across(half_min:mean_FVPTC, ~ .x * 1000)) %>%
          mutate(across(
            which(names(.) == "half_min"):ncol(.),
            ~ sapply(.x, format_ppb_value)))
      write_xlsx(PPB_version, "SF4_ppb_fulltable.xlsx")
    #- Create a pivoted version for subsequent analysis
      full_joiner <- SF4_ppm_fulltable %>%
        select(name_sub_lib_id,T001:F20) %>%
        pivot_longer(
          cols = -name_sub_lib_id,
          names_to = "sample_ID",
          values_to = "value"
        ) %>%
        mutate(
          variant = case_when(
            str_starts(sample_ID, "T") ~ "Ctrl",
            str_starts(sample_ID, "P") ~ "PTC",
            str_starts(sample_ID, "FVPTC") ~ "FV_PTC",
            str_starts(sample_ID, "F") ~ "FTC"
          )
        ) %>%
        pivot_wider(names_from = name_sub_lib_id, values_from = value) %>%
        mutate(
          tumor_vs_ctrl = ifelse(variant == "Ctrl", "Control", "Tumor"))
  #+ Create final master variant table (ST3/4)
    #- Create version of above SF4 ppm table except no filtering
      SF4_ppb_inclusive <- ppm_raw_ctrl %>%
        inner_join(ppm_raw, by = "name_sub_lib_id") %>%
        left_join(cas_key_2, by = "name_sub_lib_id") %>%
        filter(!cas %in% endog_cas) %>%
        arrange(Name) %>%
        rowwise() %>%
        mutate(
          half_min = (min(c_across(where(is.numeric)), na.rm = TRUE)) / 2,
          max_value = max(c_across(where(is.numeric)), na.rm = TRUE)
        ) %>%
        mutate(
          pct_det_ctrl = sum(!is.na(c_across(T001:T009))) / length(c_across(T001:T009)),
          pct_det_tumor = sum(!is.na(c_across(F1:F20))) / length(c_across(F1:F20))
        ) %>% 
        mutate(across(
          T001:F20,
          ~ if_else(is.na(.x), half_min, .x)
        )) %>%
        rowwise() %>%
        mutate(
          mean_ctrl_PPB = mean(c_across(T001:T009), na.rm = TRUE),
          mean_tumor_PPB = mean(c_across(F1:F20), na.rm = TRUE),
          mean_PTC_PPB = mean(c_across(all_of(ptc_cols)), na.rm = TRUE),
          mean_FTC_PPB = mean(c_across(all_of(ftc_cols)), na.rm = TRUE),
          mean_FVPTC_PPB = mean(c_across(all_of(fvptc_cols)), na.rm = TRUE)
        ) %>%
        ungroup() %>%
        arrange(cas) %>% 
        select(name_sub_lib_id,pct_det_ctrl,pct_det_tumor,mean_ctrl_PPB:mean_FVPTC_PPB) %>%
        mutate(across(mean_ctrl_PPB:mean_FVPTC_PPB, ~ .x * 1000)) # convert to ppb
    #- Join with master table, pare down to final version
      MT_final <- MTi %>%
        left_join(SF4_ppb_inclusive, by = "name_sub_lib_id") %>%
        select(short_name, cas, mode, annot_ident, p_value, highest, FTC:PTC, FTC_let:PTC_let, Carcinogenicity:GHS_var_diff_only, Potential_EDC, usage_class:Table_Class, pct_det_ctrl:mean_FVPTC_PPB)
    #- Pull superscript letters off and place on +/- columns
      #_Write function to do this
        move_superscript <- function(main, sd_raw) {
          superscript_pattern <- "[\\u1d43-\\u1d4d\\u02b0-\\u02b8\\u1d62-\\u1d6b\\u2070-\\u209f\\u2020-\\u2021]+"
          # Extract superscript from the main value
          superscript <- str_extract(main, superscript_pattern)
          # Remove superscript from main
          main_clean <- str_remove(main, superscript_pattern)
          # Format SD to 2 decimal places
          sd_clean <- format(round(as.numeric(sd_raw), 2), nsmall = 2)
          # Recombine
          paste0(main_clean, " ± ", sd_clean, ifelse(is.na(superscript), "", superscript))
        }
      #_Apply function to the relevant columns
        MT_export <- MT_final %>%
          mutate(
            FTC_let = case_when(
              mode == "quantitative" & !is.na(FTC_let) & str_detect(FTC, "±") ~
                move_superscript(FTC_let, str_trim(str_extract(FTC, "(?<=±).*"))),
              mode == "qualitative" & is.na(FTC_let) ~ FTC,
              TRUE ~ FTC_let
            ),
            FV_PTC_let = case_when(
              mode == "quantitative" & !is.na(FV_PTC_let) & str_detect(FV_PTC, "±") ~
                move_superscript(FV_PTC_let, str_trim(str_extract(FV_PTC, "(?<=±).*"))),
              mode == "qualitative" & is.na(FV_PTC_let) ~ FV_PTC,
              TRUE ~ FV_PTC_let
            ),
            PTC_let = case_when(
              mode == "quantitative" & !is.na(PTC_let) & str_detect(PTC, "±") ~
                move_superscript(PTC_let, str_trim(str_extract(PTC, "(?<=±).*"))),
              mode == "qualitative" & is.na(PTC_let) ~ PTC,
              TRUE ~ PTC_let
            )
          ) %>%
          select(-c(FTC:PTC)) %>%
          rename(FTC = FTC_let, FV_PTC = FV_PTC_let, PTC = PTC_let)
    #- Export
      write_xlsx(MT_export, "MT_final_SF4.xlsx")
  #+ Create subsetted table from master (Table 3)
    #- Bring in final info
      table_classes <- feature_metadata %>%
        select(cas, Table_Subclass, Table_Qualifier)
    #- Create table columns and format
      table_3 <- MT_final %>%
        #trouble populating OD-PABA so doing ident here manually
        mutate(
          annot_ident = if_else(cas == "21245-02-3", "Identification", annot_ident)
        ) %>%
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
          ) %>%
        left_join(table_classes, by = "cas") %>%
          mutate(
            Table_Class = ifelse(
              !is.na(Table_Subclass) & Table_Subclass != "",
              paste0(Table_Class, " (", Table_Subclass, ")"),
              Table_Class
            )
          ) %>%
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
        ) %>%
        mutate(short_name = str_replace(short_name, "NA$", "")) %>%
        select(name_sub_lib_id,short_name, cas,annot_ident, Table_Class, FTC_let:PTC_let, p_value) %>%
        # error with OD-PABA (identification) populating, so doing it manually
        mutate(
          short_name = if_else(cas == "21245-02-3", "OD-PABA", short_name),
          Table_Class = if_else(cas == "21245-02-3", "Organic UV Filter", Table_Class))
    #- Export
      write_xlsx(table_3, "table_3_export.xlsx")
  #+ Pull IARC 1 features (3D)
    #- Get the CAS numbers for IARC 1
      iarc_1 <- feature_metadata %>%
        filter(IARC_Group == "1") %>%
        pull(cas)
    #- Now pull the IARC1s in controls and tumors
      IARC_controls <- read_excel("Data and Metadata Files/primary_data.xlsx", sheet = "lib.subject.qsummary.cadaver") %>%
        filter(!is.na(comp)) %>%
        filter(pct_NA <= 0.5) %>%
        filter(cas %in% iarc_1) %>%
        select(pct_NA,name_sub_lib_id,iMean) %>%
        rename(pct_NA_ctrl = pct_NA, iMean_ctrl = iMean)
      IARC_tumors <- read_excel("Data and Metadata Files/primary_data.xlsx", sheet = "lib.subject.qsummary") %>%
        filter(cas %in% iarc_1) %>%
        filter(pct_NA <= 0.7) %>%
        select(cas, pct_NA, name, name_sub_lib_id,iMean) %>%
        arrange(name)
    #- Determine the matches, filter to 30% NA quant, pull features
      IARC_tumors_ctrl_filtered <- IARC_tumors %>%
        left_join(IARC_controls, by = "name_sub_lib_id") %>%
        select(cas, name, name_sub_lib_id,pct_NA,pct_NA_ctrl,iMean,iMean_ctrl) %>%
        filter(!is.na(pct_NA_ctrl)) %>%
        group_by(cas) %>%
        arrange(pct_NA, pct_NA_ctrl, desc(iMean)) %>%
        slice(1) %>%
        ungroup() %>%
        filter(pct_NA < 0.3, pct_NA_ctrl < 0.3) %>%
        pull(name_sub_lib_id)
    #- Now pull the actual features out, transpose for graphing
        IARC_comp <- full_joiner %>%
          select(sample_ID, any_of(IARC_tumors_ctrl_filtered)) %>%
          mutate(across(
            .cols = -c(sample_ID),
            .fns = ~ {
              col_min <- min(.x, na.rm = TRUE)
              replace(.x, is.na(.x), 0.5 * col_min)
            }
          )) %>%
          pivot_longer(-sample_ID, names_to = "name_sub_lib_id", values_to = "value") %>%
          pivot_wider(names_from = sample_ID, values_from = value) %>%
          left_join(cas_key_2, by = "name_sub_lib_id") %>%
          left_join(short_name %>% select(cas,annot_ident), by = "cas") %>%
          mutate(Name = if_else(
            annot_ident == "Annotation",
            paste0(Name, "\u2020"),
            Name
          )) %>%
          select(Name,T001:F20)
        write.csv(IARC_comp, "IARC_comp.csv")
    #! At this point, just exported and graphed these in Prism
    #- Run a ttest on log transformed data to determine signfiicant differences
      IARC_ttests <- full_joiner %>%
        select(sample_ID, any_of(IARC_tumors_ctrl_filtered)) %>%
        mutate(across(
          .cols = -c(sample_ID),
          .fns = ~ {
            col_min <- min(.x, na.rm = TRUE)
            replace(.x, is.na(.x), 0.5 * col_min)
          }
        )) %>%
        mutate(across(where(is.numeric), ~ log2(. + 1e-6))) %>%
        mutate(group = if_else(str_starts(sample_ID, "T00"), "Control", "Tumor")) %>%
        pivot_longer(cols = where(is.numeric) & !any_of("group"),
                    names_to = "chemical", values_to = "value") %>%
        group_by(chemical) %>%
        summarise(
          p_value = t.test(value ~ group)$p.value,
          mean_control = mean(value[group == "Control"]),
          mean_tumor = mean(value[group == "Tumor"]),
          .groups = "drop"
        ) %>%
        arrange(p_value)
        #! Data not shown but listed in manuscript results
  #+ Manual Spectral Validation (Supplementary)
    #- Create vector of variant significant CAS which are identifications
      sig_id <- table_3 %>%
      filter(annot_ident == "Identification") %>%
      mutate(id = str_extract(name_sub_lib_id, "(?<=_)[^_]+$")) %>%
      pull(id)
    #- Import the spectral validation data
      spectral_validation_i <- read_excel("Data and Metadata Files/primary_data.xlsx", sheet = "lib.subject.summary") %>%
        filter(id %in% sig_id) %>%
        mutate(
          mzMed = paste0(
            round(mzMed, 4)
          ),
          rtMed = paste0(
            round(rtMed, 2)
          )
        ) %>%
        select(
          name_sub_lib_id, lib, cas, name, short_display_name, id, subid, tmz,
          mzMed, trt, rtMed, libFragCnt, F1:F20
        ) %>%
        arrange(cas)
    #- Now, systematically determine the ideal sample to validate on
      #_Select candidate sample columns
        f1_index <- which(names(spectral_validation_i) == "F1")
        f20_index <- which(names(spectral_validation_i) == "F20")
        candidate_cols <- names(spectral_validation_i)[f1_index:f20_index]
      #_Summarise n and min for each candidate column per ID
        fragment_summary <- spectral_validation_i %>%
          group_by(id) %>%
          summarise(
            across(all_of(candidate_cols), list(
              n = ~ sum(!is.na(.)),
              min = ~ if (sum(!is.na(.)) >= 2) min(., na.rm = TRUE) else NA_real_
            ), .names = "{.col}__{.fn}"),
            .groups = "drop"
          )
      #_Reshape long and extract best columns by fragment count + min value
        long_summary <- fragment_summary %>%
          pivot_longer(-id, names_to = "measure", values_to = "value") %>%
          separate(measure, into = c("column", "stat"), sep = "__") %>%
          pivot_wider(names_from = stat, values_from = value)
        ideal_queries <- long_summary %>%
          filter(n >= 2 & !is.na(min)) %>%
          group_by(id) %>%
          arrange(desc(n), desc(min)) %>%
          summarise(
            ideal_query = first(column),
            ideal_query_backup = paste(head(column[-1], 5), collapse = ","),
            .groups = "drop"
          )
      #_Join back and finalize
        spectral_validation <- spectral_validation_i %>%
          left_join(ideal_queries, by = "id") %>%
          select(ideal_query, ideal_query_backup, id, everything()) %>%
          arrange(id)
    #- Now, pare this down and join with the feature library
      #_Pare down
        spectral_validation_pared <- spectral_validation %>%
          select(short_display_name,cas,ideal_query,ideal_query_backup,lib,id,subid,tmz,trt,mzMed,rtMed,libFragCnt) %>%
          arrange(id) %>%
          distinct()
      #_Get a sum of libfrags
      libfrag_summary <- spectral_validation_pared %>%
        group_by(id) %>%
        summarise(libFragCnt = sum(libFragCnt, na.rm = TRUE), .groups = "drop")
      #_Rejoin with feature library, export
        spectral_validation_wide <- spectral_validation_pared %>%
          rename(
            mz = mzMed,
            rt = rtMed
          ) %>%
          select(id, subid, mz, rt, ideal_query, ideal_query_backup) %>%
          pivot_wider(
            names_from = subid,
            values_from = c(mz, rt),
            names_glue = "Exp_{.value}{subid}"
          ) %>%
          relocate(id) %>%
          left_join(feature_lib, by = "id") %>%
          mutate(
            short_display_name = if_else(is.na(short_display_name), name, short_display_name)
          ) %>%
          left_join(libfrag_summary, by = "id") %>%
          select(
            short_display_name, cas, id, trt,
            starts_with("mz"),
            starts_with("Exp_mz"),
            starts_with("Exp_rt"),
            ideal_query, ideal_query_backup, libFragCnt
          )
      write_xlsx(spectral_validation_wide, "spectral_validation_wide.xlsx")
      #! At this point, exported to excel and formatted there