#* 5: Variant Visualization
#+ 5.1: Prepare quantitative (top 5) for visualization
#- 5.1.1: Pull the top 5 quant for visualization
top_5_quant <- summary_table |>
  arrange(p_value) |>
  slice_head(n = 5) |>
  select(name_sub_lib_id, p_value) |>
  pull(name_sub_lib_id)
#- 5.1.2: Extract CAS numbers for top 5 compounds, merge
column_renames <- tibble(name_sub_lib_id = top_5_quant) |>
  left_join(cas_key_2, by = c("name_sub_lib_id" = "name_sub_lib_id")) |>
  left_join(short_name, by = "cas") |>
  mutate(
    Short_display_name = ifelse(annot_ident == "Annotation",
      paste0(Short_display_name, "\u2020"),
      Short_display_name
    )
  ) |>
  select(name_sub_lib_id, Short_display_name) |>
  deframe()
#- 5.1.3: Subset z-scored quant data to the top 5,rename
tumors_quant_sig_top5 <- tumors_quant_sig |>
  select(variant, all_of(top_5_quant)) |>
  mutate(variant = factor(variant, levels = c("Follicular", "FV-PTC", "Papillary"))) |>
  arrange(variant) |>
  rename_with(~ column_renames[.x], .cols = all_of(names(column_renames))) |>
  rename_with(~ str_replace_all(.x, "\u2020", "\u1d43"))
#+ 5.2: Prepare qualitative data for visualization
#- 5.2.1: Filter to only the qualitative
qual_i <- quant_qual_results |>
  filter(mode == "qualitative") |>
  arrange(p_value) |>
  mutate(across(c(FTC, FV_PTC, PTC), ~ as.numeric(gsub("%", "", .)) / 100)) |>
  mutate(short_name = if_else(
    annot_ident == "Annotation",
    paste0(short_name, "\u2020"),
    short_name
  )) |>
  select(short_name, FTC:PTC)
#- 5.2.2: Make it a matrix
qual_matrix <- qual_i |>
  select(-short_name) |>
  as.matrix()
#- 5.2.3: Perform hierarchical clustering on rows
row_clusters <- hclust(dist(qual_matrix))
#- 5.2.4: Reorder the rows based on the clustering
qual_i_reordered <- qual_i |>
  slice(row_clusters$order) |> # multiply by 100
  mutate(across(c(FTC, FV_PTC, PTC), ~ round(.x * 100, 2))) |>
  mutate(order = row_number()) |>
  arrange(desc(order)) |>
  select(-order) |>
  mutate(short_name = str_replace_all(short_name, "\u2020", "\u1d43"))
#+ 5.3: Post-Hoc testing on quant
#- 5.3.1: Post-hoc on quants with CLD notation
posthoc_quant <- tumors_quant_sig |>
  mutate(variant = as.factor(variant)) |>
  mutate(variant = recode(variant,
    "Papillary" = "PTC",
    "Follicular" = "FTC",
    "FV-PTC" = "FV_PTC"
  ))
#- 5.3.2: Get the names for the quant data
posthoc_compound_names <- posthoc_quant |>
  select(where(is.numeric)) |>
  names()
#- 5.3.3: Main posthoc creation
posthoc_table <- dplyr::bind_rows(
  lapply(posthoc_compound_names, function(compound) {
    formula <- as.formula(paste0("`", compound, "` ~ variant"))
    aov_model <- aov(formula, data = posthoc_quant)
    p_value <- summary(aov_model)[[1]][["Pr(>F)"]][1]
    # Get group means
    means <- posthoc_quant |>
      group_by(variant) |>
      summarise(mean = mean(.data[[compound]], na.rm = TRUE), .groups = "drop") |>
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
  }))|> 
  left_join(cas_key_2, by = "name_sub_lib_id") |>
  select(cas, FTC, FV_PTC, PTC) |>
  rename(FTC_let = FTC, FV_PTC_let = FV_PTC, PTC_let = PTC)
#- 5.3.4: Separate superscripts from means
posthoc_decomp <- posthoc_table |>
  # Extract means and letters
  mutate(
    FTC = as.numeric(str_extract(FTC_let, "-?[0-9.]+")),
    FV_PTC = as.numeric(str_extract(FV_PTC_let, "-?[0-9.]+")),
    PTC = as.numeric(str_extract(PTC_let, "-?[0-9.]+")),
    FTC_L = convert_superscript(str_extract(FTC_let, "[ᵃᵇᶜᵈᵉᶠ]+")),
    FV_PTC_L = convert_superscript(str_extract(FV_PTC_let, "[ᵃᵇᶜᵈᵉᶠ]+")),
    PTC_L = convert_superscript(str_extract(PTC_let, "[ᵃᵇᶜᵈᵉᶠ]+"))
  ) |>
  select(cas, FTC:PTC_L) |>
  rowwise()
#- 5.3.5: Systematically determine the ranking
posthoc_ranked <- posthoc_decomp |>
rowwise() |>
mutate(
  ranking = {
    # Create a data frame for the current row
    groups <- tibble(
      group = c("FTC", "FV_PTC", "PTC"),
      mean = c(FTC, FV_PTC, PTC),
      letter = c(FTC_L, FV_PTC_L, PTC_L)
    )
    # Order by means (descending) and letters (ascending)
    ordered_groups <- groups |>
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
) |>
ungroup()
#- 5.3.6: Extract the highest group
highest_posthoc <- posthoc_ranked |>
  rowwise() |>
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
          strsplit(first_part, " = ")[[1]] |> paste(collapse = ", ")
        } else {
          first_part
        }
      }
    }
  ) |>
  select(cas, highest,ranking) |>
  ungroup()
#- 5.3.7: Join w/ original, add simple rank for qual
MTii <- quant_qual_results |>
  left_join(posthoc_table, by = "cas") |>
  left_join(highest_posthoc |> select(cas, highest), by = "cas") |>
  select(short_name, FTC:PTC, highest, everything()) |>
  # Update NA 'highest' for qualitative mode
  rowwise() |>
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
        ) |>
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
  ) |>
  ungroup()
#- 5.3.8: Compute actual p-values/stars for top 5 ANOVA
posthoc_table_pvalues <- dplyr::bind_rows(
  lapply(posthoc_compound_names, function(name_sub_lib_id) {
    formula <- as.formula(paste0("`", name_sub_lib_id, "` ~ variant"))
    aov_model <- aov(formula, data = posthoc_quant)
    p_value <- summary(aov_model)[[1]][["Pr(>F)"]][1]
    tibble(name_sub_lib_id = name_sub_lib_id, p_value = p_value)
  })
  ) |>
  arrange(p_value) |>
  rowwise() |>
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
  ) |>
  unnest(pairwise) |>
  select(name_sub_lib_id, p_value,comparison, star) |>
  left_join(cas_key_2, by = "name_sub_lib_id") |>
  select(-name_sub_lib_id) |>
  left_join(MTii |> select(cas,FTC_let:PTC_let), by = "cas") |>
  slice_head(n = 15)
#+ 5.4: Prepare Balloon Plot Data
#- 5.4.1: Organize data
balloon_data <- MTii |>
  filter(!is.na(highest)) |>
  select(cas, usage_class, highest) |>
  # Expand rows for double counting when multiple highest groups
  mutate(
    highest = strsplit(highest, ", ")
  ) |>
  unnest(highest) |>
  # Map the group names to the variant labels
  mutate(
    Variant = recode(highest,
      "FTC" = "Follicular",
      "FV_PTC" = "FV-PTC",
      "PTC" = "Papillary"
    )
  ) |>
  # Remove 'Equal' and NA values
  filter(Variant != "Equal", !is.na(Variant)) |>
  # Count occurrences for each variant within each usage class
  count(usage_class, Variant) |>
  mutate(
    usage_class = str_replace_all(usage_class, "[†‡]", "")
  )
#- 5.4.2: Determine the correct order for usage_class
usage_class_order <- balloon_data |>
  pivot_wider(names_from = Variant, values_from = n, values_fill = 0) |>
  arrange(Papillary, `FV-PTC`, Follicular) |>
  pull(usage_class)
#- 5.4.3: Apply the correct order for y-axis
balloon_data_graph <- balloon_data |>
  mutate(
    Variant = factor(Variant, levels = c("Follicular", "FV-PTC", "Papillary")),
    usage_class = factor(usage_class, levels = usage_class_order)
  )
#+ 5.5: Carcinogen class per GHS+IARC
#- 5.5.1: Apply classification function; add back subids
MTi <- MTii |>
  mutate(
    Carcinogenicity = mapply(
      classify_carcinogenicity,
      ghs = GHS_var_diff_only,
      iarc = IARC_Group
    ),
    subid = as.numeric(sub("^[^_]+_([0-9]+)_.*", "\\1", name_sub_lib_id))
  ) |>
  select(short_name, highest, Carcinogenicity, IARC_Group, GHS_var_diff_only, name_sub_lib_id, subid, everything()) |>
  arrange(Carcinogenicity)
#- 5.5.2: Leftjoin and summarize
carc_by_variant <- MTi |>
  # Expand rows for double counting when multiple highest groups
  mutate(
    highest = strsplit(highest, ", ")
  ) |>
  unnest(highest) |>
  # Map the group names to the variant labels
  mutate(
    Variant = recode(highest,
      "FTC" = "Follicular",
      "FV_PTC" = "FV-PTC",
      "PTC" = "Papillary"
    )
  ) |>
  # Count occurrences for each Variant and Carcinogenicity
  count(Variant, Carcinogenicity) |>
  # Pivot to wide format
  pivot_wider(
    names_from = Carcinogenicity,
    values_from = n,
    values_fill = 0 # Fill missing counts with 0
  ) |>
  # Reorder columns for readability
  select(
    Variant,
    "Known Carcinogen",
    "Likely Carcinogen",
    "Possible Carcinogen",
    "Uncertain Risk",
    "Unclassified",
    everything()
  ) |>
  filter(Variant != "Equal")
#+ 5.5: Create all variant graphs
#- 5.5.1: Bar plot for top 5 quant
# p3A <- plot_top5_quant(tumors_quant_sig_top5)
#- 5.5.2: Heatmap for qualitative features
# p3B <- plot_qualitative_heatmap(qual_i_reordered)
#- 5.5.3: Balloon plot for usage class vs variant
p3C <- plot_balloon(balloon_data_graph)
