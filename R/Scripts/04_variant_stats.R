#* 4: Variant Analysis Stats
#+ 4.1: ANOVA Stats (Quant)
#- 4.1.1: Full p-value vector for the FDR complementary analysis (reviewer C3)
#! Standalone on purpose: nothing downstream consumes it, so adding it cannot disturb summary_table_i, fisher_results_i, MT_final, Table 3 or ST6. Must be computed BEFORE the p < 0.05 filter -- p.adjust() needs all m p-values, and adjusting only the survivors would be anti-conservative and wrong. Adjustment happens in 4.5.
anova_all <- tumors_quant_wt |>
  pivot_longer(-variant, names_to = "name_sub_lib_id", values_to = "value") |>
  group_by(name_sub_lib_id) |>
  summarise(
    p_value = broom::tidy(aov(value ~ variant))$p.value[1],
    .groups = "drop"
  ) |>
  arrange(p_value)
#- 4.1.2: Run ANOVA, pull significant features
anova_results_sig <- tumors_quant_wt |>
  pivot_longer(-variant, names_to = "name_sub_lib_id", values_to = "value") |>
  group_by(name_sub_lib_id) |>
  summarise(
    p_value = broom::tidy(aov(value ~ variant))$p.value[1] # Extract p-value from ANOVA
  ) |>
  filter(p_value < 0.05) |>
  arrange(p_value) |>
  arrange(name_sub_lib_id) |>
  pull(name_sub_lib_id) # Extract significant feature names
#- 4.1.3: Filter tumors_quant to keep only significant compounds and apply z-score
tumors_quant_sig_i <- tumors_quant_wt |>
  select(variant, all_of(anova_results_sig)) |>
  mutate(across(-variant, ~ scale(.)[, 1]))
#- 4.1.4: Make a second cas key
cas_key_2 <- tumor_raw |>
  select(name_sub_lib_id, cas, short_display_name) |>
  rename(Name = short_display_name)
#- 4.1.5: Create summary table with reran ANOVA
summary_table_i <- tumors_quant_sig_i |>
  pivot_longer(-variant, names_to = "name_sub_lib_id", values_to = "value") |>
  group_by(name_sub_lib_id, variant) |>
  mutate(variant = if_else(variant == "FV-PTC", "FV_PTC", variant)) |>
  summarise(
    mean_sd = paste0(round(mean(value, na.rm = TRUE), 2), " ± ", round(sd(value, na.rm = TRUE), 2)),
    .groups = "drop"
  ) |>
  pivot_wider(names_from = variant, values_from = mean_sd) |>
  left_join(
    tumors_quant_sig_i |>
      pivot_longer(-variant, names_to = "name_sub_lib_id", values_to = "value") |>
      group_by(name_sub_lib_id) |>
      summarise(
        #! Full precision retained: pre-rounding to 4 dp put values like MEHP (0.009533) exactly on a 3 dp boundary, so Table 3 and ST6 rendered the same P as 0.009 vs 0.010
        p_value = summary(aov(value ~ variant))[[1]][["Pr(>F)"]][1],
        .groups = "drop"
      ),
    by = "name_sub_lib_id"
  ) |>
  select(name_sub_lib_id, Follicular, FV_PTC, Papillary, p_value) |>
  arrange(p_value) |>
  mutate(mode = "quantitative") |>
  left_join(cas_key_2, by = "name_sub_lib_id") |>
  select(Name, cas, name_sub_lib_id, mode, p_value, everything()) |>
  arrange(p_value)
#+ 4.2: Fisher's Stats (Qual)
#- 4.2.1: Full p-value vector for the FDR complementary analysis (reviewer C3)
#! Same rationale as 4.1.0 -- computed before the p < 0.05 filter, and consumed by nothing. Raw p-values only -- the adjustment is pooled across both modes in section 4.5.
fisher_all <- tumors_qual |>
  pivot_longer(-variant, names_to = "name_sub_lib_id", values_to = "detected") |>
  group_by(name_sub_lib_id) |>
  summarise(
    p_value = tryCatch(
      fisher.test(table(variant, detected))$p.value,
      error = function(e) NA_real_
    ),
    .groups = "drop"
  ) |>
  arrange(p_value)
#- 4.2.2: Significant qualitative features
fisher_results_i <- tumors_qual |>
  pivot_longer(-variant, names_to = "name_sub_lib_id", values_to = "detected") |>
  group_by(name_sub_lib_id) |>
  summarise(
    p_value = tryCatch(
      fisher.test(table(variant, detected))$p.value,
      error = function(e) NA_real_
    ),
    Follicular = (sum(detected[variant == "Follicular"]) / sum(variant == "Follicular")) * 100,
    Papillary = (sum(detected[variant == "Papillary"]) / sum(variant == "Papillary")) * 100,
    FV_PTC = (sum(detected[variant == "FV-PTC"]) / sum(variant == "FV-PTC")) * 100
  ) |>
  arrange(p_value) |>
  filter(p_value < 0.05) |>
  arrange(name_sub_lib_id) |>
  arrange(p_value) |>
  left_join(cas_key_2, by = "name_sub_lib_id") |>
  mutate(mode = "qualitative") |>
  select(Name, everything())
#+ 4.3: Quality Control
#- 4.3.1: Visually inspect for duplicates (with different fragments) between the ANOVA and Fisher's
#_Bind names and modes to check for duplicates
summary_table_i_dup <- summary_table_i |>
  select(Name, cas, name_sub_lib_id, mode)
fisher_results_dup <- fisher_results_i |>
  select(Name, cas, name_sub_lib_id, mode)
#_Add duplicate check column
dupl_check <- rbind(summary_table_i_dup, fisher_results_dup) |>
  arrange(cas) |>
    mutate(
      duplicate = case_when(
        Name %in% Name[base::duplicated(Name)] |
          cas %in% cas[base::duplicated(cas)] ~ "*",
        TRUE ~ ""
      )
    ) |>
  arrange(Name) |>
  arrange(desc(duplicate))
#_Inspect          
print(dupl_check, n = Inf)
#- 4.3.2: Duplicates with quant and qual fragments (quant chosen)
dupl_check_qual_quant <- dupl_check |>
  group_by(Name, cas) |>
  filter(n() == 2, sum(mode == "qualitative") == 1, sum(mode == "quantitative") == 1) |>
  ungroup()
#! Benz(a)anthracene_1_BP2.GC2_CP2220 (quant) chosen over Benz(a)anthracene_0_BP3.GC2_CP3027
#! 2-amino-4-nitrotoluene_1_BP3.GC2_CP3003 (quant) chosen over 2-amino-4-nitrotoluene_0_BP3.GC2_CP3003
#- 4.3.3: Duplicates with multiple qual frags (highest detection chosen)
#_Pull these features
dupl_check_qual <- dupl_check |>
  group_by(Name, cas) |>
  filter(n() > 1, all(mode == "qualitative")) |>
  ungroup()
#_Subset to only those, then systematically choose lowest detection fragment to remove
qual_only_ids <- dupl_check_qual$name_sub_lib_id
qual_dup_removed_frags <- fisher_results_i |>
  filter(name_sub_lib_id %in% qual_only_ids) |>
  mutate(sum_det = Follicular + Papillary + FV_PTC) |>
  select(Name, cas, name_sub_lib_id, sum_det) |>
  group_by(Name) |>
  slice_min(sum_det, with_ties = FALSE) |>
  ungroup()
print(qual_dup_removed_frags, n = Inf)
#! Aldicarb: aldicarb_0_BP2.GC2_CP2507 > detection than aldicarb_1_BP2.GC2_CP2507
#! Fenvalerate: Fenvalerate_1_3_BP3.GC2_CP3164 > detection than Fenvalerate_2_3_BP3.GC2_CP3165
#! Octyl-dimethyl-PABA: Octyl-dimethyl-PABA_1_BP2.GC2_CP2331 > detection than Octyl-dimethyl-PABA_2_BP2.GC2_CP2331
#! tris(tribromoneopentyl) (CP2302) removed upstream -- see 00c_FTs.R section 0c.1.20
#- 4.3.4: Pull the duplicate fragment names for removal
#_qual/quant duplicates
qual_quant_dupl_removed <- dupl_check_qual_quant |>
  filter(mode == "qualitative") |>
  pull(name_sub_lib_id)
#_qual/qual duplicates
qual_qual_dupl_removed <- qual_dup_removed_frags |>
  pull(name_sub_lib_id)
#- 4.3.5: Remove duplicate fragments with lower detection or qual when a quant is available:
qual_single_frag <- fisher_results_i |>
  filter(!name_sub_lib_id %in% qual_qual_dupl_removed) |>
  filter(!name_sub_lib_id %in% qual_quant_dupl_removed) |>
  select(cas,name_sub_lib_id,mode,everything()) |>
  arrange(Name)
print(qual_single_frag, n = Inf)
#+ 4.4: Compile Final Results
#- 4.4.1: Pull short name from feature metadata
short_name <- read_excel(config$paths$chemical_metadata, sheet = "feature_metadata") |>
  select(cas, name, Potential_EDC, IARC_Group, GHS_var_diff_only,Short_display_name, Graph_Class, Superclass, Table_Class, Annotation_or_Identification) |>
  rename(annot_ident = Annotation_or_Identification) |>
  #! re-read from source, so re-apply the 00c exclusion; harmless today (left-join target only) but keeps the guarantee structural
  drop_excluded() |>
  #! Display-casing correction pending a source fix in the feature_metadata sheet; root name is Title-Case in tables per convention
  mutate(Short_display_name = str_replace(Short_display_name, "^o-aminoazotoluene", "o-Aminoazotoluene")) |>
#! Table 3 and ST6 take their names from THIS object, not from ST1_import, so the shared nomenclature fixes must be re-applied on this read as well.
  mutate(across(c(name, Short_display_name), fix_chem_nomenclature)) |>
  #! Same reason as drop_excluded() above: this is a second read of the raw sheet, so the 00c "Dye intermediates" -> "Dye Intermediates" normalisation must be re-applied here too. Table_Class reaches Table 3 through THIS object, not through feature_metadata, so without it the 00c fix never lands.
  mutate(across(any_of(c("Table_Class", "Graph_Class")),
                \(x) str_replace(x, "^Dye intermediates$", "Dye Intermediates")))
#- 4.4.2: Bind with quant features
#_Filter quant to be safe
summary_table <- summary_table_i |>
  filter(!name_sub_lib_id %in% qual_qual_dupl_removed) |>
  filter(!name_sub_lib_id %in% qual_quant_dupl_removed)
#_Bind and name
quant_qual_results <- rbind(qual_single_frag, summary_table) |>
  left_join(short_name, by = "cas") |>
  rename(usage_class = Graph_Class) |>
  select(cas, name, Short_display_name, name_sub_lib_id, mode, p_value, Follicular, FV_PTC, Papillary, everything()) |>
  rename(FTC = Follicular, PTC = Papillary) |>
  mutate(across(c(FTC, FV_PTC, PTC),
    ~ if_else(mode == "qualitative", paste0(.x, "%"), as.character(.x)),
    .names = "{.col}"
  )) |>
  rename(short_name = Short_display_name)
#+ 4.5: FDR Complementary Analysis (reviewer C3)
#- 4.5.1: Bind together all ANOVA and Fisher's p-values
#! Pooled across BOTH modes, not adjusted per mode. Table 3 reports the quantitative and qualitative findings as a single set, so the correction family is the full set of tests; controlling each mode separately and then reporting the union gives no FDR guarantee on the combined set, and is the more permissive choice. Pooling is also strictly more conservative (larger m tightens every threshold).
fdr_all <- bind_rows(
  anova_all  |> mutate(mode = "quantitative"),
  fisher_all |> mutate(mode = "qualitative")
) |>
  mutate(q_BH = p.adjust(p_value, method = "BH")) |>
  arrange(p_value)
#- 4.5.2: Emit every number the manuscript claims, so the prose can be checked against the pipeline
{
  .m    <- sum(!is.na(fdr_all$p_value))
  .sig  <- sum(fdr_all$p_value < 0.05, na.rm = TRUE)
  .surv <- sum(fdr_all$q_BH   < 0.05, na.rm = TRUE)
  cat("\n-- FDR complementary analysis (pooled Benjamini-Hochberg) --\n")
  cat(sprintf("  features tested              : %s  (quant %d + qual %d)\n",
              format(.m, big.mark = ","), nrow(anova_all), nrow(fisher_all)))
  cat(sprintf("  nominally significant P<0.05 : %d  (quant %d + qual %d)\n", .sig,
              sum(anova_all$p_value < 0.05, na.rm = TRUE),
              sum(fisher_all$p_value < 0.05, na.rm = TRUE)))
  cat(sprintf("  expected by chance at P<0.05 : %.1f\n", .m * 0.05))
  cat(sprintf("  smallest P                   : %.3g\n", min(fdr_all$p_value, na.rm = TRUE)))
  cat(sprintf("  surviving pooled BH q<0.05   : %d\n", .surv))
  cat(sprintf("  unique chemicals after dedup : %d\n", nrow(quant_qual_results)))
  if (.surv > 0) print(as.data.frame(fdr_all |> filter(q_BH < 0.05)), row.names = FALSE)
  rm(.m, .sig, .surv)
}
