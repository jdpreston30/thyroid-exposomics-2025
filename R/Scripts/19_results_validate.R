#* 19: Numerical Claim Validation
#! Emits every numeric claim made in the manuscript prose, in the order it appears, next to the value the pipeline actually produces. PROSE values are hardcoded from the submitted text on purpose -- that is what makes drift detectable; when a claim legitimately changes, update the PROSE argument here in the same commit as the manuscript edit. Read-only: derives nothing new and writes no files.
#+ 19.1: Reporting Helpers
#! Every getter is wrapped, so a partial pipeline degrades to "??" rather than aborting the run and losing the session-info step that follows. On any non-PASS the expression itself is echoed, so a wrong object or column name is immediately diagnosable instead of looking like a failed claim.
.vfail <- 0L; .vskip <- 0L
.numify <- function(x) suppressWarnings(as.numeric(gsub("[^0-9.eE+-]", "", as.character(x))))
.chk <- function(label, prose, expr, tol = 0) {
  src <- paste(deparse(substitute(expr)), collapse = " ")
  got <- tryCatch(suppressWarnings(expr), error = \(e) NULL, warning = \(w) NULL)
  if (is.null(got) || length(got) != 1L || all(is.na(got))) {
    .vskip <<- .vskip + 1L
    cat(sprintf("  [ ?? ] %-50s prose=%-11s  <-- could not evaluate: %s\n", label, prose, src))
    return(invisible(NULL))
  }
  p <- .numify(prose); q <- .numify(got)
  ok <- if (is.na(p) || is.na(q)) identical(as.character(prose), as.character(got)) else abs(p - q) <= tol
  if (!ok) .vfail <<- .vfail + 1L
  cat(sprintf("  [%s] %-50s prose=%-11s pipeline=%s%s\n",
              if (ok) "PASS" else "FAIL", label, prose, got,
              if (ok) "" else sprintf("   <-- %s", src)))
}
.sec <- function(x) cat(sprintf("\n-- %s --\n", x))
cat("\n================= NUMERICAL CLAIM VALIDATION =================\n")
cat("PROSE = as written in the submitted manuscript | PIPELINE = computed now\n")
#+ 19.2: Abstract
.sec("Abstract")
.chk("samples per type (n = 20)",            20,   min(table(tumors_qual$variant)))
.chk("total DTC samples (N = 60)",           60,   nrow(tumors_qual))
.chk("cadaver thyroids (n = 8)",             8,    nrow(cadaver_metadata))
.chk("library screened (709 compounds)",     709,  dplyr::n_distinct(ST1_import$cas))
.chk("chemicals annotated (442)",            442,  nrow(feature_metadata))
.chk("potential EDCs, abstract (40%)",       40,
     round(100 * sum(feature_metadata$Potential_EDC == "Y", na.rm = TRUE) / nrow(feature_metadata)))
.chk("chemicals differing between types (29)", 29, nrow(ancova_summary))
#! The abstract's "Two IARC Group 1 carcinogens" is an output of the >= 0.7 detection filter at 10_post_validation_clean.R:54, not of significance -- exactly the kind of claim that would drift silently if that threshold ever moved.
.chk("IARC Group 1 elevated in tumors (2)", 2, dplyr::n_distinct(IARC_combined$short_name))
#+ 19.3: Methods
.sec("Methods")
.chk("features tested, total (1,476)",       1476, nrow(anova_all) + nrow(fisher_all))
.chk("  quantitative features (435)",        435,  nrow(anova_all))
.chk("  qualitative features (1,041)",       1041, nrow(fisher_all))
.chk("sent to manual validation (62)",       62,   nrow(quant_qual_results))
.chk("retained after validation (29)",       29,   nrow(ancova_summary))
#! lib_composition is the standalone pre-filter read added at 00c_FTs.R section 0c.1.22; ST1_import cannot answer these because it is already filtered on Disposition != "Endogenous".
.chk("standards in internal library (737)", 737, sum(lib_composition$n_chemicals))
.chk("  endogenous metabolites (28)",       28,
     lib_composition$n_chemicals[lib_composition$Disposition == "Endogenous"])
.chk("  exogenous chemicals (699)",         699,
     lib_composition$n_chemicals[lib_composition$Disposition == "Exogenous"])
.chk("  mixed-origin chemicals (10)",       10,
     lib_composition$n_chemicals[lib_composition$Disposition == "Exogenous and Endogenous"])
#! The 23 are the annotated endogenous chemicals dropped at 00c section 0c.4.2, a different set from the 28 purely endogenous library standards above.
.chk("endogenous chemicals removed (23)",   23,   dplyr::n_distinct(endog_cas))
#! feature_metadata, not ST1_import: the 10 in lib_composition are the dual-disposition entries in the LIBRARY, while section 2.4's 7 are the ones actually annotated in samples. Two different sets, and only the annotated one belongs to this claim.
.chk("dual-disposition chemicals retained (7)", 7,
     sum(feature_metadata$Disposition == "Exogenous and Endogenous", na.rm = TRUE))
#! The IARC validation set holds 10 distinct chemicals, but section 2.7 says "an ADDITIONAL n ... beyond the 62", and o-toluidine sits in both sets, so the additional count is 9 and 62 + 9 = 71 chemicals reached manual validation. Three traps here: iv_wide is keyed on library id so one chemical spans several rows (pentachlorophenol is CP1016 and CP2242); validated_iarc is the set that PASSED, not the set submitted; and short_name carries a trailing * for Level 2 ids in some objects, so it must be stripped or the overlap silently under-counts. The PCB exclusion the prose names is already applied at 07_validation_prep.R.
#! Prose was 11 as submitted -- reproducible as neither 10 nor 9, so it was corrected to 9 on 2026-09-03.
.chk("additional IARC Group 1 validated (9)", 9,
     sum(!unique(sub("\\*$", "", iv_wide$short_name)) %in%
           unique(sub("\\*$", "", quant_qual_results$short_name))))
.chk("4-aminobiphenyl fragments (3)", 3, sum(grepl("^4-aminobiphenyl_", IARC_namesub_pull)))
#+ 19.4: Results 3.1 -- Demographics
.sec("Results 3.1 - demographics")
.chk("mean age (60.5 y)",  60.5, round(mean(clinical_data$Age, na.rm = TRUE), 1), tol = 0.05)
.chk("female (68%)",       68,   round(mean(clinical_data$Sex == "Female", na.rm = TRUE) * 100))
.chk("  Stage I (63%)",    63,   round(mean(clinical_data$Pathologic_Stage == "I")   * 100))
.chk("  Stage II (25%)",   25,   round(mean(clinical_data$Pathologic_Stage == "II")  * 100))
.chk("  Stage III (8%)",   8,    round(mean(clinical_data$Pathologic_Stage == "III") * 100))
.chk("  Stage IV (3%)",    3,    round(mean(clinical_data$Pathologic_Stage == "IV")  * 100))
#+ 19.5: Results 3.2-3.3 -- Annotation Landscape
.sec("Results 3.2-3.3 - annotation landscape")
#! Table_Class, not Graph_Class: Graph_Class collapses small-n categories for figure legibility and yields 23, while the usage-class taxonomy the prose describes is the 26 Class-level rows in ST2.
#! detection_no_endog is built from tumor_raw, which already has drop_excluded() applied at 00c, so this count is post-exclusion and the pipeline value is correct by construction. It is exactly the check that would catch a future exclusion silently shifting the median.
#! Both updated 2026-09-03 from 303 / 0.397. Those were stale from the TTBNP (CP2302) exclusion of 2026-08-03: dropping one widely detected chemical lowered the per-sample median by one and shifted the Kruskal-Wallis test computed on those same counts. Figure 1B had already self-corrected because plot_detection_scatter() is passed kw_result directly, so the figure and the prose had disagreed since August. Caught by this script on its first run.
.chk("median annotated per sample (302)", 302, median(detection_no_endog$total_detected, na.rm = TRUE))
.chk("annotated per type, Kruskal-Wallis P (0.413)", 0.413, round(kw_result, 3), tol = 0.0005)
.chk("usage classes (26)",                   26,   dplyr::n_distinct(feature_metadata$Table_Class))
#! Table_Header is the three-way broad grouping that Table 2 and Figure 2A key off; the prose counts and percentages come straight from it.
.chk("  pollutants and industrial (222)", 222, sum(feature_metadata$Table_Header == "Pollutants and Industrial Chemicals", na.rm = TRUE))
.chk("  pollutants and industrial (50%)", 50,  round(100 * mean(feature_metadata$Table_Header == "Pollutants and Industrial Chemicals", na.rm = TRUE)))
.chk("  agrochemicals (161)",             161, sum(feature_metadata$Table_Header == "Agrochemicals", na.rm = TRUE))
.chk("  agrochemicals (36%)",             36,  round(100 * mean(feature_metadata$Table_Header == "Agrochemicals", na.rm = TRUE)))
.chk("  other chemicals (59)",            59,  sum(feature_metadata$Table_Header == "Other Chemicals", na.rm = TRUE))
.chk("  other chemicals (13%)",           13,  round(100 * mean(feature_metadata$Table_Header == "Other Chemicals", na.rm = TRUE)))
.chk("superclasses (10)", 10, dplyr::n_distinct(feature_metadata$Superclass, na.rm = TRUE))
.chk("classes (59)",      59, dplyr::n_distinct(feature_metadata$Class,      na.rm = TRUE))
#! superclass_pct / class_pct are built in 03_classes.R and carry pct to one decimal; the prose rounds to whole percent.
.chk("  benzenoids (53%)",               53, round(superclass_pct$pct[superclass_pct$Superclass == "Benzenoids"]))
.chk("  organoheterocyclic (14%)",       14, round(superclass_pct$pct[superclass_pct$Superclass == "Organoheterocyclic compounds"]))
.chk("  organic acids and deriv. (11%)", 11, round(superclass_pct$pct[superclass_pct$Superclass == "Organic acids and derivatives"]))
.chk("  benzene and subst. deriv. (41%)",41, round(class_pct$pct[class_pct$Class == "Benzene and substituted derivatives"]))
.chk("  fatty acyls (7%)",                7, round(class_pct$pct[class_pct$Class == "Fatty Acyls"]))
.chk("  organic phosphoric acids (4%)",   4, round(class_pct$pct[class_pct$Class == "Organic phosphoric acids and derivatives"]))
.chk("IARC-evaluated, count (133)",          133,
     sum(!is.na(feature_metadata$IARC_Group) & feature_metadata$IARC_Group != ""))
.chk("IARC-evaluated, percent (30.1%)",      30.1,
     round(100 * sum(!is.na(feature_metadata$IARC_Group) & feature_metadata$IARC_Group != "") /
             nrow(feature_metadata), 1), tol = 0.05)
.chk("  IARC Group 1 (23)",                  23,   sum(feature_metadata$IARC_Group == "1",  na.rm = TRUE))
.chk("  IARC Group 2A (13)",                 13,   sum(feature_metadata$IARC_Group == "2A", na.rm = TRUE))
.chk("  IARC Group 2B (52)",                 52,   sum(feature_metadata$IARC_Group == "2B", na.rm = TRUE))
.chk("  IARC Group 3 (45)",                  45,   sum(feature_metadata$IARC_Group == "3",  na.rm = TRUE))
.chk("potential EDCs, Discussion (40.5%)",   40.5,
     round(100 * sum(feature_metadata$Potential_EDC == "Y", na.rm = TRUE) / nrow(feature_metadata), 1),
     tol = 0.05)
#+ 19.6: Results 3.4 -- Exposome-wide Association
.sec("Results 3.4 - exposome-wide association")
.chk("statistical differences, chemicals (62)", 62, nrow(quant_qual_results))
.chk("  quantitative of those (33)",         33,   sum(quant_qual_results$mode == "quantitative"))
.chk("  qualitative of those (29)",          29,   sum(quant_qual_results$mode == "qualitative"))
.chk("validated chemicals (29)",             29,   nrow(ancova_summary))
.chk("  quantitative validated (17)",        17,   sum(ancova_summary$mode == "quant"))
.chk("  qualitative validated (12)",         12,   sum(ancova_summary$mode == "qual"))
.chk("nominally significant features (67)",  67,   sum(fdr_all$p_value < 0.05, na.rm = TRUE))
.chk("expected by chance, P<0.05 (~74)",     74,   round(sum(!is.na(fdr_all$p_value)) * 0.05), tol = 1)
.chk("surviving pooled BH (0)",              0,    sum(fdr_all$q_BH < 0.05, na.rm = TRUE))
.chk("  Level 1 identifications (25)", 25, sum(MT_final$quality == 1, na.rm = TRUE))
.chk("  Level 2 identifications (4)",   4, sum(MT_final$quality == 2, na.rm = TRUE))
#! anova_all carries only name_sub_lib_id and p_value, so DEET is located by name rather than joined.
.chk("DEET, P (0.039)", 0.039,
     round(anova_all$p_value[grepl("DEET|toluamide", anova_all$name_sub_lib_id, ignore.case = TRUE)][1], 3),
     tol = 0.0005)
#! carc_summary is already counted by Variant, so the distinct-chemical figure is re-derived from MTi here.
.chk("carcinogenic-evidence chemicals (10)", 10,
     dplyr::n_distinct(MTi$short_name[MTi$cas %in% MT_final_cas_list & MTi$Carcinogenicity != "Unclassified"]))
#+ 19.7: Results 3.5 -- Robustness
.sec("Results 3.5 - robustness")
.chk("significant under both adjustments (22/29)", 22, sum(ancova_summary$survives_year, na.rm = TRUE))
.chk("significant, binned collection year (14/29)", 14, sum(ancova_summary$survives_binned, na.rm = TRUE))
#+ 19.8: Results 3.6 -- Quantitative Estimates
.sec("Results 3.6 - quantitative estimates")
#! The prose claims significant elevation without printing a value, so the claim under test is the
#! threshold itself. Table 4's literature columns are transcribed from cited papers and are external
#! to the pipeline -- only the study-side ppb means are derived here, and those are covered by Table 4's
#! own build. A TRUE/FALSE prose value falls through .numify() to the string comparison branch.
.chk("o-toluidine elevated in tumors (P < 0.05)",     TRUE, toluidine_p < 0.05)
.chk("4-aminobiphenyl elevated in tumors (P < 0.05)", TRUE, aminobiphenyl_0_p < 0.05)
#+ 19.9: Discussion
.sec("Discussion")
.chk("phthalate species annotated (~20)", 20,
     sum(grepl("phthalate", feature_metadata$name, ignore.case = TRUE)), tol = 2)
.chk("not evaluated as carcinogens (~70%)", 70,
     round(100 - 100 * sum(!is.na(feature_metadata$IARC_Group) & feature_metadata$IARC_Group != "") /
             nrow(feature_metadata)), tol = 1)
#+ 19.10: Verdict
cat(sprintf("\n  %s   mismatches: %d   not-evaluable: %d\n",
            if (.vfail == 0L) "ALL EVALUABLE CLAIMS MATCH." else "*** REVIEW NEEDED ***",
            .vfail, .vskip))
cat("=============================================================\n\n")
rm(.chk, .sec, .numify, .vfail, .vskip)
