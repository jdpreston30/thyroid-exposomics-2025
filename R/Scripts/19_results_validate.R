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
#+ 19.3: Methods
.sec("Methods")
.chk("features tested, total (1,476)",       1476, nrow(anova_all) + nrow(fisher_all))
.chk("  quantitative features (435)",        435,  nrow(anova_all))
.chk("  qualitative features (1,041)",       1041, nrow(fisher_all))
.chk("sent to manual validation (62)",       62,   nrow(quant_qual_results))
.chk("retained after validation (29)",       29,   nrow(ancova_summary))
#+ 19.4: Results 3.2-3.3 -- Annotation Landscape
.sec("Results 3.2-3.3 - annotation landscape")
#! Table_Class, not Graph_Class: Graph_Class collapses small-n categories for figure legibility and yields 23, while the usage-class taxonomy the prose describes is the 26 Class-level rows in ST2.
#! detection_no_endog is built from tumor_raw, which already has drop_excluded() applied at 00c, so this count is post-exclusion and the pipeline value is correct by construction. It is exactly the check that would catch a future exclusion silently shifting the median.
.chk("median annotated per sample (303)", 303, median(detection_no_endog$total_detected, na.rm = TRUE))
.chk("usage classes (26)",                   26,   dplyr::n_distinct(feature_metadata$Table_Class))
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
#+ 19.5: Results 3.4 -- Exposome-wide Association
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
#+ 19.6: Results 3.5 -- Robustness
.sec("Results 3.5 - robustness")
.chk("significant under both adjustments (22/29)", 22, sum(ancova_summary$survives_year, na.rm = TRUE))
.chk("significant, binned collection year (14/29)", 14, sum(ancova_summary$survives_binned, na.rm = TRUE))
#+ 19.7: Verdict
cat(sprintf("\n  %s   mismatches: %d   not-evaluable: %d\n",
            if (.vfail == 0L) "ALL EVALUABLE CLAIMS MATCH." else "*** REVIEW NEEDED ***",
            .vfail, .vskip))
cat("=============================================================\n\n")
rm(.chk, .sec, .numify, .vfail, .vskip)
