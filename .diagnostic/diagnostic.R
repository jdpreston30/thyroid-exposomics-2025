#* 1: Data Landscape Probe — structures we need for the confounding investigation
#! Run AFTER the pipeline (through ~script 12) so the objects exist. For each object:
#! class, dims, column names, key/join columns, and a small preview. Goal: figure out
#! how to attach collection YEAR + tumor TYPE onto the quant/qual chemical data so we
#! can run covariate screens + targeted ANCOVA on the headline chemicals.
#! (Covariate-balance code now lives in covariates.R; this file is the structure probe.)

#+ 1.0: Describe helper — prints structure without flooding on wide chemical tables
desc <- function(nm) {
  if (!exists(nm, envir = .GlobalEnv)) { cat(sprintf("  [absent]  %s\n", nm)); return(invisible()) }
  x <- get(nm, envir = .GlobalEnv)
  cl <- paste(class(x), collapse = "/")
  if (is.data.frame(x)) {
    cat(sprintf("\n▸ %-26s <%s>  %d rows x %d cols\n", nm, cl, nrow(x), ncol(x)))
    cn <- names(x)
    show <- if (length(cn) > 16) c(cn[1:10], sprintf("...(+%d)...", length(cn) - 13), cn[(length(cn) - 2):length(cn)]) else cn
    cat("   cols:     ", paste(show, collapse = ", "), "\n", sep = "")
    key <- grep("variant|type|sample|patient|_id$|^id$|\\bID\\b|cas|name_sub|mode|^p_|pval|p.value|year",
                cn, value = TRUE, ignore.case = TRUE)
    cat("   key cols: ", paste(utils::head(key, 14), collapse = ", "), "\n", sep = "")
    types <- vapply(x[seq_len(min(6, ncol(x)))], function(c) class(c)[1], "")
    cat("   types:    ", paste(sprintf("%s<%s>", names(types), types), collapse = "  "), "\n", sep = "")
    print(utils::head(as.data.frame(x)[, seq_len(min(5, ncol(x))), drop = FALSE], 3))
  } else if (is.list(x)) {
    cat(sprintf("\n▸ %-26s <%s>  length %d; names: %s\n", nm, cl, length(x),
                paste(utils::head(names(x), 12), collapse = ", ")))
  } else {
    cat(sprintf("\n▸ %-26s <%s>  length %d; head: %s\n", nm, cl, length(x),
                paste(utils::head(x, 10), collapse = ", ")))
  }
  invisible()
}

cat("\n###### A. PER-SAMPLE DATA (rows = samples; need sample-ID + variant + year join) ######\n")
for (nm in c("clinical_data", "tumor_pathology_raw",
             "tumors_quant_wt", "tumors_quant", "tumors_qual",
             "tumors_quant_sig_i", "tumors_quant_sig")) desc(nm)

cat("\n\n###### B. PER-CHEMICAL RESULTS (rows = chemicals/fragments; p-values) ######\n")
for (nm in c("anova_results_sig", "fisher_results_i", "fisher_results",
             "posthoc_table_pvalues_i", "posthoc_table_pvalues", "summary_table_i")) desc(nm)

cat("\n\n###### C. VALIDATED / MASTER TABLES (the 30 validated + metadata) ######\n")
for (nm in c("MT_final_i", "MT_final", "MTi",
             "validated_compounds", "validated_variant", "validated_iarc", "validated_list",
             "MT_final_cas_list", "MT_final_namesub_list")) desc(nm)

cat("\n\n###### D. TUMOR vs CADAVER (IARC absolute comparison: o-toluidine / 4-ABP) ######\n")
for (nm in c("full_joiner", "IARC_ttests", "carc_by_variant", "carc_summary",
             "ppm_full_table", "ppm_ppb_inclusive")) desc(nm)

cat("\n\n###### E. CHEMICAL METADATA / NAME MAPS ######\n")
for (nm in c("feature_metadata", "short_name", "MT")) desc(nm)

cat("\n\n###### F. What sample IDs look like (to plan the year join) ######\n")
for (nm in c("clinical_data", "tumor_pathology_raw", "tumors_quant_wt")) {
  if (exists(nm, envir = .GlobalEnv)) {
    x <- get(nm, envir = .GlobalEnv)
    idcol <- intersect(c("Patient_ID", "ID", "sample_ID", "sample", "patient_ID"), names(x))
    if (length(idcol)) cat(sprintf("  %-22s %s: %s\n", nm, idcol[1],
                                    paste(utils::head(unique(x[[idcol[1]]]), 8), collapse = ", ")))
    else cat(sprintf("  %-22s (no obvious sample-ID column; cols: %s)\n", nm,
                     paste(utils::head(names(x), 6), collapse = ", ")))
  }
}

cat("\n\n###### G. RECOVER sample IDs for per-sample tables + attach YEAR (VERIFIED) ######\n")
#! tumors_quant_wt/qual/sig kept only `variant`; their row order was set by the column
#! order of tumor_column when 00c pivoted it. Recover that, join tumor_pathology_raw for
#! year/age/sex, and VERIFY by matching the recovered variant to tumors_quant_wt$variant.
cat("-- ID-bearing candidate objects --\n")
for (nm in c("tumor_column", "tumor", "tumor_seq", "conc_raw")) desc(nm)

if (exists("tumor_column") && exists("tumors_quant_wt") && exists("tumor_seq") && exists("tumor_pathology_raw")) {
  ids <- setdiff(names(tumor_column), "name_sub_lib_id")     # row order of the per-sample tables
  cat(sprintf("\ntumor_column sample cols: n=%d (need 60) | first: %s | last: %s\n",
              length(ids), paste(utils::head(ids, 6), collapse = ", "),
              paste(utils::tail(ids, 3), collapse = ", ")))

  # VERIFY via tumor_seq (already in the quant tables' Follicular/FV-PTC/Papillary scheme).
  # Interleaved order => a 60/60 match is a strong fingerprint of correct row alignment.
  variant_seq <- as.character(tumor_seq$variant[match(ids, tumor_seq$ID)])
  q_ok <- length(ids) == nrow(tumors_quant_wt) &&
          identical(variant_seq, as.character(tumors_quant_wt$variant))
  ql_ok <- !exists("tumors_qual") ||
           identical(variant_seq, as.character(tumors_qual$variant))
  n_match <- sum(variant_seq == as.character(tumors_quant_wt$variant), na.rm = TRUE)
  cat(sprintf("VERIFY (tumor_seq vs tumors_quant_wt, interleaved): %d / %d -> %s | tumors_qual aligns: %s\n",
              n_match, nrow(tumors_quant_wt),
              if (q_ok) "PASS" else "FAIL", if (ql_ok) "yes" else "NO"))

  pr <- tumor_pathology_raw[match(ids, tumor_pathology_raw$Patient_ID), ]
  sample_meta <- data.frame(
    row       = seq_along(ids),
    sample_ID = ids,
    variant   = as.character(tumors_quant_wt$variant),                  # quant scheme
    Variant   = as.character(pr$Variant),                              # FTC/FV-PTC/PTC
    year      = suppressWarnings(as.numeric(pr$year)),
    Age       = suppressWarnings(as.numeric(pr$Age)),
    Sex       = ifelse(pr$Sex == 1, "Female", ifelse(pr$Sex == 0, "Male", NA_character_)),
    stringsAsFactors = FALSE
  )
  cat(sprintf("Patient_ID matched: %d/%d | year non-missing: %d/%d\n",
              sum(!is.na(match(ids, tumor_pathology_raw$Patient_ID))), length(ids),
              sum(!is.na(sample_meta$year)), length(ids)))
  cat("Recovered per-row metadata (head):\n"); print(utils::head(sample_meta, 6))
  if (q_ok && ql_ok) {
    assign("sample_meta", sample_meta, envir = .GlobalEnv)
    cat("\n=> PASS: sample_meta stashed. Row i aligns to row i of tumors_quant_wt / tumors_qual.\n",
        "   Year + Age + Sex are now attachable to every per-sample chemical value. UNBLOCKED.\n", sep = "")
  } else {
    cat("\n=> alignment not clean; do not proceed on a positional join.\n")
  }
} else {
  cat("Missing one of: tumor_column / tumor_seq / tumors_quant_wt / tumor_pathology_raw\n")
}
