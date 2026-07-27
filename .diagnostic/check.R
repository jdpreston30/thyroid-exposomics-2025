#* 1: ST3 Supplementary-Table Verification
#! Run AFTER sourcing R/Scripts/16_supplementary_tables.R (so ST3_tern, demo_S3_data,
#! build_ST3, etc. exist). Prints the exact data structures the ST3 build depends on
#! so they can be pasted back for review. Read-only; writes nothing except stdout.
cat("\n================ ST3 DIAGNOSTIC ================\n")

#+ 1.1: Presence of required objects
.need <- c("ST3_tern", "demo_S3_data", "cadaver_metadata", "tumor_pathology_raw", "build_ST3")
cat("\n--- 1.1 object presence ---\n")
for (o in .need) cat(sprintf("  %-22s %s\n", o, if (exists(o)) "OK" else "*** MISSING ***"))
if (!all(vapply(.need, exists, logical(1)))) {
  cat("\nMissing objects — source R/Scripts/16_supplementary_tables.R (and u()) first.\n")
}

#+ 1.2: ST3_tern — exact shape the emitter keys off
if (exists("ST3_tern")) {
  cat("\n--- 1.2 ST3_tern column names (dput) ---\n"); dput(names(ST3_tern))
  cat("\n--- 1.2 ST3_tern .indent column (dput) ---\n")
  if (".indent" %in% names(ST3_tern)) dput(ST3_tern$.indent) else cat("  NO .indent column (indent_info_column=TRUE not applied?)\n")
  cat("\n--- 1.2 ST3_tern full content (as data.frame) ---\n")
  print(as.data.frame(ST3_tern, check.names = FALSE), right = FALSE)
  cat("\n--- 1.2 ST3_tern column classes ---\n"); print(sapply(ST3_tern, function(x) base::class(x)[1]))  # base:: dodges leaked `class` var from script 16 ST2 loop
}

#+ 1.3: build_ST3 output on the REAL object (what gets written to ST3.tex)
if (exists("ST3_tern") && exists("build_ST3")) {
  cat("\n--- 1.3 build_ST3(ST3_tern) LaTeX ---\n")
  cat(tryCatch(build_ST3(ST3_tern), error = function(e) paste("ERROR:", conditionMessage(e))), "\n")
}

#+ 1.4: Cohort frame — group Ns and level coding
if (exists("demo_S3_data")) {
  cat("\n--- 1.4 demo_S3_data: n per cohort ---\n"); print(table(demo_S3_data$Cohort, useNA = "ifany"))
  cat("\n--- 1.4 Sex x Cohort ---\n"); print(table(demo_S3_data$Sex, demo_S3_data$Cohort, useNA = "ifany"))
  cat("\n--- 1.4 Sample Collection Timing x Cohort ---\n"); print(table(demo_S3_data$`Sample Collection Timing`, demo_S3_data$Cohort, useNA = "ifany"))
  cat("\n--- 1.4 Age summary by cohort (mean/sd/min/max) ---\n")
  print(aggregate(Age ~ Cohort, demo_S3_data, function(x) c(mean = round(mean(x), 1), sd = round(sd(x), 1), min = min(x), max = max(x))))
}

#+ 1.5: Categorical test ternG used for Sex (Fisher vs chi-square) — verify caption wording
if (exists("demo_S3_data")) {
  cat("\n--- 1.5 sex test cross-check (caption says Fisher) ---\n")
  .tab <- table(demo_S3_data$Sex, demo_S3_data$Cohort)
  cat("  min expected cell count:", round(min(suppressWarnings(chisq.test(.tab)$expected)), 2), "(<5 => Fisher is appropriate)\n")
  cat("  fisher p =", signif(fisher.test(.tab)$p.value, 3),
      "| chisq p =", signif(suppressWarnings(chisq.test(.tab)$p.value), 3), "\n")
}

#+ 1.6: Source-data sanity — DTC filter recovers 60; cadaver = 8
if (exists("tumor_pathology_raw")) {
  cat("\n--- 1.6 tumor_pathology_raw: Sex coding + year range + DTC filter count ---\n")
  cat("  unique Sex values:", paste(sort(unique(as.character(tumor_pathology_raw$Sex))), collapse = ", "), "\n")
  cat("  year range:", paste(range(suppressWarnings(as.numeric(tumor_pathology_raw$year)), na.rm = TRUE), collapse = "-"), "\n")
  .n_dtc <- sum(grepl("^(F|P|FVPTC)\\d+$", tumor_pathology_raw$Patient_ID))
  cat("  Patient_IDs matching ^(F|P|FVPTC)\\d+$ :", .n_dtc, "(expect 60)\n")
}
if (exists("cadaver_metadata")) {
  cat("\n--- 1.6 cadaver_metadata (expect n=8, all 2019) ---\n"); print(as.data.frame(cadaver_metadata))
}

#+ 1.7: Confirm Table 1 (table1) unaffected by the new ternG call
cat("\n--- 1.7 Table 1 object still present & intact ---\n")
if (exists("table1")) { cat("  table1 dims:", paste(dim(as.data.frame(table1)), collapse = " x "), "\n") } else cat("  table1 not in env (fine if 15 not sourced this session)\n")
cat("\n================ END DIAGNOSTIC ================\n")
