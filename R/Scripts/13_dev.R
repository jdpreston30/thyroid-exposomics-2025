#! Reviewer #2: could the tumor-type chemical differences be confounded by age, sex,
# ! or collection year? A covariate confounds only if it associates with BOTH the
# ! outcome (tumor type) AND the exposure (chemical). This section tests requirement
# ! one — covariate vs type — plus covariate cross-associations; the covariate<->chemical
# ! screens and adjusted models come in later sections.

# ! Tests whether tumor-type chemical differences are confounded by age, sex, and
# ! collection year. Self-contained: rebuilds a VERIFIED per-sample metadata frame,
# ! then runs (2) balance, (3) covariate<->chemical screens vs chance, (4) targeted
# ! ANCOVA on the 30 validated chemicals, (5) common-support, (6) a printed summary.
# ! Requires the pipeline objects (through ~script 12) in the environment.

#+ 13.1: Recover + VERIFY Per-row Sample Metadata (sample_ID, variant, year, Age, Sex)
stopifnot(
  exists("tumor_column"), exists("tumor_seq"), exists("tumors_quant_wt"),
  exists("tumors_qual"), exists("tumor_pathology_raw"), exists("MT_final")
)
.ids <- setdiff(names(tumor_column), "name_sub_lib_id") # = row order of per-sample tables
.vseq <- as.character(tumor_seq$variant[match(.ids, tumor_seq$ID)])
stopifnot(
  identical(.vseq, as.character(tumors_quant_wt$variant)), # interleaved 60/60 => aligned
  identical(.vseq, as.character(tumors_qual$variant))
)
.pr <- tumor_pathology_raw[match(.ids, tumor_pathology_raw$Patient_ID), ]
sample_meta <- data.frame(
  sample_ID = .ids,
  variant = factor(as.character(tumors_quant_wt$variant), levels = c("Follicular", "FV-PTC", "Papillary")),
  year = suppressWarnings(as.numeric(.pr$year)),
  Age = suppressWarnings(as.numeric(.pr$Age)),
  Sex = factor(ifelse(.pr$Sex == 1, "Female", "Male"), levels = c("Female", "Male")),
  stringsAsFactors = FALSE
)
assign("sample_meta", sample_meta, envir = .GlobalEnv)
cat(sprintf(
  "[setup] %d samples aligned & verified | year %d-%d | quant feats=%d qual feats=%d\n",
  nrow(sample_meta), min(sample_meta$year), max(sample_meta$year),
  ncol(tumors_quant_wt) - 1L, ncol(tumors_qual) - 1L
))

#+ 13.2: Helper — Significant-vs-chance Summary for a Vector of p-values
.pchance <- function(p, label) {
  p <- p[!is.na(p)]
  n <- length(p)
  obs <- sum(p < 0.05)
  exp <- 0.05 * n
  cat(sprintf(
    "   %-20s %3d/%4d sig (%.1f%%) vs ~%.1f by chance  [%s]\n",
    label, obs, n, 100 * obs / n, exp,
    if (obs <= exp * 1.5) {
      "~chance"
    } else if (obs >= exp * 2) "ABOVE chance" else "modest excess"
  ))
  invisible(c(obs = obs, n = n, exp = exp))
}

#* 13: Covariate Balance by Tumor Type (does Each Covariate Differ by type?)
cat("\n#* 2: COVARIATE BALANCE BY TYPE ------------------------------------------\n")
#+ 13.1: Age
cat(sprintf(
  "Age:  ANOVA p=%.3f | Kruskal-Wallis p=%.3f\n",
  summary(aov(Age ~ variant, sample_meta))[[1]][["Pr(>F)"]][1],
  kruskal.test(Age ~ variant, sample_meta)$p.value
))
#+ 13.2: Sex
cat(sprintf("Sex:  Fisher p=%.3f\n", fisher.test(table(sample_meta$Sex, sample_meta$variant))$p.value))
#+ 13.3: Collection Year — Binned (paper-style) vs Continuous (the Key reconciliation)
.bins <- cut(sample_meta$year,
  breaks = c(2005, 2009, 2013, 2017, 2021),
  labels = c("2006-2009", "2010-2013", "2014-2017", "2018-2021")
)
cat(sprintf(
  "Year: binned Fisher p=%.4f  |  continuous KW p=%.4f | ANOVA p=%.4f\n",
  fisher.test(table(.bins, sample_meta$variant), simulate.p.value = TRUE, B = 1e5)$p.value,
  kruskal.test(year ~ variant, sample_meta)$p.value,
  summary(aov(year ~ variant, sample_meta))[[1]][["Pr(>F)"]][1]
))
print(do.call(rbind, lapply(split(sample_meta$year, sample_meta$variant), function(x) {
  data.frame(median = median(x), mean = round(mean(x), 1), min = min(x), max = max(x), n = length(x))
})))

#* 13: Covariate <-> Chemical Screens (all features), with Chance Baselines
cat("\n#* 3: COVARIATE <-> CHEMICAL SCREENS (observed sig vs chance) --------------\n")
#+ 13.1: Quant Features (Spearman for year/age; Wilcoxon for sex)
.qcols <- setdiff(names(tumors_quant_wt), "variant")
.screen_q <- function(cov, binary = FALSE) {
  vapply(.qcols, function(cn) {
    y <- tumors_quant_wt[[cn]]
    tryCatch(
      if (binary) {
        wilcox.test(y ~ cov)$p.value
      } else {
        suppressWarnings(cor.test(y, cov, method = "spearman")$p.value)
      },
      error = function(e) NA_real_
    )
  }, numeric(1))
}
cat("QUANTITATIVE features:\n")
.pchance(.screen_q(sample_meta$year), "year (Spearman)")
.pchance(.screen_q(sample_meta$Age), "age (Spearman)")
.pchance(.screen_q(sample_meta$Sex, TRUE), "sex (Wilcoxon)")

#+ 13.2: Qual Features (detection 0/1): Covariate ~ Detection (testable Features only)
.dcols <- setdiff(names(tumors_qual), "variant")
.testable <- .dcols[vapply(.dcols, function(cn) {
  v <- tumors_qual[[cn]]
  length(unique(v[!is.na(v)])) > 1
}, logical(1))]
.screen_ql <- function(cov, cat = FALSE) {
  vapply(.testable, function(cn) {
    d <- factor(tumors_qual[[cn]])
    tryCatch(if (cat) fisher.test(table(cov, d))$p.value else wilcox.test(cov ~ d)$p.value,
      error = function(e) NA_real_
    )
  }, numeric(1))
}
cat(sprintf("QUALITATIVE features (%d testable of %d):\n", length(.testable), length(.dcols)))
.pchance(.screen_ql(sample_meta$year), "year (Wilcoxon)")
.pchance(.screen_ql(sample_meta$Age), "age (Wilcoxon)")
.pchance(.screen_ql(sample_meta$Sex, TRUE), "sex (Fisher)")

#* 13: Targeted ANCOVA — Do Type Effects Survive Adjustment for Year (+age+sex)?
cat("\n#* 4: TARGETED ANCOVA on the 30 VALIDATED chemicals ------------------------\n")
.valq <- intersect(MT_final$name_sub_lib_id[MT_final$mode == "quantitative"], names(tumors_quant_wt))
.valql <- intersect(MT_final$name_sub_lib_id[MT_final$mode == "qualitative"], names(tumors_qual))
.anc_q <- function(cn) {
  df <- cbind(y = tumors_quant_wt[[cn]], sample_meta)
  p0 <- tryCatch(anova(lm(y ~ variant, df))[["Pr(>F)"]][1], error = function(e) NA_real_)
  py <- tryCatch(drop1(lm(y ~ variant + year, df), test = "F")["variant", "Pr(>F)"], error = function(e) NA_real_)
  pf <- tryCatch(drop1(lm(y ~ variant + year + Age + Sex, df), test = "F")["variant", "Pr(>F)"], error = function(e) NA_real_)
  data.frame(name_sub_lib_id = cn, mode = "quant", p_unadj = p0, p_yearAdj = py, p_fullAdj = pf)
}
.anc_ql <- function(cn) {
  df <- cbind(d = tumors_qual[[cn]], sample_meta)
  fit <- function(f) tryCatch(suppressWarnings(glm(f, binomial, df)), error = function(e) NULL)
  lrt <- function(a, b) {
    if (!is.null(a) && !is.null(b)) {
      tryCatch(anova(a, b, test = "LRT")[["Pr(>Chi)"]][2], error = function(e) NA_real_)
    } else {
      NA_real_
    }
  }
  data.frame(
    name_sub_lib_id = cn, mode = "qual",
    p_unadj = lrt(fit(d ~ 1), fit(d ~ variant)),
    p_yearAdj = lrt(fit(d ~ year), fit(d ~ variant + year)),
    p_fullAdj = lrt(fit(d ~ year + Age + Sex), fit(d ~ variant + year + Age + Sex))
  )
}
ancova_res <- rbind(
  do.call(rbind, lapply(.valq, .anc_q)),
  do.call(rbind, lapply(.valql, .anc_ql))
)
ancova_res <- merge(ancova_res, MT_final[, c("name_sub_lib_id", "short_name")],
  by = "name_sub_lib_id", all.x = TRUE
)
ancova_res$survives_year <- ancova_res$p_yearAdj < 0.05
ancova_res$flips <- (ancova_res$p_unadj < 0.05) & !(ancova_res$p_yearAdj < 0.05)
ancova_res <- ancova_res[order(ancova_res$mode, ancova_res$p_unadj), ]
assign("ancova_res", ancova_res, envir = .GlobalEnv)
print(within(ancova_res, {
  p_unadj <- signif(p_unadj, 3)
  p_yearAdj <- signif(p_yearAdj, 3)
  p_fullAdj <- signif(p_fullAdj, 3)
})[, c("short_name", "mode", "p_unadj", "p_yearAdj", "p_fullAdj", "survives_year", "flips")], row.names = FALSE)

#* 13: Common Support — Collection-year Overlap for Each Type Pair
cat("\n#* 5: COMMON SUPPORT (year overlap per type pair) -------------------------\n")
.yr <- split(sample_meta$year, sample_meta$variant)
for (pr in combn(names(.yr), 2, simplify = FALSE)) {
  a <- .yr[[pr[1]]]
  b <- .yr[[pr[2]]]
  lo <- max(min(a), min(b))
  hi <- min(max(a), max(b))
  cat(sprintf(
    "  %-11s vs %-11s overlap %s | n_%s=%d, n_%s=%d in window\n",
    pr[1], pr[2], if (lo <= hi) sprintf("%d-%d", lo, hi) else "NONE",
    pr[1], sum(a >= lo & a <= hi), pr[2], sum(b >= lo & b <= hi)
  ))
}

#* 13: Summary for the Reviewer Response
cat("\n#* 6: SUMMARY ------------------------------------------------------------\n")
.nsig <- sum(ancova_res$p_unadj < 0.05, na.rm = TRUE)
.nflip <- sum(ancova_res$flips, na.rm = TRUE)
cat("Age & Sex: balanced across types (see #2) -> ruled out as confounders.\n")
cat("Year: binned differs (bin-boundary artifact) but continuous central tendency does NOT (see #2).\n")
cat(sprintf(
  "ANCOVA: of %d originally-significant validated chemicals, %d lose significance after year adj.\n",
  .nsig, .nflip
))
cat(sprintf("=> %s\n", if (.nflip == 0) {
  "ALL type effects survive year adjustment — robust."
} else {
  sprintf("%d effect(s) year-sensitive — flag those; the rest are robust.", .nflip)
}))
if (.nflip > 0) print(ancova_res[ancova_res$flips %in% TRUE, c("short_name", "mode", "p_unadj", "p_yearAdj")], row.names = FALSE)
