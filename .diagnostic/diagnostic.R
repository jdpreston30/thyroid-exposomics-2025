#* 1: Classification-path audit — is the flowchart's 2nd-to-last diamond needed?
#! Question: can we collapse "q_iarc3 (IARC3 & no GHS -> Uncertain)" + "q_none (no IARC
#! & no GHS -> Unclassified, else catch-all Uncertain)" into ONE diamond? Only safe if
#! NO chemical reaches the CATCH-ALL (iarc==3 or iarc==NA WHILE a GHS carcinogen statement
#! is present but non-qualifying, e.g. H350 listed at 0%). This traces the real path each
#! classified chemical took. Requires MTi (from script 05).
stopifnot(exists("MTi"))

#+ 1.1: Path-tracing mirror of classify_carcinogenicity() (returns the terminal + WHY)
classify_path <- function(ghs, iarc) {
  if (is.na(ghs)  || ghs  == "NA" || ghs  == "no ghs carcinogen statement") ghs  <- NA
  if (is.na(iarc) || iarc == "NA" || iarc == "Not classified")               iarc <- NA
  h350  <- as.numeric(stringr::str_extract(tolower(ghs), "(?<=h350 \\()\\d+\\.?\\d*"))
  h350i <- as.numeric(stringr::str_extract(tolower(ghs), "(?<=h350i \\()\\d+\\.?\\d*"))
  h351  <- as.numeric(stringr::str_extract(tolower(ghs), "(?<=h351 \\()\\d+\\.?\\d*"))
  if (!is.na(iarc) && iarc == "1")                    return("Known|iarc1")
  if (!is.na(iarc) && iarc == "2A")                   return("Likely|iarc2A")
  if (!is.na(h350)  && h350  >= 50)                   return("Likely|h350>=50")
  if (!is.na(h350i) && h350i >= 50)                   return("Likely|h350i>=50")
  if (!is.na(iarc) && iarc == "2B")                   return("Possible|iarc2B")
  if (!is.na(h350)  && h350  > 0 && h350  < 50)       return("Possible|h350 0-50")
  if (!is.na(h350i) && h350i > 0 && h350i < 50)       return("Possible|h350i 0-50")
  if (!is.na(h351)  && h351  > 0)                      return("Possible|h351>0")
  if (!is.na(iarc) && iarc == "3" && is.na(h350) && is.na(h350i) && is.na(h351))
    return("Uncertain|q_iarc3 (IARC3, no GHS)")
  if (is.na(iarc) && is.na(ghs))                      return("Unclassified|q_none (nothing)")
  return("Uncertain|CATCH-ALL (needs 2nd diamond)")
}

paths <- mapply(classify_path, MTi$GHS_var_diff_only, MTi$IARC_Group, USE.NAMES = FALSE)

#+ 1.2: Path frequencies
cat(sprintf("\n=== Classification paths across %d classified chemicals ===\n", length(paths)))
print(sort(table(paths), decreasing = TRUE))

#+ 1.3: Sanity — does the mirror reproduce the pipeline's Carcinogenicity?
term_map <- c(Known = "Known Carcinogen", Likely = "Likely Carcinogen",
              Possible = "Possible Carcinogen", Uncertain = "Uncertain Risk",
              Unclassified = "Unclassified")
my_term <- unname(term_map[sub("\\|.*", "", paths)])
cat(sprintf("Sanity: mirror terminals match pipeline Carcinogenicity in %d/%d rows%s\n",
            sum(my_term == MTi$Carcinogenicity, na.rm = TRUE), nrow(MTi),
            if (all(my_term == MTi$Carcinogenicity, na.rm = TRUE)) " (PASS)" else " (MISMATCH!)"))

#+ 1.4: The verdict — do any chemicals need the 2nd diamond?
catchall <- grepl("CATCH-ALL", paths)
uncls    <- grepl("q_none",    paths)
cat(sprintf("\nCATCH-ALL cases (simplification would be WRONG): %d\n", sum(catchall)))
cat(sprintf("q_none 'Unclassified' cases:                     %d\n", sum(uncls)))
if (any(catchall)) {
  cat("\n=> KEEP the two-diamond version. These chemicals require it:\n")
  print(data.frame(short_name = MTi$short_name[catchall],
                   IARC = MTi$IARC_Group[catchall],
                   GHS  = MTi$GHS_var_diff_only[catchall],
                   result = MTi$Carcinogenicity[catchall]), row.names = FALSE)
} else {
  cat("\n=> No catch-all cases in this data. The 1-diamond simplification would give\n",
      "   IDENTICAL results HERE (still not equivalent to the general function, but\n",
      "   accurate for this dataset). Your instinct holds for the actual data.\n", sep = "")
}
