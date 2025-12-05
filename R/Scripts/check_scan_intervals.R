# Check scan intervals in mzML file
library(mzR)
library(yaml)

# Load config
config <- yaml::read_yaml("All_Run/config_dynamic.yaml")

# Detect which computer
if (Sys.info()["user"] == "jdp2019") {
  computer <- config$computers$laptop
} else {
  computer <- config$computers$desktop
}

# Build paths
base_data_path <- gsub("\\{user_home\\}", computer$user_home, config$paths$base_data_path)
base_data_path <- gsub("\\{onedrive_path\\}", computer$onedrive_path, base_data_path)
tumor_raw_dir <- gsub("\\{base_data_path\\}", base_data_path, config$paths$tumor_raw_dir)

mzml_dir <- file.path(tumor_raw_dir, "mzML_validation")
mzml_path <- file.path(mzml_dir, "BL_12082022_037.mzML")

cat("Looking for file:", mzml_path, "\n")

if (!file.exists(mzml_path)) {
  stop("File not found!")
}

# Open mzML file
ms_data <- mzR::openMSfile(mzml_path)
header_info <- mzR::header(ms_data)

# Get MS1 scans only
ms1_scans <- header_info[header_info$msLevel == 1, ]

# Calculate time intervals between consecutive scans
time_diffs <- diff(ms1_scans$retentionTime)

# Convert to milliseconds and minutes
time_diffs_ms <- time_diffs * 1000
time_diffs_min <- time_diffs / 60

cat("\n=== Scan Interval Statistics (MS1 scans only) ===\n")
cat(sprintf("  Mean: %.3f seconds (%.1f ms, %.5f min)\n", mean(time_diffs), mean(time_diffs_ms), mean(time_diffs_min)))
cat(sprintf("  Median: %.3f seconds (%.1f ms, %.5f min)\n", median(time_diffs), median(time_diffs_ms), median(time_diffs_min)))
cat(sprintf("  Min: %.3f seconds (%.1f ms)\n", min(time_diffs), min(time_diffs_ms)))
cat(sprintf("  Max: %.3f seconds (%.1f ms)\n", max(time_diffs), max(time_diffs_ms)))
cat(sprintf("  Total MS1 scans: %d\n", nrow(ms1_scans)))

# Show first few intervals
cat("\nFirst 10 scan intervals (seconds):\n")
print(head(time_diffs, 10))

# Check around the RT range from your test
rt_window <- ms1_scans[ms1_scans$retentionTime >= 14.97*60 & ms1_scans$retentionTime <= 15.02*60, ]
cat(sprintf("\n=== Scans in RT window 14.97-15.02 minutes ===\n"))
cat(sprintf("  Number of scans: %d\n", nrow(rt_window)))
if (nrow(rt_window) > 1) {
  window_diffs <- diff(rt_window$retentionTime)
  cat(sprintf("  Mean interval: %.3f seconds (%.1f ms)\n", mean(window_diffs), mean(window_diffs)*1000))
}

mzR::close(ms_data)
