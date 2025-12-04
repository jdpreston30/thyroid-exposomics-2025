#* Manual Spectral Validation
#+ Data
mz_reference_table
cadaver_top5_iarc
tumor_top5_iarc

renv::snapshot()

base_dir <- config$paths$tumor_raw_dir  # or use cadaver_raw_dir

# File to process
raw_file <- "BP2-1_1.raw"

# Full path to raw file
raw_path <- file.path(base_dir, raw_file)

# Target m/z values for chromatogram extraction
target_mzs <- c(282.98, 139.05, 204.02, 254.98)

# m/z tolerance (in Da)
mz_tolerance <- 0.5

# RT range (in minutes) - set to NULL for full range
rt_range <- c(10.1, 10.9)  # e.g., c(5, 15) for 5-15 minutes

# Output directory for converted files
output_dir <- file.path(dirname(raw_path), "converted_mzML")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# =============================================================================
# Convert .raw to mzML using ThermoRawFileParser
# =============================================================================

cat("Converting .raw file to mzML format...\n")

# Output mzML file path
mzml_file <- file.path(output_dir, gsub("\\.raw$", ".mzML", raw_file, ignore.case = TRUE))

# Paths
mono_path <- "/Library/Frameworks/Mono.framework/Versions/Current/Commands/mono"
trfp_exe <- path.expand("~/bin/ThermoRawFileParser/ThermoRawFileParser.exe")

# Check if files exist
if (!file.exists(mono_path)) {
  stop(sprintf("Mono not found at: %s\nPlease ensure Mono is installed.", mono_path))
}
if (!file.exists(trfp_exe)) {
  stop(sprintf("ThermoRawFileParser not found at: %s\nPlease ensure it's installed correctly.", trfp_exe))
}

# Convert .raw to mzML
conversion_cmd <- sprintf('"%s" "%s" -i="%s" -o="%s" -f=1', 
                          mono_path, trfp_exe, raw_path, output_dir)
cat(sprintf("Running conversion...\n"))
system(conversion_cmd)

# Check if conversion was successful
if (!file.exists(mzml_file)) {
  stop(sprintf("Conversion failed. Expected output file not found: %s", mzml_file))
}

cat(sprintf("Conversion complete: %s\n", basename(mzml_file)))

# =============================================================================
# Read mzML file
# =============================================================================

cat("\nReading mzML file...\n")
ms_data <- mzR::openMSfile(mzml_file)

# Get file header information
header_info <- mzR::header(ms_data)
cat(sprintf("Total scans: %d\n", nrow(header_info)))
cat(sprintf("RT range: %.2f - %.2f minutes\n", 
            min(header_info$retentionTime) / 60, 
            max(header_info$retentionTime) / 60))

cat("\nFile loaded successfully!\n")

# =============================================================================
# Extract Chromatograms for Target m/z Values
# =============================================================================

cat("\nExtracting chromatograms for target m/z values...\n")

# Function to extract XIC (extracted ion chromatogram)
extract_xic <- function(ms_data, header_info, target_mz, mz_tol, rt_range = NULL) {
  
  # Filter MS1 scans only
  ms1_scans <- header_info[header_info$msLevel == 1, ]
  
  # Apply RT filter if specified
  if (!is.null(rt_range)) {
    ms1_scans <- ms1_scans[ms1_scans$retentionTime >= rt_range[1] * 60 & 
                           ms1_scans$retentionTime <= rt_range[2] * 60, ]
  }
  
  # Extract intensity for target m/z across all scans
  xic_data <- lapply(ms1_scans$seqNum, function(scan_num) {
    # Get spectrum
    spectrum <- mzR::peaks(ms_data, scan_num)
    
    # Find peaks within m/z tolerance
    mz_match <- abs(spectrum[, 1] - target_mz) <= mz_tol
    
    # Sum intensities if multiple peaks match
    if (sum(mz_match) > 0) {
      intensity <- sum(spectrum[mz_match, 2])
    } else {
      intensity <- 0
    }
    
    data.frame(
      scan = scan_num,
      rt = ms1_scans$retentionTime[ms1_scans$seqNum == scan_num] / 60,
      intensity = intensity,
      mz = target_mz
    )
  })
  
  # Combine all scans
  do.call(rbind, xic_data)
}

# Extract XICs for all target m/z values
all_xics <- lapply(target_mzs, function(mz) {
  cat(sprintf("  Extracting m/z %.2f...\n", mz))
  extract_xic(ms_data, header_info, mz, mz_tolerance, rt_range)
})

# Combine into single data frame
chromatogram_data <- do.call(rbind, all_xics)
chromatogram_data$mz_label <- sprintf("m/z %.2f", chromatogram_data$mz)

cat(sprintf("\nExtracted %d data points across %d m/z values\n", 
            nrow(chromatogram_data), length(target_mzs)))

# =============================================================================
# Plot Chromatograms
# =============================================================================

cat("\nGenerating chromatogram plot...\n")

# Create plot
p <- ggplot(chromatogram_data, aes(x = rt, y = intensity, color = mz_label)) +
  geom_line(linewidth = 0.8) +
  labs(
    title = sprintf("Extracted Ion Chromatograms - %s", raw_file),
    x = "Retention Time (minutes)",
    y = "Intensity",
    color = "m/z"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "right",
    plot.title = element_text(hjust = 0.5, face = "bold")
  ) +
  scale_color_brewer(palette = "Set1")

# Display plot
print(p)

# Save plot
plot_file <- file.path(output_dir, gsub("\\.raw$", "_chromatogram.png", raw_file, ignore.case = TRUE))
ggsave(plot_file, p, width = 10, height = 6, dpi = 300)
cat(sprintf("Plot saved: %s\n", plot_file))

# Close mzML file
mzR::close(ms_data)

cat("\nAnalysis complete!\n")