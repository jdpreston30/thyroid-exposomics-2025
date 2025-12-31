# Test script to verify ST3 generation works
source("R/Utilities/Helpers/load_dynamic_config.R")
config <- load_dynamic_config(computer = "auto", config_path = "All_Run/config_dynamic.yaml")

# Load just the utilities we need
source("R/Utilities/Helpers/update.R")
u()

# Check if build_ST3 function exists and what columns it expects
cat("\n=== Testing build_ST3 function ===\n")
cat("Function exists:", exists("build_ST3"), "\n")

# Check the actual function body for column references
build_ST3_body <- deparse(build_ST3)
has_old_column <- any(grepl("Mean Non-Cancer Thyroid", build_ST3_body, fixed = TRUE))
cat("Contains old 'Mean Non-Cancer Thyroid' column:", has_old_column, "\n")

if (has_old_column) {
  cat("\n❌ ERROR: build_ST3 still contains reference to deleted column!\n")
  cat("Lines containing the old column:\n")
  print(grep("Mean Non-Cancer", build_ST3_body, value = TRUE))
} else {
  cat("\n✅ build_ST3 function is clean - no references to deleted column\n")
}

cat("\nColumns expected by build_ST3:\n")
print(grep("columns = ", build_ST3_body, value = TRUE))
