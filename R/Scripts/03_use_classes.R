#* 3: Use Classes and Detection Distribution
#+ 3.1: Use Classes (Figure 2A)
feature_metadata <- read_excel(config$paths$chemical_metadata, sheet = "feature_metadata")
all <- feature_metadata %>%
  select(Graph_Class) %>%
  count(Graph_Class) %>%
  arrange(desc(n))
write.csv(all, file.path(config$paths$output, "Tables/all.csv"))

#+ 3.2: Detection Distribution Plot (Figure 2A)
