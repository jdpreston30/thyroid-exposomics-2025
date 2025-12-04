#* 0c: Load and Clean Clinical Data
#+ 0c.1: Import clinical demographics data; clean
demographics <- read_excel(config$paths$clinical_metadata, sheet = "Manifest") |>
  select(ID, variant, sex, age_collection, year_collection) |>
  mutate(year_bin = cut(as.numeric(year_collection),
    breaks = seq(2006, 2022, length.out = 5),
    labels = c("2006-2009", "2010-2013", "2014-2017", "2018-2021"),
    include.lowest = TRUE
  ))
