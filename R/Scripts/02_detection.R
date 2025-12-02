#* 2: Detection Analysis
#+ 2.1: Group by CAS and get total unique which were detected
detection <- tumor_raw %>%
  mutate(cas = as.character(cas)) %>%
  mutate(across(where(is.numeric), ~ ifelse(is.na(.), 0, 1))) %>%
  select(cas, "F1":"F20") %>%
  group_by(cas) %>%
  summarise(across(where(is.numeric), max)) %>%
  pivot_longer(-cas, names_to = "ID", values_to = "Value") %>%
  pivot_wider(names_from = cas, values_from = Value) %>%
  left_join(tumor_seq, by = "ID") %>%
  select(ID, variant, everything(), -ID) %>%
  arrange(variant)
#+ 2.2: Remove endogenous features from detection analysis
#- 2.2.1: Import the list which is purely endogenous which was compiled externally from above
endog_cas <- read_excel(config$paths$chemical_metadata, sheet = "endogenous_excluded_features") %>%
  pull(cas)
#- 2.2.2: Remove endogenous features and calculate total detected per sample
detection_no_endog <- detection %>%
  select(-any_of(endog_cas)) %>%
  mutate(total_detected = rowSums(across(-variant))) %>%
  select(variant, total_detected) %>%
  arrange(variant, total_detected)
#- Now, make a frequency distribution with bins for visualization
freq_dist_bins <- detection_no_endog %>%
  mutate(Bin = cut(total_detected, breaks = seq(260, 330, by = 10), right = FALSE)) %>%
  count(Bin, variant) %>%
  pivot_wider(names_from = variant, values_from = n, values_fill = 0) %>%
  arrange(Bin) %>%
  select(Bin, Follicular, everything())
#+ 2.3: Statistical testing
#- 2.3.1: Kruskal-Wallis test comparing medians across variants
kw_result <- kruskal.test(total_detected ~ variant, data = detection_no_endog)$p.value
#+ 2.4: Create Figures
p1A <- plot_detection_distribution(freq_dist_bins)
p1B <- plot_detection_scatter(detection_no_endog) 
