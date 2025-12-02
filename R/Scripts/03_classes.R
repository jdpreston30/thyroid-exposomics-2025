#* 3: Use Classes and Detection Distribution
#+ 3.1: Summarize Usage Classes
all <- feature_metadata |>
  select(Graph_Class) |>
  count(Graph_Class) |>
  arrange(desc(n))
#+ 3.2: IARC and EDC Class Summaries
#- 3.2.1: IARC Group Summary
IARC <- feature_metadata |>
  select(IARC_Group) |>
  count(IARC_Group) |>
  mutate(IARC_Group = ifelse(is.na(IARC_Group), "Not Classified", paste0("Group ", IARC_Group))) |>
  mutate(IARC_Group = factor(IARC_Group, levels = c("Group 1", "Group 2A", "Group 2B", "Group 3", "Not Classified"))) |>
  arrange(IARC_Group)
#- 3.2.2: EDC Summary
EDC <- feature_metadata |>
  select(Potential_EDC) |>
  count(Potential_EDC) |>
  arrange(desc(n)) |>
  mutate(Potential_EDC = ifelse(is.na(Potential_EDC) | Potential_EDC != "Y", "Non-EDC", "Potential EDC")) |>
  mutate(Potential_EDC = factor(Potential_EDC, levels = c("Non-EDC", "Potential EDC")))
#+ 3.3: Chemical Classes and Superclasses
#- 3.3.1: Superclass Summary
superclasses <- feature_metadata |>
  select(Superclass) |>
  mutate(Superclass = ifelse(is.na(Superclass), "Other", Superclass)) |>
  count(Superclass) |>
  mutate(Superclass = str_to_sentence(Superclass)) |>
  arrange(Superclass == "Other", desc(n), Superclass) |>
  rename(Graph_Class = Superclass)
#- 3.3.2: Class Summary
classes <- feature_metadata |>
  select(Class) |>
  mutate(Class = ifelse(is.na(Class), "Other", Class)) |>
  count(Class) |>
  mutate(Class = ifelse(n <= 9, "Other", Class)) |>
  group_by(Class) |>
  summarise(n = sum(n), .groups = "drop") |>
  mutate(Class = str_to_sentence(Class)) |>
  arrange(Class == "Other", desc(n), Class) |>
  rename(Graph_Class = Class)
#+ 3.4: Create Figures
#- 3.4.1: Use Class Distribution Plot
p2A <- plot_class_distribution(all, x_max = 100)
#- 3.4.2: IARC Donut Chart
p2B <- plot_iarc_donut(IARC)
#- 3.4.3: EDC Donut Chart
p2C <- plot_edc_donut(EDC)
#- 3.4.4: Class and Superclass Distribution Plots
s1A <- plot_class_distribution(superclasses, x_max = 250, sup = TRUE)
s1B <- plot_class_distribution(classes, x_max = 200, sup = TRUE)
