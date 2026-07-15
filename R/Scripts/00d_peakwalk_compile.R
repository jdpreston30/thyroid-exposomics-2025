#* 0d: PeakWalk Compilation
#+ 0d.1: Import Data
#- 0d.1.1: Tumors intensity tables from PeakWalk
bp1_tumor <- read_csv(file.path(config$paths$peakwalk_tumor_dir, "BP1.GC2/feature.sample.i.csv"))
bp2_tumor <- read_csv(file.path(config$paths$peakwalk_tumor_dir, "BP2.GC2/feature.sample.i.csv"))
bp3_tumor <- read_csv(file.path(config$paths$peakwalk_tumor_dir, "BP3.GC2/feature.sample.i.csv"))
#- 0d.1.2: Cadaver intensity tables from PeakWalk
bp1_cadaver <- read_csv(file.path(config$paths$peakwalk_cadaver_dir, "BP1.GC2/feature.sample.i.csv"))
bp2_cadaver <- read_csv(file.path(config$paths$peakwalk_cadaver_dir, "BP2.GC2/feature.sample.i.csv"))
bp3_cadaver <- read_csv(file.path(config$paths$peakwalk_cadaver_dir, "BP3.GC2/feature.sample.i.csv"))
#+ 0d.2: Bind Tables; Create a name_sub_lib_id Column
#- 0d.2.1: Tumors
combined_peakwalk_tumor <- bind_rows(bp1_tumor, bp2_tumor, bp3_tumor) |>
  select(-...1) |>
  mutate(id_subid = paste0(id, "_", subid))
#- 0d.2.2: Cadaver
combined_peakwalk_cadaver <- bind_rows(bp1_cadaver, bp2_cadaver, bp3_cadaver) |> 
  select(-...1) |>
  mutate(id_subid = paste0(id, "_", subid))
#+ 0d.3: Import RT Data
#- 0d.3.1: Tumors RT tables from PeakWalk
bp1_tumor_rt <- read_csv(file.path(config$paths$peakwalk_tumor_dir, "BP1.GC2/feature.sample.rt.csv"))
bp2_tumor_rt <- read_csv(file.path(config$paths$peakwalk_tumor_dir, "BP2.GC2/feature.sample.rt.csv"))
bp3_tumor_rt <- read_csv(file.path(config$paths$peakwalk_tumor_dir, "BP3.GC2/feature.sample.rt.csv"))
#- 0d.3.2: Cadaver RT tables from PeakWalk
bp1_cadaver_rt <- read_csv(file.path(config$paths$peakwalk_cadaver_dir, "BP1.GC2/feature.sample.rt.csv"))
bp2_cadaver_rt <- read_csv(file.path(config$paths$peakwalk_cadaver_dir, "BP2.GC2/feature.sample.rt.csv"))
bp3_cadaver_rt <- read_csv(file.path(config$paths$peakwalk_cadaver_dir, "BP3.GC2/feature.sample.rt.csv"))
#+ 0d.4: Bind RT Tables; Create a name_sub_lib_id Column
#- 0d.4.1: Tumors
combined_peakwalk_tumor_rt <- bind_rows(bp1_tumor_rt, bp2_tumor_rt, bp3_tumor_rt) |>
  select(-...1) |>
  mutate(id_subid = paste0(id, "_", subid))
#- 0d.4.2: Cadaver
combined_peakwalk_cadaver_rt <- bind_rows(bp1_cadaver_rt, bp2_cadaver_rt, bp3_cadaver_rt) |> 
  select(-...1) |>
  mutate(id_subid = paste0(id, "_", subid))
#+ 0d.5: Simplify RT Tables for Lookup
#- 0d.5.1: Tumors - pivot to long format with id_subid, file, and rt
tumor_rt_long <- combined_peakwalk_tumor_rt |>
  select(id_subid, starts_with("BL_"), starts_with("BP")) |>
  pivot_longer(
    cols = -id_subid,
    names_to = "file",
    values_to = "rt"
  ) |>
  filter(!is.na(rt) & rt != 0)
#- 0d.5.2: Cadaver - pivot to long format with id_subid, file, and rt
cadaver_rt_long <- combined_peakwalk_cadaver_rt |>
  select(id_subid, starts_with("BL_"), starts_with("BP")) |>
  pivot_longer(
    cols = -id_subid,
    names_to = "file",
    values_to = "rt"
  ) |>
  filter(!is.na(rt) & rt != 0)
