# Download and preprocess data
source("r/header.R")

read_rds("https://github.com/cdmuir/solanum-aa/raw/refs/heads/main/data/accession-info.rds") |>
  write_rds("data/accession_info.rds")

read_rds("https://github.com/cdmuir/solanum-aa/raw/refs/heads/main/data/plant-info.rds") |>
  write_rds("data/plant_info.rds")

read_rds("https://github.com/cdmuir/solanum-aa/raw/refs/heads/main/data/phylogeny.rds") |>
  write_rds("data/phylogeny.rds")

read_rds("https://github.com/cdmuir/solanum-aa/raw/refs/heads/main/data/trimmed_rh_curves.rds") |>
  mutate(ci = as.numeric(as.factor(curve))) |>
  mutate(t_sec = elapsed - min(elapsed), .by = ci) |>
  write_rds("data/rh_curves.rds")

read_rds("https://github.com/cdmuir/solanum-aa/raw/refs/heads/main/data/stomata.rds") |>
  select(-contains("pavement")) |>
  # Add estimate of average epidermal cell area based on Sack and Buckley (2016)
  mutate(
    lower_guard_cell_area_um2 = lower_guard_cell_length_um^2 / 2,
    upper_guard_cell_area_um2 = upper_guard_cell_length_um^2 / 2,
    lower_stomatal_density_um2 = lower_stomatal_density_mm2 / 1e6,
    upper_stomatal_density_um2 = upper_stomatal_density_mm2 / 1e6,
    lower_epidermal_cell_area_um2 = lower_stomatal_index * (1 - lower_stomatal_density_um2 * lower_guard_cell_area_um2) / (lower_stomatal_density_um2 * (1 - lower_stomatal_index)),
    upper_epidermal_cell_area_um2 = upper_stomatal_index * (1 - upper_stomatal_density_um2 * upper_guard_cell_area_um2) / (upper_stomatal_density_um2 * (1 - upper_stomatal_index))
  ) |>
  write_rds("data/stomata.rds")

file.copy("../../data/adaptive-amphistomy/processed-data/guard-cell-length.rds",
          "data/guard-cell-length.rds")

read_rds("https://github.com/cdmuir/solanum-aa/raw/refs/heads/main/data/df_germ_summary.rds") |>
  write_rds("data/df_germ_summary.rds")

read_rds("https://github.com/cdmuir/solanum-aa/raw/refs/heads/main/data/df_growth_summary.rds") |>
  write_rds("data/df_growth_summary.rds")
