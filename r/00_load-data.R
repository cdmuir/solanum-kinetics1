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
  select(-contains("pavement"), -contains("index")) |>
  write_rds("data/stomata.rds")

read_rds("https://github.com/cdmuir/solanum-aa/raw/refs/heads/main/data/df_germ_summary.rds") |>
  write_rds("data/df_germ_summary.rds")

read_rds("https://github.com/cdmuir/solanum-aa/raw/refs/heads/main/data/df_growth_summary.rds") |>
  write_rds("data/df_growth_summary.rds")
