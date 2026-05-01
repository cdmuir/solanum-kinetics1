# Decompose variance of log(gcl) into:
# - within leaf + error
# - among individual within population
# - among population (nonphylogenetic)
# - among population (phylogeny)

source("r/header.R")

plant_info = read_rds("data/plant_info.rds") |>
  select(accession, replicate, light_treatment)

gcl = read_rds("data/guard-cell-length.rds") |>
  ungroup() |>
  mutate(loggcl = log(guard_cell_length_um)) 

phy = read_rds("data/phylogeny.rds")
A = vcv(phy, corr = TRUE)

# Filter guard cell data to only include accessions and replicates present in
# joined_data
joined_data = read_rds("data/joined-data.rds") |>
  select(accession = acc, replicate = id) |>
  distinct() |>
  mutate(use = TRUE)

gcl1 = gcl |>
  full_join(joined_data, by = join_by(accession, replicate)) |>
  filter(!is.na(use)) |>
  select(-use) |>
  left_join(plant_info, by = join_by(accession, replicate)) |>
  unite("acc_id",
        accession,
        replicate,
        sep = "-",
        remove = FALSE) |>
  mutate(phy = accession) |>
  # Trim most extreme 1%
  filter(loggcl > quantile(loggcl, 0.01 / 2),
         loggcl < quantile(loggcl, 1 - 0.01 / 2))

thin = 4

fit_gcl = brm(
  loggcl ~ light_treatment * surface + 
    (1 + | acc_id) + 
    (1 + | accession) + 
    (1 + | gr(phy, cov = A)),
  data = gcl1,
  data2 = list(A = A),
  cores = 4,
  chains = 4,
  iter = thin * 2e3,
  thin = thin,
  refresh = thin * 1e2,
  # control = list(adapt_delta = 0.9),
  backend = "cmdstanr",
  # family = student(),
  seed = 613135062
) |> 
  add_criterion("loo")

write_rds(fit_gcl, "objects/fit-gcl.rds")
