# Analyze effect of anatomy on tau in treated (pseudohypostomatous) leaves
source("r/header.R")

# Use same model structure as amphi for direct comparison
fit_amphi_high = read_rds("objects/fit_amphi_high.rds")
fit_amphi_low = read_rds("objects/fit_amphi_low.rds")

joined_summary = read_rds("data/joined-summary.rds") |>
  prepare_tau_anatomy_data(logtau_threshold)

phy = read_rds("data/phylogeny.rds")
A = vcv(phy, corr = TRUE)
thin = 3

# Fit model under high measurement light intensity
fit_pseudohypo_high = brm(
  formula(fit_amphi_high),
  data = joined_summary |>
    filter(curve_type == "pseudohypo", lightintensity == "high") |>
    mutate(phy = accession),
  data2 = list(A = A),
  cores = 1,
  chains = 1,
  iter = thin * 2e3,
  thin = thin,
  refresh = thin * 1e2,
  backend = "cmdstanr",
  seed = 850704602,
)

# Fit model under low measurement light intensity
fit_pseudohypo_low = brm(
  formula(fit_amphi_low),
  data = joined_summary |>
    filter(curve_type == "pseudohypo", lightintensity == "low") |>
    mutate(phy = accession),
  data2 = list(A = A),
  cores = 1,
  chains = 1,
  iter = thin * 2e3,
  thin = thin,
  refresh = thin * 1e2,
  backend = "cmdstanr",
  seed = 789804965
)

write_rds(fit_pseudohypo_high, "objects/fit_pseudohypo_high.rds")
write_rds(fit_pseudohypo_low, "objects/fit_pseudohypo_low.rds")
