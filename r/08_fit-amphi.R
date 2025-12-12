# Analyze effect of anatomy on tau in untreated (amphistomatous) leaves
# Note: I did some preliminary model comparisons suggesting no random effects 
# of light_intensity, light_treatment, log_gcl, or log_fgmax
source("r/header.R")

joined_data = read_rds("data/joined-data.rds") |>
  # need to do more careful outlier removal
  filter(log_tau_mean < 7) |>
  mutate(loggcl = log(guard_cell_length_um),
         logfgmax = log(f_gmax)) |>
  rename(
    logtaumean = log_tau_mean,
    logtausd = log_tau_sd,
    accid = acc_id,
    lightintensity = light_intensity,
    lighttreatment = light_treatment
  )

phy = read_rds("data/phylogeny.rds")
A = vcv(phy, corr = TRUE)
thin = 1

# Model 0 ----

# Define formula
bf1 = bf(logtaumean | se(logtausd, sigma = TRUE) ~ lighttreatment + (1|a|accession) + (1|b|gr(phy, cov = A)))
bf2 = bf(loggcl ~ lighttreatment + (1|a|accession) + (1|b|gr(phy, cov = A)))
bf3 = bf(logfgmax ~ lighttreatment + (1|a|accession) + (1|b|gr(phy, cov = A)))

# Fit model under high measurement light intensity
fit_amphi_high0 = brm(
  bf1 + bf2 + bf3 + set_rescor(TRUE),
  data = joined_data |>
    filter(curve_type == "amphi", lightintensity == "high") |>
    mutate(phy = accession),
  data2 = list(A = A),
  cores = 1,
  chains = 1,
  iter = thin * 2e3,
  thin = thin,
  refresh = thin * 1e2,
  backend = "cmdstanr",
  seed = 613135062,
) |>
  add_criterion("loo")

# Fit model under low measurement light intensity
fit_amphi_low0 = brm(
  bf1 + bf2 + bf3 + set_rescor(TRUE),
  data = joined_data |>
    filter(curve_type == "amphi", lightintensity == "low") |>
    mutate(phy = accession),
  data2 = list(A = A),
  cores = 1,
  chains = 1,
  iter = thin * 2e3,
  thin = thin,
  refresh = thin * 1e2,
  backend = "cmdstanr",
  seed = 261582805
) |>
  add_criterion("loo")

# Model 1 ----

# Define formula
bf1 = bf(logtaumean | se(logtausd, sigma = TRUE) ~ lighttreatment + (1|b|gr(phy, cov = A)))
bf2 = bf(loggcl ~ lighttreatment + (1|b|gr(phy, cov = A)))
bf3 = bf(logfgmax ~ lighttreatment + (1|b|gr(phy, cov = A)))

# Fit model under high measurement light intensity
fit_amphi_high1 = brm(
  bf1 + bf2 + bf3 + set_rescor(TRUE),
  data = joined_data |>
    filter(curve_type == "amphi", lightintensity == "high") |>
    mutate(phy = accession),
  data2 = list(A = A),
  cores = 1,
  chains = 1,
  iter = thin * 2e3,
  thin = thin,
  refresh = thin * 1e2,
  backend = "cmdstanr",
  seed = 613135062,
) |>
  add_criterion("loo")

# Fit model under low measurement light intensity
fit_amphi_low1 = brm(
  bf1 + bf2 + bf3 + set_rescor(TRUE),
  data = joined_data |>
    filter(curve_type == "amphi", lightintensity == "low") |>
    mutate(phy = accession),
  data2 = list(A = A),
  cores = 1,
  chains = 1,
  iter = thin * 2e3,
  thin = thin,
  refresh = thin * 1e2,
  backend = "cmdstanr",
  seed = 261582805
) |>
  add_criterion("loo")

# Model 2 ----

# Define formula
bf1 = bf(logtaumean | se(logtausd, sigma = TRUE) ~ lighttreatment + (1|a|accession))
bf2 = bf(loggcl ~ lighttreatment + (1|a|accession))
bf3 = bf(logfgmax ~ lighttreatment + (1|a|accession))

# Fit model under high measurement light intensity
fit_amphi_high2 = brm(
  bf1 + bf2 + bf3 + set_rescor(TRUE),
  data = joined_data |>
    filter(curve_type == "amphi", lightintensity == "high") |>
    mutate(phy = accession),
  data2 = list(A = A),
  cores = 1,
  chains = 1,
  iter = thin * 2e3,
  thin = thin,
  refresh = thin * 1e2,
  backend = "cmdstanr",
  seed = 613135062,
) |>
  add_criterion("loo")

# Fit model under low measurement light intensity
fit_amphi_low2 = brm(
  bf1 + bf2 + bf3 + set_rescor(TRUE),
  data = joined_data |>
    filter(curve_type == "amphi", lightintensity == "low") |>
    mutate(phy = accession),
  data2 = list(A = A),
  cores = 1,
  chains = 1,
  iter = thin * 2e3,
  thin = thin,
  refresh = thin * 1e2,
  backend = "cmdstanr",
  seed = 261582805
) |>
  add_criterion("loo")

# Model 3 ----

# Define formula
bf1 = bf(logtaumean | se(logtausd, sigma = TRUE) ~ lighttreatment)
bf2 = bf(loggcl ~ lighttreatment)
bf3 = bf(logfgmax ~ lighttreatment)

# Fit model under high measurement light intensity
fit_amphi_high3 = brm(
  bf1 + bf2 + bf3 + set_rescor(TRUE),
  data = joined_data |>
    filter(curve_type == "amphi", lightintensity == "high") |>
    mutate(phy = accession),
  data2 = list(A = A),
  cores = 1,
  chains = 1,
  iter = thin * 2e3,
  thin = thin,
  refresh = thin * 1e2,
  backend = "cmdstanr",
  seed = 613135062,
) |>
  add_criterion("loo")

# Fit model under low measurement light intensity
fit_amphi_low3 = brm(
  bf1 + bf2 + bf3 + set_rescor(TRUE),
  data = joined_data |>
    filter(curve_type == "amphi", lightintensity == "low") |>
    mutate(phy = accession),
  data2 = list(A = A),
  cores = 1,
  chains = 1,
  iter = thin * 2e3,
  thin = thin,
  refresh = thin * 1e2,
  backend = "cmdstanr",
  seed = 261582805
) |>
  add_criterion("loo")

# Write model objects
loo_compare(fit_amphi_high0, fit_amphi_high1, fit_amphi_high2, fit_amphi_high3)
loo_compare(fit_amphi_low0, fit_amphi_low1, fit_amphi_low2, fit_amphi_low3)

# Save model with phylogenetic effects only
write_rds(fit_amphi_high1, "objects/fit_amphi_high.rds")
write_rds(fit_amphi_low1, "objects/fit_amphi_low.rds")
