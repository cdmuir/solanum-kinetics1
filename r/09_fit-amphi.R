# Analyze effect of anatomy on tau in untreated (amphistomatous) leaves
# Note: I did some preliminary model comparisons suggesting no random effects 
# of light_intensity, light_treatment, log_gcl, or log_fgmax
source("r/header.R")

joined_summary = read_rds("data/joined-summary.rds") |>
  prepare_tau_anatomy_data(logtau_threshold)

phy = read_rds("data/phylogeny.rds")
A = vcv(phy, corr = TRUE)
thin = 4

# Model 0 ----
# fixed effect of lighttreatment, random effects of accession and phylogeny

# Define formula
bf1 = bf(loglambdamean | se(loglambdasd, sigma = TRUE) ~ lighttreatment + (1|a|accession) + (1|b|gr(phy, cov = A)))
bf2 = bf(logtaumean | se(logtausd, sigma = TRUE) ~ lighttreatment + (1|a|accession) + (1|b|gr(phy, cov = A)))
bf3 = bf(loggcl ~ lighttreatment + (1|a|accession) + (1|b|gr(phy, cov = A)))
bf4 = bf(logfgmax ~ lighttreatment + (1|a|accession) + (1|b|gr(phy, cov = A)))

# Fit model under high measurement light intensity
fit_amphi_high0 = brm(
  bf1 + bf2 + bf3 + bf4 + set_rescor(TRUE),
  data = joined_summary |>
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
  bf1 + bf2 + bf3 + bf4 + set_rescor(TRUE),
  data = joined_summary |>
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
# fixed effect of lighttreatment, random effect of phylogeny only

# Define formula
bf1 = bf(loglambdamean | se(loglambdasd, sigma = TRUE) ~ lighttreatment + (1|b|gr(phy, cov = A)))
bf2 = bf(logtaumean | se(logtausd, sigma = TRUE) ~ lighttreatment + (1|b|gr(phy, cov = A)))
bf3 = bf(loggcl ~ lighttreatment + (1|b|gr(phy, cov = A)))
bf4 = bf(logfgmax ~ lighttreatment + (1|b|gr(phy, cov = A)))

# Fit model under high measurement light intensity
fit_amphi_high1 = brm(
  bf1 + bf2 + bf3 + bf4 + set_rescor(TRUE),
  data = joined_summary |>
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
  bf1 + bf2 + bf3 + bf4 + set_rescor(TRUE),
  data = joined_summary |>
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
# fixed effect of lighttreatment, random effect of accession only

# Define formula
bf1 = bf(loglambdamean | se(loglambdasd, sigma = TRUE) ~ lighttreatment + (1|a|accession))
bf2 = bf(logtaumean | se(logtausd, sigma = TRUE) ~ lighttreatment + (1|a|accession))
bf3 = bf(loggcl ~ lighttreatment + (1|a|accession))
bf4 = bf(logfgmax ~ lighttreatment + (1|a|accession))

# Fit model under high measurement light intensity
fit_amphi_high2 = brm(
  bf1 + bf2 + bf3 + bf4 + set_rescor(TRUE),
  data = joined_summary |>
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
  bf1 + bf2 + bf3 + bf4 + set_rescor(TRUE),
  data = joined_summary |>
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
# fixed effect of lighttreatment only

# Define formula
bf1 = bf(loglambdamean | se(loglambdasd, sigma = TRUE) ~ lighttreatment)
bf2 = bf(logtaumean | se(logtausd, sigma = TRUE) ~ lighttreatment)
bf3 = bf(loggcl ~ lighttreatment)
bf4 = bf(logfgmax ~ lighttreatment)

# Fit model under high measurement light intensity
fit_amphi_high3 = brm(
  bf1 + bf2 + bf3 + bf4 + set_rescor(TRUE),
  data = joined_summary |>
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
  bf1 + bf2 + bf3 + bf4 + set_rescor(TRUE),
  data = joined_summary |>
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

# Save model with phylogenetic effects only - models 0-2 all fit similarly, but
# model 1 most clearly addresses hypotheses. Model 3 fits much worse, so discard.
write_rds(fit_amphi_high1, "objects/fit_amphi_high.rds")
write_rds(fit_amphi_low1, "objects/fit_amphi_low.rds")


# This model actually fits better according to LOOIC, but not sure if clarifies anything. Keep here for now, but delete if not using.

bf1 = bf(logtaumean | se(logtausd, sigma = TRUE) ~ lighttreatment + (lighttreatment|a|accession) + (lighttreatment|b|gr(phy, cov = A)))
bf2 = bf(loggcl ~ lighttreatment + (lighttreatment|a|accession) + (lighttreatment|b|gr(phy, cov = A)))
bf3 = bf(logfgmax ~ lighttreatment + (lighttreatment|a|accession) + (lighttreatment|b|gr(phy, cov = A)))

# Fit model under high measurement light intensity
# fit_amphi_highX = brm(
#   bf1 + bf2 + bf3 + set_rescor(TRUE),
#   data = joined_summary |>
#     filter(curve_type == "amphi", lightintensity == "high") |>
#     mutate(phy = accession),
#   data2 = list(A = A),
#   cores = 1,
#   chains = 1,
#   iter = thin * 2e3,
#   thin = thin,
#   refresh = thin * 1e2,
#   backend = "cmdstanr",
#   seed = 613135062,
# ) |>
#   add_criterion("loo")
# 
# fit_amphi_lowX = brm(
#   bf1 + bf2 + bf3 + set_rescor(TRUE),
#   data = joined_summary |>
#     filter(curve_type == "amphi", lightintensity == "low") |>
#     mutate(phy = accession),
#   data2 = list(A = A),
#   cores = 1,
#   chains = 1,
#   iter = thin * 2e3,
#   thin = thin,
#   refresh = thin * 1e2,
#   backend = "cmdstanr",
#   seed = 613135062,
# ) |>
#   add_criterion("loo")
# 
#   