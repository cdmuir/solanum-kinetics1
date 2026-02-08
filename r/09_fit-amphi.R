# Analyze effect of anatomy on tau in untreated (amphistomatous) leaves
# Note: I did some preliminary model comparisons suggesting no random effects 
# of light_intensity, light_treatment, log_gcl, or log_fgmax
source("r/header.R")

joined_summary = read_rds("data/joined-summary.rds") |>
  prepare_tau_anatomy_data(logtau_threshold)

phy = read_rds("data/phylogeny.rds")
A = vcv(phy, corr = TRUE)
thin = 4

# Define formula
bf1 = bf(loglambdamean | se(loglambdasd, sigma = TRUE) ~ lighttreatment + (1|a|accession) + (1|b|gr(phy, cov = A)))
bf2 = bf(logtaumean | se(logtausd, sigma = TRUE) ~ lighttreatment + (1|a|accession) + (1|b|gr(phy, cov = A)))
bf3 = bf(loggcl ~ lighttreatment + (1|a|accession) + (1|b|gr(phy, cov = A)))
bf4 = bf(logfgmax ~ lighttreatment + (1|a|accession) + (1|b|gr(phy, cov = A)))

# Fit model under high measurement light intensity
fit_amphi_high = brm(
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
) 

# Fit model under low measurement light intensity
fit_amphi_low = brm(
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
) 

write_rds(fit_amphi_high, "objects/fit_amphi_high.rds")
write_rds(fit_amphi_low, "objects/fit_amphi_low.rds")
