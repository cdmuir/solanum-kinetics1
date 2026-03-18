# Analyze effect of anatomy on tau in untreated (amphistomatous) leaves
# Note: I did some preliminary model comparisons suggesting no random effects
# of light_intensity, light_treatment, log_gcl, or log_fgmax
source("r/header.R")

joined_summary = read_rds("data/joined-summary.rds") |>
  prepare_tau_anatomy_data(logtau_threshold)

phy = read_rds("data/phylogeny.rds")
A = vcv(phy, corr = TRUE)
thin = 1 #12

# Define formula
bf_lambda0 = bf(
  loglambdamean | se(loglambdasd, sigma = TRUE) ~
    lighttreatment +
    lightintensity +
    curve_type +
    loggcl +
    logeca +
    logitfgmax +
    (1 | accid) +
    (1 | a | accession) +
    (1 | b | gr(phy, cov = A))
)
bf_lambda1 = update(bf_lambda0, . ~ . - loggcl)
bf_lambda2 = update(bf_lambda0, . ~ . - logitfgmax)
bf_lambda3 = update(bf_lambda0, . ~ . - loggcl - logitfgmax)

bf_tau0 = bf(
  logtaumean | se(logtausd, sigma = TRUE) ~
    lighttreatment +
    lightintensity +
    curve_type +
    loggcl +
    logeca +
    logitfgmax +
    (1 | accid) +
    (1 | a | accession) +
    (1 | b | gr(phy, cov = A))
)

bf_tau1 = update(bf_tau0, . ~ . - loggcl)
bf_tau2 = update(bf_tau0, . ~ . - logitfgmax)
bf_tau3 = update(bf_tau0, . ~ . - loggcl - logitfgmax)

bf_gcl = bf(loggcl ~
              lighttreatment +
              curve_type +
              (1 | a | accession) +
              (1 | b | gr(phy, cov = A)))

bf_eca = bf(logeca | mi() ~
              lighttreatment +
              curve_type +
              (1 | a | accession) +
              (1 | b | gr(phy, cov = A)))

bf_fgmax = bf(
  logitfgmax ~
    lighttreatment +
    lightintensity +
    curve_type +
    (1 | accid) +
    (1 | a | accession) +
    (1 | b | gr(phy, cov = A))
)

fit = brm(
  formula = bf_tau1 + bf_lambda1 + bf_gcl + bf_eca + bf_fgmax + set_rescor(TRUE) ,
  data = joined_summary |>
    mutate(phy = accession),
  data2 = list(A = A),
  cores = 4,
  chains = 4,
  iter = thin * 2e3,
  thin = thin,
  refresh = thin * 1e2,
  control = list(adapt_delta = 0.9),
  backend = "cmdstanr",
  family = student(),
  seed = 613135062
) |> add_criterion("loo")


write_rds(fit, "objects/fit-test.rds")
