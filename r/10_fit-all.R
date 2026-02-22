# Analyze effect of anatomy on tau in untreated (amphistomatous) leaves
# Note: I did some preliminary model comparisons suggesting no random effects 
# of light_intensity, light_treatment, log_gcl, or log_fgmax
source("r/header.R")

plan(multisession, workers = 16)

joined_summary = read_rds("data/joined-summary.rds") |>
  prepare_tau_anatomy_data(logtau_threshold)

write_rds(list(logtau_threshold = logtau_threshold, n_removed = attr(joined_summary, "n_removed")), "objects/n_removed.rds")

phy = read_rds("data/phylogeny.rds")
A = vcv(phy, corr = TRUE)
thin = 9

# Define formula
bf_lambda0 = bf(loglambdamean | se(loglambdasd, sigma = TRUE) ~ 
                  lighttreatment + 
                  lightintensity +
                  curve_type +
                  loggcl + 
                  logfgmax + 
                  (1|accid) +
                  (1|a|accession) +
                  (1|b|gr(phy, cov = A)))
bf_lambda1 = update(bf_lambda0, . ~ . - loggcl)
bf_lambda2 = update(bf_lambda0, . ~ . - logfgmax)
bf_lambda3 = update(bf_lambda0, . ~ . - loggcl - logfgmax)

bf_tau0 = bf(logtaumean | se(logtausd, sigma = TRUE) ~ 
               lighttreatment + 
               lightintensity +
               curve_type +
               loggcl + 
               logfgmax + 
               (1|accid) +
               (1|a|accession) +
               (1|b|gr(phy, cov = A)))

bf_tau1 = update(bf_tau0, . ~ . - loggcl)
bf_tau2 = update(bf_tau0, . ~ . - logfgmax)
bf_tau3 = update(bf_tau0, . ~ . - loggcl - logfgmax)

bf_gcl = bf(loggcl ~ 
              lighttreatment + 
              curve_type +
              (1|a|accession) + 
              (1|b|gr(phy, cov = A)))

bf_fgmax = bf(logfgmax ~ 
                lighttreatment + 
                lightintensity + 
                curve_type +
                (1|accid) +
                (1|a|accession) + 
                (1|b|gr(phy, cov = A)))

fits = crossing(
  bf_lambda = list(bf_lambda0, bf_lambda1, bf_lambda2, bf_lambda3),
  bf_tau = list(bf_tau0, bf_tau1, bf_tau2, bf_tau3)
) |>
  mutate(
    form = map2(bf_lambda, bf_tau, ~ .x + .y + bf_gcl + bf_fgmax + set_rescor(TRUE)),
    fit = future_map2(form, row_number(), \(.form, i) {
      brm(
        formula = .form,
        data = joined_summary |>
          mutate(phy = accession),
        data2 = list(A = A),
        cores = 1,
        chains = 1,
        iter = thin * 2e3,
        thin = thin,
        refresh = thin * 1e2,
        control = list(adapt_delta = 0.9),
        backend = "cmdstanr",
        seed = 613135062 + i
      ) |> add_criterion("loo")
    }, .options = furrr_options(seed = TRUE), .progress = TRUE)
  )

write_rds(fits, "objects/fits.rds")
