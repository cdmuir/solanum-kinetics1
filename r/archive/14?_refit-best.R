# Refit best model with VPD as covariate - not currently implementing
source("r/header.R")

plan(multisession, workers = 16)

joined_summary = read_rds("data/joined-summary.rds") |>
  prepare_tau_anatomy_data(logtau_threshold)

write_rds(list(logtau_threshold = logtau_threshold, n_removed = attr(joined_summary, "n_removed")), "objects/n_removed.rds")

curve_vpd = read_rds("objects/curve-vpd.rds")

joined_summary1 = joined_summary |>
  left_join(curve_vpd, by = join_by(curvetype, lightintensity, accid, lighttreatment)) |>
  mutate(vpd = scale(VPDleaf)[,1])

assert_true(all(!is.na(joined_summary1$loglambdamean)))
assert_true(all(!is.na(joined_summary1$loglambdasd)))
assert_true(all(!is.na(joined_summary1$logtausd)))
assert_true(all(!is.na(joined_summary1$logtausd)))
assert_true(all(!is.na(joined_summary1$loggcl)))
assert_true(all(!is.na(joined_summary1$logitfgmax)))

assert_true(all(!is.na(joined_summary1$lighttreatment)))
assert_true(all(!is.na(joined_summary1$lightintensity)))
assert_true(all(!is.na(joined_summary1$curvetype)))
assert_true(all(!is.na(joined_summary1$VPDleaf)))

phy = read_rds("data/phylogeny.rds")
A = vcv(phy, corr = TRUE)
thin = 6 #12

# Define formula
bf_lambda0 = bf(loglambdamean | se(loglambdasd, sigma = TRUE) ~ 
                  lighttreatment + 
                  lightintensity +
                  curvetype +
                  vpd +
                  loggcl + 
                  logitfgmax + 
                  (1|accid) +
                  (1|a|accession) +
                  (1|b|gr(phy, cov = A)))
bf_lambda1 = update(bf_lambda0, . ~ . - loggcl)
bf_lambda2 = update(bf_lambda0, . ~ . - logitfgmax)
bf_lambda3 = update(bf_lambda0, . ~ . - loggcl - logitfgmax)

bf_tau0 = bf(logtaumean | se(logtausd, sigma = TRUE) ~ 
               lighttreatment + 
               lightintensity +
               curvetype +
               vpd +
               loggcl + 
               logitfgmax + 
               (1|accid) +
               (1|a|accession) +
               (1|b|gr(phy, cov = A)))

bf_tau1 = update(bf_tau0, . ~ . - loggcl)
bf_tau2 = update(bf_tau0, . ~ . - logitfgmax)
bf_tau3 = update(bf_tau0, . ~ . - loggcl - logitfgmax)

bf_gcl = bf(loggcl ~ 
              lighttreatment + 
              curvetype +
              (1|a|accession) + 
              (1|b|gr(phy, cov = A)))

bf_fgmax = bf(logitfgmax ~ 
                lighttreatment + 
                lightintensity + 
                curvetype +
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
        data = joined_summary1 |>
          mutate(phy = accession),
        data2 = list(A = A),
        cores = 1,
        chains = 4,
        iter = thin * 2e3,
        thin = thin,
        refresh = thin * 1e2,
        control = list(adapt_delta = 0.9),
        backend = "cmdstanr",
        family = student(),
        seed = 613135062 + i
      ) |> add_criterion("loo")
    }, .options = furrr_options(seed = TRUE), .progress = TRUE)
  )

write_rds(fits, "objects/fits.rds")
