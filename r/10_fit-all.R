# --- Fit the multiresponse phylogenetic model relating anatomy to kinetics ---
#
# Note: preliminary model comparisons found no support for random effects of
# light_intensity, light_treatment, log_gcl, log_gi, or log_gmax, so none are
# included here.

source("r/header.R")

plan(multisession, workers = 16)

joined_summary = read_rds("data/joined-summary.rds") |>
  prepare_tau_anatomy_data(logtau_threshold)

write_rds(
  list(
    logtau_threshold = logtau_threshold,
    n_removed = attr(joined_summary, "n_removed")
  ),
  "objects/n_removed.rds"
)

assert_true(all(!is.na(joined_summary$loglambdamean)))
assert_true(all(!is.na(joined_summary$loglambdasd)))
assert_true(all(!is.na(joined_summary$logtausd)))
assert_true(all(!is.na(joined_summary$logtausd)))
assert_true(all(!is.na(joined_summary$loggcl)))
assert_true(all(!is.na(joined_summary$loggi)))
assert_true(all(!is.na(joined_summary$loggmax)))

assert_true(all(!is.na(joined_summary$lighttreatment)))
assert_true(all(!is.na(joined_summary$lightintensity)))
assert_true(all(!is.na(joined_summary$leaftype)))

phy = read_rds("data/phylogeny.rds")
A = vcv(phy, corr = TRUE)
thin = 6

# Define formula
bf_lambda0 = bf(
  loglambdamean | se(loglambdasd, sigma = TRUE) ~
    lighttreatment +
    lightintensity +
    leaftype +
    loggcl +
    loggi +
    loggmax +
    (1 | accid) +
    (1 | a | accession) +
    (1 | b | gr(phy, cov = A))
)

bf_lambda1 = update(bf_lambda0, . ~ . - loggcl)
bf_lambda2 = update(bf_lambda0, . ~ . - loggi)
bf_lambda3 = update(bf_lambda0, . ~ . - loggmax)
bf_lambda4 = update(bf_lambda0, . ~ . - loggcl - loggi)
bf_lambda5 = update(bf_lambda0, . ~ . - loggcl - loggmax)
bf_lambda6 = update(bf_lambda0, . ~ . - loggi - loggmax)
bf_lambda7 = update(bf_lambda0, . ~ . - loggcl - loggi - loggmax)

bf_tau0 = bf(
  logtaumean | se(logtausd, sigma = TRUE) ~
    lighttreatment +
    lightintensity +
    leaftype +
    loggcl +
    loggi +
    loggmax +
    (1 | accid) +
    (1 | a | accession) +
    (1 | b | gr(phy, cov = A))
)

bf_tau1 = update(bf_tau0, . ~ . - loggcl)
bf_tau2 = update(bf_tau0, . ~ . - loggi)
bf_tau3 = update(bf_tau0, . ~ . - loggmax)
bf_tau4 = update(bf_tau0, . ~ . - loggcl - loggi)
bf_tau5 = update(bf_tau0, . ~ . - loggcl - loggmax)
bf_tau6 = update(bf_tau0, . ~ . - loggi - loggmax)
bf_tau7 = update(bf_tau0, . ~ . - loggcl - loggi - loggmax)

bf_gcl = bf(loggcl ~
              lighttreatment +
              leaftype +
              (1 | a | accession) +
              (1 | b | gr(phy, cov = A)))

bf_gi = bf(
  loggi ~
    lighttreatment +
    lightintensity +
    leaftype +
    (1 | accid) +
    (1 | a | accession) +
    (1 | b | gr(phy, cov = A))
)

bf_gmax = bf(loggmax ~
               lighttreatment +
               leaftype +
               (1 | a | accession) +
               (1 | b | gr(phy, cov = A)))

df_forms = crossing(
  bf_lambda = list(
    bf_lambda0,
    bf_lambda1,
    bf_lambda2,
    bf_lambda3,
    bf_lambda4,
    bf_lambda5,
    bf_lambda6,
    bf_lambda7
  ),
  bf_tau = list(
    bf_tau0,
    bf_tau1,
    bf_tau2,
    bf_tau3,
    bf_tau4,
    bf_tau5,
    bf_tau6,
    bf_tau7
  )
) |>
  mutate(
    form = map2(
      bf_lambda,
      bf_tau,
      ~ .x + .y + bf_gcl + bf_gi + bf_gmax + set_rescor(TRUE)
    ),
    file = glue("objects/fits/fit_{n}.rds", n = str_pad(row_number(), 2, "left", 0))
  )

future_walk2(
  df_forms$form,
  df_forms$file,
  \(.form, .file) {
    fit = brm(
      formula = .form,
      data = joined_summary |>
        mutate(phy = accession),
      data2 = list(A = A),
      cores = 1,
      chains = 3,
      iter = thin * 2e3,
      thin = thin,
      refresh = thin * 1e2,
      control = list(adapt_delta = 0.9),
      backend = "cmdstanr",
      family = student(),
      seed = 613135062 + as.numeric(str_extract(.file, "[0-9]+"))
    ) |> add_criterion("loo")
    
    write_rds(fit, .file)
    
  },
  .options = furrr_options(seed = TRUE),
  .progress = TRUE
)

write_rds(df_forms, "objects/df_forms.rds")
