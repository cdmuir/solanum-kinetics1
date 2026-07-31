# Validate that the fast nls() re-fit used in r/32_null-sim-fgmax-tau.R is a
# reasonable approximation to the production Bayesian pipeline
# (fit_rh1()/bform_cdweibull/pri from r/02_fit-weibull.R and r/functions.R),
# by re-fitting a subsample of simulated curves both ways and comparing.

source("r/header.R")
library(brms)

n_validation = 50
validation_seed = 9182

example_fits = read_rds("objects/null-sim-fgmax-tau-example.rds")

# Regenerate the same example simulated dataset used to produce
# null-sim-fgmax-tau-example.rds (see r/32_null-sim-fgmax-tau.R) so we have
# the raw simulated (t_sec, gsw) observations, not just the nls point
# estimates, to refit with brms.
pars_summary = read_rds("objects/pars-summary.rds")
joined_data = read_rds("data/joined-data.rds")

real_pars = pars_summary |>
  filter(variable %in% c("b_logtau_Intercept", "b_loglambda_Intercept", "sigma", "ginit", "gfinal")) |>
  select(id, variable, mean) |>
  pivot_wider(names_from = variable, values_from = mean) |>
  transmute(id, tau = exp(b_logtau_Intercept), lambda = exp(b_loglambda_Intercept),
            sigma, gi = ginit, gf = gfinal)

time_design = joined_data |> select(curve, t_sec)
median_gi = median(real_pars$gi, na.rm = TRUE)
median_gf = median(real_pars$gf, na.rm = TRUE)

set.seed(4471) # same seed used for the example dataset in r/32
example_sim_data = real_pars |>
  left_join(time_design, by = c("id" = "curve")) |>
  mutate(
    gsw_true = median_gf + (median_gi - median_gf) * exp(-(t_sec / tau) ^ lambda),
    gsw_sim = gsw_true + rnorm(n(), 0, sigma)
  )

set.seed(validation_seed)
validation_ids = sample(unique(example_sim_data$id), n_validation)

# --- Bayesian re-fit (same functional form/priors as production) ---------

form_cdweibull = gsw ~ exp(loggf) + exp(logdg) * exp(-(t_sec / exp(logtau))^exp(loglambda))
bform_cdweibull = bf(form_cdweibull, loggf ~ 1, logdg ~ 1, logtau ~ 1, loglambda ~ 1, nl = TRUE)

pri = c(
  prior(normal(log(0.1), 2), nlpar = "loggf", lb = log(0.001), ub = log(2)),
  prior(normal(log(0.5), 2), nlpar = "logdg", lb = log(0.001), ub = log(2)),
  prior(normal(log(300), 1), nlpar = "logtau", lb = log(10), ub = log(10000)),
  prior(normal(log(2), 1), nlpar = "loglambda", lb = log(0.1), ub = log(10))
)

# Fit the first curve fully (compiles the Stan model), then reuse the
# compiled model via update() for the rest -- much faster than recompiling
# per curve.
base_df = example_sim_data |> filter(id == validation_ids[1]) |> rename(gsw = gsw_sim)
base_fit = brm(
  formula = bform_cdweibull, data = base_df, prior = pri,
  iter = 2000, thin = 1, chains = 4, cores = 4,
  backend = "cmdstanr", control = list(adapt_delta = 0.8),
  seed = 1234, refresh = 0
)

bayes_fits = map(validation_ids, \(cid) {
  df = example_sim_data |> filter(id == cid) |> rename(gsw = gsw_sim)
  tryCatch(update(base_fit, newdata = df, recompile = FALSE, refresh = 0), error = function(e) NULL)
})
names(bayes_fits) = validation_ids

bayes_summary = map_dfr(names(bayes_fits), \(cid) {
  f = bayes_fits[[cid]]
  if (is.null(f)) return(tibble(id = cid, gi_bayes = NA, gf_bayes = NA, tau_bayes = NA, lambda_bayes = NA))
  cf = fixef(f)
  tibble(
    id = cid,
    gi_bayes = exp(cf["loggf_Intercept", "Estimate"]) + exp(cf["logdg_Intercept", "Estimate"]),
    gf_bayes = exp(cf["loggf_Intercept", "Estimate"]),
    tau_bayes = exp(cf["logtau_Intercept", "Estimate"]),
    lambda_bayes = exp(cf["loglambda_Intercept", "Estimate"])
  )
})

# --- Compare to nls() point estimates -------------------------------------

nls_summary = example_fits |>
  filter(id %in% validation_ids) |>
  select(id, gi_hat, gf_hat, tau_hat, lambda_hat)

nls_vs_bayes = bayes_summary |>
  left_join(nls_summary, by = "id") |>
  filter(!is.na(gi_hat), !is.na(gi_bayes))

write_rds(nls_vs_bayes, "objects/null-sim-nls-vs-bayes-validation.rds")

# Note: gi/gf show restricted-range correlations here because the *true*
# gi/gf are, by design of the null simulation, identical (fixed at the
# median) for every curve -- so both methods estimate values in a narrow
# band and small absolute estimation noise dominates the correlation
# coefficient despite close absolute agreement. tau and lambda vary across
# the real empirical range and are the parameters of primary interest for
# the null test, so their agreement is the most informative check.
validation_summary = nls_vs_bayes |>
  summarize(
    n = n(),
    cor_gi = cor(gi_bayes, gi_hat),
    cor_gf = cor(gf_bayes, gf_hat),
    cor_tau = cor(tau_bayes, tau_hat),
    cor_lambda = cor(lambda_bayes, lambda_hat),
    mean_abs_pct_diff_gi = mean(abs(gi_bayes - gi_hat) / gi_bayes) * 100,
    mean_abs_pct_diff_tau = mean(abs(tau_bayes - tau_hat) / tau_bayes) * 100,
    mean_abs_pct_diff_lambda = mean(abs(lambda_bayes - lambda_hat) / lambda_bayes) * 100
  )

write_rds(validation_summary, "objects/null-sim-nls-vs-bayes-summary.rds")

message(glue::glue(
  "nls vs. Bayesian validation ({validation_summary$n} curves): ",
  "cor(tau) = {formatC(validation_summary$cor_tau, format = 'f', digits = 3)}, ",
  "cor(lambda) = {formatC(validation_summary$cor_lambda, format = 'f', digits = 3)}, ",
  "mean abs %% diff tau = {formatC(validation_summary$mean_abs_pct_diff_tau, format = 'f', digits = 2)}%%."
))
