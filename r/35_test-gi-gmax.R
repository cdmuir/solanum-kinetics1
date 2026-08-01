# Reviewer comment R1.3 (latter part): fgmax = gi/gmax may not be an
# independent aperture measurement, since it is a ratio of a curve-fit
# parameter (gi) and an anatomically-derived quantity (gmax). To check
# whether fgmax carries information beyond its components, this script
# refits the selected tau/lambda model twice, replacing the logitfgmax fixed
# effect with (1) log(gi) alone and (2) log(gmax) alone, and compares model
# fit (LOOIC) and fixed-effect estimates against the original fgmax model.
#
# As in r/29_refit-vpd.R, brms does not support update()-ing the formula of
# a multivariate model, so each alternative is refit from scratch with
# brm() using the same sampler settings as r/10_fit-all.R, starting from the
# formula structure of the previously selected model (objects/best_model.rds).

source("r/header.R")

best_model = read_rds("objects/best_model.rds")

gi_gmax_data = read_rds("data/joined-summary.rds") |>
  prepare_tau_anatomy_data(logtau_threshold) |>
  select(accid, lighttreatment, lightintensity, leaftype, ginit_mean, gmax) |>
  mutate(loggi = log(ginit_mean), loggmax = log(gmax))

# --- Build data with gi/gmax joined on -----------------------------------

model_data = best_model$data |>
  left_join(gi_gmax_data,
            by = join_by(accid, lighttreatment, lightintensity, leaftype))

assert_true(all(!is.na(model_data$loggi)))
assert_true(all(!is.na(model_data$loggmax)))

# --- Build alternative formulas: replace logitfgmax with loggi or loggmax
# in the tau/lambda formulas only. gcl and fgmax formulas are unchanged
# (fgmax remains a jointly-modeled response either way).

bf_gcl = best_model$formula$forms$loggcl
bf_fgmax = best_model$formula$forms$logitfgmax

bf_lambda_gi = update(best_model$formula$forms$loglambdamean, . ~ . - logitfgmax + loggi)
bf_tau_gi = update(best_model$formula$forms$logtaumean, . ~ . - logitfgmax + loggi)
form_gi = bf_lambda_gi + bf_tau_gi + bf_gcl + bf_fgmax + set_rescor(TRUE)

bf_lambda_gmax = update(best_model$formula$forms$loglambdamean, . ~ . - logitfgmax + loggmax)
bf_tau_gmax = update(best_model$formula$forms$logtaumean, . ~ . - logitfgmax + loggmax)
form_gmax = bf_lambda_gmax + bf_tau_gmax + bf_gcl + bf_fgmax + set_rescor(TRUE)

# --- Fit both alternatives -------------------------------------------------

thin = 6

fit_gi = brm(
  formula = form_gi,
  data = model_data,
  data2 = list(A = best_model$data2$A),
  cores = 4,
  chains = 4,
  iter = thin * 2e3,
  thin = thin,
  refresh = thin * 1e2,
  control = list(adapt_delta = 0.9),
  backend = "cmdstanr",
  family = student(),
  seed = 8402911
) |> add_criterion("loo")

assert_true(check_convergence(fit_gi, convergence_criteria))
write_rds(fit_gi, "objects/best_model_gi.rds")

fit_gmax = brm(
  formula = form_gmax,
  data = model_data,
  data2 = list(A = best_model$data2$A),
  cores = 4,
  chains = 4,
  iter = thin * 2e3,
  thin = thin,
  refresh = thin * 1e2,
  control = list(adapt_delta = 0.9),
  backend = "cmdstanr",
  family = student(),
  seed = 8402912
) |> add_criterion("loo")

assert_true(check_convergence(fit_gmax, convergence_criteria))
write_rds(fit_gmax, "objects/best_model_gmax.rds")

# --- Compare LOOIC across the three predictor specifications --------------
#
# brms's loo_compare.brmsfit() takes brmsfit objects as separate arguments
# (not a list); model_names controls the row labels in the output.

loo_comparison = loo_compare(
  best_model, fit_gi, fit_gmax,
  criterion = "loo",
  model_names = c("fgmax", "gi", "gmax")
) |>
  as.data.frame() |>
  rownames_to_column("model")

write_rds(loo_comparison, "objects/tbl-comparison-gi-gmax.rds")

# --- Compare fixed-effect estimates for the predictor of interest ---------

extract_effect = function(fit, par_tau, par_lambda, predictor_name) {
  draws = as_draws_df(fit, variable = c(par_tau, par_lambda))
  bind_rows(
    tibble(
      predictor = predictor_name,
      response = "logtaumean",
      estimate = mean(draws[[par_tau]]),
      lower = quantile(draws[[par_tau]], 0.025),
      upper = quantile(draws[[par_tau]], 0.975)
    ),
    tibble(
      predictor = predictor_name,
      response = "loglambdamean",
      estimate = mean(draws[[par_lambda]]),
      lower = quantile(draws[[par_lambda]], 0.025),
      upper = quantile(draws[[par_lambda]], 0.975)
    )
  )
}

effect_comparison = bind_rows(
  extract_effect(best_model, "b_logtaumean_logitfgmax", "b_loglambdamean_logitfgmax", "fgmax"),
  extract_effect(fit_gi, "b_logtaumean_loggi", "b_loglambdamean_loggi", "gi"),
  extract_effect(fit_gmax, "b_logtaumean_loggmax", "b_loglambdamean_loggmax", "gmax")
)

write_rds(effect_comparison, "objects/tbl-effect-comparison-gi-gmax.rds")

message(glue::glue(
  "LOOIC comparison (higher elpd_diff = better fit, relative to best):\n",
  "{paste(capture.output(print(loo_comparison)), collapse = '\\n')}"
))
