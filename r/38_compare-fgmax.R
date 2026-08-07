# Compare best_model to best_model_gi and best_model_gmax
source("r/header.R")

best_model = read_rds("objects/best_model.rds")
best_model_gi = read_rds("objects/best_model_gi.rds")
best_model_gmax = read_rds("objects/best_model_gmax.rds")

# --- Compare LOOIC across the three predictor specifications --------------
#
# brms's loo_compare.brmsfit() takes brmsfit objects as separate arguments
# (not a list); model_names controls the row labels in the output.

loo_comparison = loo_compare(
  best_model, best_model_gi, best_model_gmax,
  criterion = "loo",
  model_names = c("fgmax", "gi", "gmax")
)

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
  extract_effect(best_model_gi, "b_logtaumean_loggi", "b_loglambdamean_loggi", "gi"),
  extract_effect(best_model_gmax, "b_logtaumean_loggmax", "b_loglambdamean_loggmax", "gmax")
)

write_rds(effect_comparison, "objects/tbl-effect-comparison-gi-gmax.rds")
