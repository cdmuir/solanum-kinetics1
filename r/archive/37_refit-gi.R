# Refit best model with gi as a covariate in place of fgmax
source("r/header.R")

best_model = read_rds("objects/best_model.rds")

gi_data = read_rds("data/joined-summary.rds") |>
  prepare_tau_anatomy_data(logtau_threshold) |>
  select(accid, lighttreatment, lightintensity, leaftype, ginit_mean) |>
  mutate(loggi = log(ginit_mean))

# --- Build data with gi/gmax joined on -----------------------------------

model_data = best_model$data |>
  left_join(gi_data,
            by = join_by(accid, lighttreatment, lightintensity, leaftype))

assert_true(all(!is.na(model_data$loggi)))

# --- Build alternative formulas: replace logitfgmax with loggi in the 
# tau/lambda formulas only. gcl and fgmax formulas are unchanged (fgmax remains 
# a jointly-modeled response either way).

bf_gcl = best_model$formula$forms$loggcl
bf_fgmax = best_model$formula$forms$logitfgmax

bf_lambda_gi = update(best_model$formula$forms$loglambdamean, . ~ . - logitfgmax + loggi)
bf_tau_gi = update(best_model$formula$forms$logtaumean, . ~ . - logitfgmax + loggi)
form_gi = bf_lambda_gi + bf_tau_gi + bf_gcl + bf_fgmax + set_rescor(TRUE)

# --- Fit model ----------------------------------------------------------

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

# assert_true(check_convergence(fit_gi, convergence_criteria))

write_rds(fit_gi, "objects/best_model_gi.rds")
