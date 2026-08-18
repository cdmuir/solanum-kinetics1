# --- Refit the kinetics model with a curve-level VPD covariate ---------------
#
# Checks whether the gi -> tau effect is confounded by the realized VPD
# stimulus. brms does not support update()-ing the formula of a multivariate
# model (see ?update.brmsfit), so the model is refit from scratch with brm(),
# using the same sampler settings as r/10_fit-all.R and the formula structure
# of the previously selected model (objects/selected_model.rds).

source("r/header.R")

selected_model = read_rds("objects/selected_model.rds")
curve_vpd = read_rds("objects/curve-vpd-summary.rds") |>
  select(accid, lighttreatment, lightintensity, leaftype, final_vpd)

# --- Build data with VPD covariates joined on ---------------------------

joined_summary_amphi = selected_model$data |>
  left_join(curve_vpd,
            by = join_by(lighttreatment, lightintensity, leaftype, accid))

# --- Build updated formulas ---------------------------------------------

bf_lambda_vpd = update(selected_model$formula$forms$loglambdamean, . ~ . + final_vpd)

bf_tau_vpd = update(selected_model$formula$forms$logtaumean, . ~ . +
                      final_vpd)

bf_gcl_vpd = selected_model$formula$forms$loggcl

bf_gi_vpd = selected_model$formula$forms$loggi

bf_gmax_vpd = selected_model$formula$forms$loggmax

form_vpd = bf_lambda_vpd + bf_tau_vpd + bf_gcl_vpd + bf_gi_vpd + bf_gmax_vpd +
  set_rescor(TRUE)

# --- Fit ------------------------------------------------------------------

thin = selected_model$fit@sim$thin

fit_vpd = brm(
  formula = form_vpd,
  data = joined_summary_amphi |> mutate(phy = accession),
  data2 = list(A = selected_model$data2$A),
  cores = nchains(selected_model),
  chains = nchains(selected_model),
  iter = thin * 2e3,
  thin = thin,
  refresh = thin * 1e2,
  control = list(adapt_delta = 0.9),
  backend = "cmdstanr",
  family = student(),
  seed = 2716043
) |> add_criterion("loo")

write_rds(fit_vpd, "objects/selected_model_vpd.rds")

assert_true(check_convergence(fit_vpd, convergence_criteria))
