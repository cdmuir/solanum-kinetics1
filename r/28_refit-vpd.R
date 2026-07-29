# Reviewer comment R1.1: refit the tau/lambda model with curve-level VPD
# covariates (final VPD and VPD saturation rate, with measurement error) to
# check whether the fgmax -> tau effect survives accounting for realized VPD
# exposure.
#
# Restricted to amphistomatous (amphi) curves only, because the VPD
# covariates come from data/rh_curves.rds and we want a single, unambiguous
# leaf-type context for interpreting them (leaftype cannot be included as a
# fixed effect once the data are restricted to one leaf type).
#
# final_vpd and me(sat_rate, sat_rate_se) are each given a separate
# coefficient within every combination of lighttreatment x lightintensity
# (4 groups), because the realized VPD trajectory -- and its relationship to
# tau/lambda -- may differ by treatment combination. This only changes the
# tau and lambda formulas; the gcl and fgmax formulas are unchanged apart
# from dropping leaftype (which is otherwise a constant in the amphi-only
# subset and cannot be estimated).
#
# brms does not support update()-ing the formula of a multivariate model
# (see ?update.brmsfit), so the model is refit from scratch with brm() using
# the same sampler settings as r/10_fit-all.R, starting from the formula
# structure of the previously selected model (objects/best_model.rds).

source("r/header.R")

best_model = read_rds("objects/best_model.rds")
curve_vpd = read_rds("objects/curve-vpd-summary.rds") |>
  mutate(
    lighttreatment = case_when(
      light_treatment == "low" ~ "shade",
      light_treatment == "high" ~ "sun"
    ),
    lightintensity = case_when(
      light_intensity == "150" ~ "low",
      light_intensity == "2000" ~ "high"
    ),
    leaftype = case_when(
      curve_type == "2-sided RH" ~ "amphi",
      curve_type == "1-sided RH" ~ "pseudohypo"
    )
  ) |>
  rename(
    accid = acc_id
  ) |>
  select(
    accid, lighttreatment, lightintensity, leaftype, final_vpd, sat_rate, sat_rate_se
  )

# --- Build amphi-only data with VPD covariates joined on ---------------

joined_summary_amphi = best_model$data |>
  filter(leaftype == "amphi") |>
  left_join(curve_vpd,
            by = join_by(lighttreatment, lightintensity, leaftype, accid))

# --- Build updated formulas ---------------------------------------------
# Add final_vpd and me(sat_rate, sat_rate_se), each with a separate slope
# per lighttreatment:lightintensity combination, to the tau and lambda
# formulas only. Drop leaftype everywhere (amphi-only data).

bf_lambda_vpd = update(best_model$formula$forms$loglambdamean,
                       . ~ . - leaftype +
                         finalvpd)

bf_tau_vpd = update(best_model$formula$forms$logtaumean, . ~ . - leaftype +
                      finalvpd)

bf_gcl_vpd = update(best_model$formula$forms$loggcl, . ~ . - leaftype)

bf_fgmax_vpd = update(best_model$formula$forms$logitfgmax, . ~ . - leaftype)

form_vpd = bf_lambda_vpd + bf_tau_vpd + bf_gcl_vpd + bf_fgmax_vpd +
  set_rescor(TRUE)

# --- Fit ------------------------------------------------------------------

thin = 6

fit_vpd = brm(
  formula = form_vpd,
  data = joined_summary_amphi |> mutate(phy = accession),
  data2 = list(A = best_model$data2$A),
  cores = 4,
  chains = 4,
  iter = thin * 2e3,
  thin = thin,
  refresh = thin * 1e2,
  control = list(adapt_delta = 0.9),
  backend = "cmdstanr",
  family = student(),
  seed = 2716043
) |> add_criterion("loo")

write_rds(fit_vpd, "objects/best_model_vpd.rds")

assert_true(check_convergence(fit_vpd, convergence_criteria))
