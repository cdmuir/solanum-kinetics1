# --- Validate the fast nls() refit against the full Bayesian pipeline --------
#
# Confirms that the fast nls() re-fit used in the null simulation
# (r/33_simulate-null.R) is a reasonable approximation to the production
# Bayesian pipeline (fit_rh1()/bform_cdweibull/pri from r/02_fit-weibull.R and
# r/functions.R), by refitting a subsample of simulated curves both ways and
# comparing.

source("r/header.R")

joined_data = read_rds("data/joined-data.rds") |>
  rename(gsw_sim = gsw, accession = acc) |>
  mutate(leaftype = recode_leaftype(curve_type),
         lightintensity = recode_lightintensity(str_extract(curve, "150|2000$")))

brms_estimates = read_csv("tables/tbl-estimates-curve.csv", show_col_types = FALSE) |>
  select(
    accession,
    id,
    leaftype,
    lightintensity,
    loglambda_brms = loglambdamean,
    logtau_brms = logtaumean
  ) |>
  mutate(
    leaftype = factor(leaftype, levels = levels(joined_data$leaftype)),
    lightintensity = factor(lightintensity, levels = levels(joined_data$lightintensity))
  )

nls_estimates = split(joined_data, . ~ curve) |>
  map_dfr(\(.x) {
    bind_cols(summarize(.x, across(
      c(accession, id, leaftype, lightintensity), first
    )), fit_nls_one(.x))
  }) |>
  mutate(loglambda_nls = log(lambda_hat),
         logtau_nls = log(tau_hat)) |>
  select(accession,
         id,
         leaftype,
         lightintensity,
         loglambda_nls,
         logtau_nls)

nls_vs_bayes = left_join(brms_estimates,
                         nls_estimates,
                         by = join_by(accession, id, leaftype, lightintensity))

write_rds(nls_vs_bayes, "objects/nls_vs_bayes.rds")
