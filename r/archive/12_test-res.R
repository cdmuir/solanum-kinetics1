# GAVE UP - too complicated
# Compare models to test which random effects should be included in the final
# model
source("r/header.R")

plan(multisession, workers = 9)

joined_summary = read_rds("data/joined-summary.rds") |>
  prepare_tau_anatomy_data(logtau_threshold)

phy = read_rds("data/phylogeny.rds")
A = vcv(phy, corr = TRUE)
thin = 1

nesting(
  response = c(
    "loglambdamean | se(loglambdasd, sigma = TRUE)",
    "logtaumean | se(logtausd, sigma = TRUE)",
    "logitfgmax",
    "loggcl"
  ),
  fixed_effects = rep(
    c(
      "lighttreatment + lightintensity + curve_type",
      "lighttreatment + curve_type"
    ),
    c(3, 1)
  )
)

crossing(
  response = "loggcl",
  fixed_effects = "lighttreatment + curve_type",
 
  crossing(
   crossing(
    group = "accid",
    form = c("1", "1 + lighttreatment", "1 + curve_type", "1 + lighttreatment + curve_type")
  ) |>
    unite("accid_form", group, form, sep = "_"),
  
  crossing(
    group = "accession",
    form = c("1", "1 + lighttreatment", "1 + lightintensity", "1 + curve_type", "1 + lighttreatment + lightintensity", "1 + lighttreatment + curve_type", "1 + lightintensity + curve_type", "1 + lighttreatment + lightintensity + curve_type")
  ) |>
    unite("accession_form", group, form, sep = "_")
  )
  
)
  

fixef = 
# Define formula
bf_lambda0 = bf(loglambdamean | se(loglambdasd, sigma = TRUE) ~ 
                  lighttreatment + 
                  lightintensity +
                  curve_type +
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
               curve_type +
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
              curve_type +
              (1|a|accession) + 
              (1|b|gr(phy, cov = A)))

bf_fgmax = bf(logitfgmax ~ 
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
