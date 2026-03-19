# Fit CDWeibulll model to each curve separately
source("r/header.R")

sk_dir = c("objects/weibull/")
if (!dir.exists(sk_dir)) {
  dir.create(sk_dir)
}

joined_data = read_rds("data/joined-data.rds")

plan(multisession, workers = 19)

# Alternative formulation of CDWeibull
form_cdweibull = gsw ~ exp(loggf) + exp(logdg) * exp(-(t_sec / exp(logtau))^exp(loglambda))

bform_cdweibull = bf(form_cdweibull, loggf ~ 1, logdg ~ 1, logtau ~ 1, loglambda ~ 1, nl = TRUE)

pri = c(
  prior(
    normal(log(0.1), 2),
    nlpar = "loggf",
    lb = log(0.001),
    ub = log(2)
  ),
  prior(
    normal(log(0.5), 2),
    nlpar = "logdg",
    lb = log(0.001),
    ub = log(2)
  ),
  prior(
    normal(log(300), 1),
    nlpar = "logtau",
    lb = log(10),
    ub = log(10000)
  ),
  prior(
    normal(log(2), 1),
    nlpar = "loglambda",
    lb = log(0.1),
    ub = log(10)
  )
)

joined_data |>
  split( ~ curve) |>
  future_iwalk(\(df, curve_id) {
    
    file = paste0(sk_dir, curve_id, ".rds")

    fit = fit_rh1(
      formula = bform_cdweibull,
      data = df,
      prior = pri,
      thin = 2,
      adapt_delta = 0.8,
      seed = 360036340 + round(1000 * df$gsw[1])
    )
    
    write_rds(fit, file)

  }, .progress = TRUE, .options = furrr_options(seed = TRUE))

  