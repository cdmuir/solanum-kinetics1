# Fit model 3 to each curve separately
source("r/header.R")

sk_dir = c("objects/sk-curves3/")
if (!dir.exists(sk_dir)) {
  dir.create(sk_dir)
}

rh_curves = read_rds("data/rh_curves.rds")

plan(multisession, workers = 9)
# a_exact = a_max + (a_star - a_max) * R ^ exp(-t/tau_P)

form3 = gsw ~ exp(loggmax) + (exp(loggf) - exp(loggmax)) * R ^ exp(-(t_sec / exp(logtau))^exp(loglambda))

bform3 = bf(form3, loggf ~ 1, R ~ 1, loggmax ~ 1, logtau ~ 1, loglambda ~ 1, nl = TRUE)

pri = c(
  prior(
    normal(log(0.1), 2),
    nlpar = "loggf",
    lb = log(0.001),
    ub = log(2)
  ),
  prior(
    normal(0, 1),
    nlpar = "R"
  ),
  prior(
    normal(log(300), 1),
    nlpar = "logtau",
    lb = log(10),
    ub = log(10000)
  ),
  prior(
    normal(log(0.5), 2),
    nlpar = "loggmax",
    lb = log(0.001),
    ub = log(2)
  ),
  prior(
    normal(log(2), 1),
    nlpar = "loglambda",
    lb = log(0.1),
    ub = log(10)
  )
)

rh_curves |>
  split(~ curve) |>
  magrittr::extract(1:9) |>
  future_iwalk(\(df, curve_id) {
    file = paste0(sk_dir, curve_id, ".rds")
    
    fit = fit_rh1(
      formula = bform3,
      data = df,
      prior = pri,
      thin = 2,
      adapt_delta = 0.8,
      seed = 360036340 + df$ci[1]
    )
    
    write_rds(fit, file)
    
  }, .progress = TRUE, .options = furrr_options(seed = TRUE))

