# Fit model 3 (power function) to each curve separately
# this is not quite correct since it's using raw g, not g normalized distance to gmax
source("r/header.R")

sk_dir = c("objects/sk-curves3/")
if (!dir.exists(sk_dir)) {
  dir.create(sk_dir)
}

rh_curves = read_rds("data/rh_curves.rds")

plan(multisession, workers = 9)
# g_final ^ (1/k) + (g_init ^ (1/k) - g_final ^ (1/k)) * exp(-t/tau_P)
# let lambda = -1/k

form3 = gsw ~ exp(loggf) ^ -exp(loglambda) + (exp(loggi) ^ -exp(loglambda) - exp(loggf) ^ -exp(loglambda)) * exp(-t_sec / exp(logtau)) 

bform3 = bf(form3, loggf ~ 1, loggi ~ 1, logtau ~ 1, loglambda ~ 1, nl = TRUE)

pri = c(
  prior(
    normal(log(0.1), 2),
    nlpar = "loggf",
    lb = log(0.001),
    ub = log(2)
  ),
  prior(
    normal(log(0.1), 2),
    nlpar = "loggi",
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

