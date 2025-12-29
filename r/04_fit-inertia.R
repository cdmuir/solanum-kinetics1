# Fit inertia model to each curve separately
source("r/header.R")

sk_dir = c("objects/inertia/")
if (!dir.exists(sk_dir)) {
  dir.create(sk_dir)
}

joined_data = read_rds("data/joined-data.rds")

plan(multisession, workers = 19)
# gmin = minimum stomatal conductance
# gmax = anatomical maximum conductance
# ginit = initial conductance
# gstar = target conductance

form_inertia = gsw ~ gmax + (gmin - gmax) * (((gstar - gmax) / (gmin - gmax)) ^ ik + (((ginit - gmax) / (gmin - gmax)) ^ ik - ((gstar - gmax) / (gmin - gmax)) ^ ik) * exp(-t_sec / exp(logtau))) ^ (1 / ik) 

bform_inertia = bf(form_inertia, gmin ~ 1, gstar ~ 1, ginit ~ 1, ik ~ 1, logtau ~ 1, nl = TRUE)

joined_data |>
  split(~ curve) |>
  magrittr::extract(1:9) |>
  future_iwalk(\(df, curve_id) {
    file = paste0(sk_dir, curve_id, ".rds")
    
    fit = fit_rh1(
      formula = bform_inertia,
      data = df,
      prior = get_interia_prior(df),
      thin = 2,
      adapt_delta = 0.8,
      seed = 360036340 + round(1e3 * df$gsw[1])
    )
    
    write_rds(fit, file)
    
  }, .progress = TRUE, .options = furrr_options(seed = TRUE))

