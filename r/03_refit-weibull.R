# Refit each CDWeibull model curve to achieve convergence if needed
source("r/header.R")

plan(multisession, workers = 15)

sk_dir = c("objects/weibull/")
list.files(sk_dir, full.names = TRUE) |>
  future_walk(
    refit_rh,
    convergence_criteria = convergence_criteria,
    .progress = TRUE,
    .options = furrr_options(seed = TRUE)
  )

# zip
zipr("objects/weibull.zip",
     list.files(sk_dir, full.names = TRUE),
     recurse = TRUE)
