# Refit each curve to achieve convergence
source("r/header.R")

plan(multisession, workers = 16)

sk_dir = c("objects/sk-curves/")
list.files(sk_dir, full.names = TRUE) |>
  future_walk(refit_rh, convergence_criteria = convergence_criteria)
