# Estimate Bayesian R2 for all curves
source("r/header.R")

sk_dir = "objects/sk-curves"

plan(multisession, workers = 10)

r2 = list.files(sk_dir) |>
  str_remove(".rds$") |>
  future_map_dfr(\(.x) {
    fit = read_rds(file.path(sk_dir, paste0(.x, ".rds")))
    fit |>
      bayes_R2() |>
      as_tibble() |>
      mutate(id = .x,
             converged = check_convergence(fit, convergence_criteria))
  }, .progress = TRUE)

write_rds(r2, "objects/r2.rds")
