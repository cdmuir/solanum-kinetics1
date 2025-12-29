# Estimate Bayesian R2 for all curves
source("r/header.R")

assert_true(setequal(
  list.files("objects/weibull"),
  list.files("objects/inertia")
))

sk_dir1 = "objects/weibull"
sk_dir2 = "objects/inertia"

plan(multisession, workers = 19)

r2 = list.files(sk_dir1) |>
  future_map_dfr(\(.x) {
    fit_weibull = read_rds(file.path(sk_dir1, .x))
    fit_inertia = read_rds(file.path(sk_dir2, .x))
    
    bind_rows(
      fit_weibull |>
        bayes_R2() |>
        as_tibble() |>
        mutate(
          id = str_remove(.x, "\\.rds$"),
          converged = check_convergence(fit_weibull, convergence_criteria),
          model = "weibull"
        ),
      fit_inertia |>
        bayes_R2() |>
        as_tibble() |>
        mutate(
          id = str_remove(.x, "\\.rds$"),
          converged = check_convergence(fit_inertia, convergence_criteria),
          model = "inertia"
        )
    )
    
  }, .progress = TRUE)

write_rds(r2, "objects/r2.rds")
