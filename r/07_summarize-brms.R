# Summarize parameter estimates for all curves
source("r/header.R")

assert_true(setequal(
  list.files("objects/weibull"),
  list.files("objects/inertia")
))

sk_dir1 = "objects/weibull"
sk_dir2 = "objects/inertia"

plan(multisession, workers = 19)

brms_summary = list.files(sk_dir1) |>
  future_map_dfr(\(.x) {
    fit_weibull = read_rds(file.path(sk_dir1, .x))
    fit_inertia = read_rds(file.path(sk_dir2, .x))
    
    fit_weibull |>
      as_draws_df() |>
      mutate(log_gi = log(exp(b_loggf_Intercept) + exp(b_logdg_Intercept))) |>
      summarize_draws() |>
      mutate(id = str_remove(.x, "\\.rds$"), model = "weibull")
    fit_inertia |>
      as_draws_df() |>
      summarize_draws() |>
      mutate(id = str_remove(.x, "\\.rds$"), model = "inertia")
  }, .progress = TRUE)

write_rds(brms_summary, "objects/brms-summary.rds")
