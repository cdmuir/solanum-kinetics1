# Summarize parameter estimates for all curves
source("r/header.R")

sk_dir1 = "objects/weibull"

plan(multisession, workers = 19)

pars_summary = list.files(sk_dir1) |>
  future_map_dfr(\(.x) {
    fit_weibull = read_rds(file.path(sk_dir1, .x))
    
    fit_weibull |>
      as_draws_df() |>
      mutate(ginit =
               exp(b_loggf_Intercept) + exp(b_logdg_Intercept)) |>
      summarize_draws() |>
      mutate(id = str_remove(.x, "\\.rds$"))
    
  }, .progress = TRUE)

write_rds(pars_summary, "objects/pars-summary.rds")
