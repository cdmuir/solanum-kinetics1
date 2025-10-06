# Summarize parameter estimates for all curves
source("r/header.R")

sk_dir = "objects/sk-curves"

plan(multisession, workers = 10)

brms_summary = list.files(sk_dir) |>
  str_remove(".rds$") |>
  future_map_dfr(\(.x) {
    fit = read_rds(file.path(sk_dir, paste0(.x, ".rds")))
    fit |>
      as_draws_df() |>
      summarize_draws() |>
      mutate(id = .x)
  }, .progress = TRUE)

write_rds(brms_summary, "objects/brms-summary.rds")
