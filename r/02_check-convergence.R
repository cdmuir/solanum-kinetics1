source("r/header.R")

df1 = list.files("objects/sk-curves", full.names = TRUE) |>
  map(read_rds) |>
  map(as_draws_df) |>
  map(summarise_draws) |>
  map(as_tibble) |>
  map_dfr(filter, variable != "lprior") |>
  filter(rhat >= 1.05 | ess_bulk <= 400)

assert_true(nrow(df1) == 0)
