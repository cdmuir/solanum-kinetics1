# Compare estimated to observed maximum gsw
source("r/header.R")

sk_dir = "objects/sk-curves"

plan(multisession, workers = 10)

df1 = list.files(sk_dir) |>
  str_remove(".rds$") |>
  future_map_dfr(\(.x) {
    fit = read_rds(file.path(sk_dir, paste0(.x, ".rds")))
    fit |>
      as_draws_df() |>
      mutate(gi = exp(b_loggf_Intercept) + exp(b_logdg_Intercept)) |>
      summarize_draws() |>
      filter(variable == "gi") |>
      mutate(id = .x,
             gsw_max = max(fit$data$gsw, na.rm = TRUE))
  }, .progress = TRUE)

ggplot(df1, aes(mean, gsw_max)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0, color = "red") +
  theme_minimal() +
  labs(x = "Estimated initial GSW",
       y = "Maximum observed GSW")
  