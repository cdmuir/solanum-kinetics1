# Compare estimated to observed maximum gsw
source("r/header.R")

assert_true(setequal(
  list.files("objects/weibull"),
  list.files("objects/inertia")
))

sk_dir1 = "objects/weibull"
sk_dir2 = "objects/inertia"

plan(multisession, workers = 19)

df1 = list.files(sk_dir1) |>
  future_map_dfr(\(.x) {
    fit_weibull = read_rds(file.path(sk_dir1, .x))
    fit_inertia = read_rds(file.path(sk_dir2, .x))
    
    bind_rows(
      fit_weibull |>
        as_draws_df() |>
        mutate(gi = exp(b_loggf_Intercept) + exp(b_logdg_Intercept)) |>
        summarize_draws() |>
        filter(variable == "gi") |>
        mutate(model = "weibull"),
      
      fit_inertia |>
        as_draws_df() |>
        summarize_draws() |>
        filter(variable == "b_ginit_Intercept") |>
        mutate(model = "intertia")
    ) |>
      mutate(id = .x,
             gsw_max = max(fit_weibull$data$gsw, na.rm = TRUE))
    
  }, .progress = TRUE)

ggplot(df1, aes(mean, gsw_max)) +
  facet_grid(. ~ model) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0, color = "red") +
  theme_minimal() +
  labs(x = "Estimated initial GSW",
       y = "Maximum observed GSW")
  