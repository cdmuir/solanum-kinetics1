# Compare estimated to observed maximum gsw
source("r/header.R")

sk_dir1 = "objects/weibull"

plan(multisession, workers = 19)

df1 = list.files(sk_dir1) |>
  future_map_dfr(\(.x) {
    fit_weibull = read_rds(file.path(sk_dir1, .x))
    
    fit_weibull |>
      as_draws_df() |>
      mutate(gi = exp(b_loggf_Intercept) + exp(b_logdg_Intercept)) |>
      summarize_draws() |>
      filter(variable == "gi") |>
      mutate(id = .x,
             gsw_max = max(fit_weibull$data$gsw, na.rm = TRUE))
    
  }, .progress = TRUE)

ggplot(df1, aes(mean, gsw_max)) +
  geom_point() +
  geom_abline(slope = 1,
              intercept = 0,
              color = "tomato") +
  scale_x_log10() +
  scale_y_log10() +
  theme_cowplot() +
  coord_equal() +
  labs(x = expression(Maximum ~ italic(g)[sw] ~ (estimated)),
       y = expression(Maximum ~ italic(g)[sw] ~ (observed)))

ggsave("figures/compare-gsw.pdf",
       width = 4,
       height = 4)
