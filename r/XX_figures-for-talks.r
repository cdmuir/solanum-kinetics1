source("r/header.R")

fit_amphi_high = read_rds("objects/fit_amphi_high.rds")
fit_amphi_low = read_rds("objects/fit_amphi_low.rds")

df_high = fit_amphi_high$data

df_high |>
  summarize(loggcl = mean(loggcl), 
            logtau = mean(logtaumean),
            .by = c(phy, lighttreatment)) |>
  rename(`Light treatment` = lighttreatment) |>
  ggplot(aes(x = exp(loggcl), y = exp(logtau), color = `Light treatment`)) +
  scale_color_manual(values = c("tomato", "tomato4")) +
  geom_point(size = 3, alpha = 0.7) +
  scale_x_log10() +
  scale_y_log10() +
  xlab("Guard cell length (um)") +
  ylab("Stomatal closing time constant (seconds)")
ggsave("figures/amphi_high_gcl_vs_tau.png", width = 6, height = 4)


df_high |>
  rename(`Light treatment` = lighttreatment) |>
  ggplot(aes(x = exp(logfgmax) * 1e3, y = exp(logtaumean), color = `Light treatment`)) +
  scale_color_manual(values = c("tomato", "tomato4")) +
  geom_point(size = 3, alpha = 0.7) +
  scale_x_log10() +
  scale_y_log10() +
  xlab("Fraction maximum opening") +
  ylab("Stomatal closing time constant (seconds)")
ggsave("figures/amphi_high_fgmax_vs_tau.png", width = 6, height = 4)
