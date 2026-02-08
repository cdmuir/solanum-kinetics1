# Estimate phylogenetic h2 from fitted models
source("r/header.R")

fit_amphi_high = read_rds("objects/fit_amphi_high.rds")
fit_amphi_low = read_rds("objects/fit_amphi_low.rds")

bind_rows(
  get_phy_h2(fit_amphi_high) |>
    mutate(light_intensity = "high"),
  get_phy_h2(fit_amphi_low) |>
    mutate(light_intensity = "low")
) |>
  mutate(
    Trait = case_when(
      resp == "loglambdamean" ~ "log(lambda)",
      resp == "logtaumean" ~ "log(tau)",
      resp == "loggcl" ~ "log(guard~cell~length)",
      resp == "logfgmax" ~ "log(italic(f)[gmax])"
    )
  ) |>
  ggplot(aes(
    x = Trait,
    y = estimate,
    ymin = lowerCI,
    ymax = upperCI
  )) +
  facet_wrap( ~ light_intensity) +
  geom_pointrange() +
  scale_x_discrete(labels = scales::label_parse()) +
  labs(y = "Phylogenetic heritability") +
  ylim(0, 1) +
  theme(axis.text.x = element_text(angle = 45, h  = 1))

