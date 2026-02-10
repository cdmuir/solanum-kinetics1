# Estimate phylogenetic h2 from fitted models
source("r/header.R")

fit_amphi = read_rds("objects/best_amphi_model.rds")

get_phy_h2(fit_amphi)

get_phy_h2(fit_amphi) |>
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
  geom_pointrange() +
  scale_x_discrete(labels = scales::label_parse()) +
  labs(y = "Phylogenetic heritability") +
  ylim(0, 1) +
  theme(axis.text.x = element_text(angle = 45, h  = 1))

