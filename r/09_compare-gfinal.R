# Compare minimum gsw to tau within treatments
# See if any confounding between RH_cham and tau
source("r/header.R")

joined_summary = read_rds("data/joined-summary.rds") |>
  filter(logtau_mean < logtau_threshold)

joined_summary  |>
  ggplot(aes(gfinal_mean, logtau_mean)) +
  geom_point() +
  facet_grid(light_intensity ~ light_treatment + curve_type) +
  scale_x_log10()

joined_summary |>
  ggplot(aes(ginit_mean, logtau_mean)) +
  geom_point() +
  facet_grid(light_intensity ~ light_treatment + curve_type, scales = "free") +
  scale_x_log10()

joined_summary |>
  ggplot(aes(ginit_mean, loglambda_mean)) +
  geom_point() +
  facet_grid(light_intensity ~ light_treatment + curve_type, scales = "free") +
  scale_x_log10()
