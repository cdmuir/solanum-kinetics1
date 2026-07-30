# --- Relationship between tau and VPD saturation rate ---------------------
#
# The VPD saturation rate k is strongly correlated with the fitted kinetics
# time-constant tau, which is why k is not used as a covariate for tau (that
# would be close to circular). This script quantifies and plots that
# relationship directly.

source("r/header.R")

best_model_vpd = read_rds("objects/best_model_vpd.rds")
curve_vpd = read_rds("objects/curve-vpd-summary.rds") |>
  filter(!exclude_sat_rate) |>
  select(accid, lighttreatment, lightintensity, leaftype, sat_rate)

tau_vpd_rate = best_model_vpd$data |>
  select(lighttreatment, lightintensity, accid, logtaumean) |>
  left_join(curve_vpd, by = join_by(lighttreatment, lightintensity, accid))

p_tau_vs_rate = ggplot(tau_vpd_rate, aes(1/sat_rate, exp(logtaumean))) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "lm", se = FALSE) +
  labs(x = expression(1 / VPD~saturation~rate * "," ~italic(k^-1)~"(s, log scale)"),
       y = expression(tau~"(s, log scale)")) +
  scale_x_log10() +
  scale_y_log10() +
  coord_equal()

ggsave("figures/tau-vpd-rate.pdf", width = 5, height = 4)
