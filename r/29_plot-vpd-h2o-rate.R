# --- Check that the VPD saturation rate is driven by the sample, not reference, airstream ---
#
# The same saturating model fit to VPD (r/27_summarize-vpd.R) is also fit to
# the two raw humidity channels it is built from: H2O_s (sample cell, i.e. the
# air actually exchanged with the leaf) and H2O_r (reference cell, i.e.
# incoming air upstream of the leaf). H2O_s should behave like VPD (both
# reflect declining leaf transpiration as stomata close), while H2O_r should be
# comparatively flat (it mostly reflects the controlled/scrubbed incoming air,
# not the leaf response). This is a check, not a covariate source -- see bottom
# of script.

source("r/header.R")

rh_curves = read_rds("data/rh_curves.rds")
curve_vpd = read_rds("objects/curve-vpd-summary.rds")

# --- Raw H2O_r and H2O_s trajectories -------

h2os_saturating = rh_curves |>
  group_by(curve) |>
  group_modify( ~ fit_saturating(.x, "H2O_s")) |>
  ungroup()

h2or_saturating = rh_curves |>
  group_by(curve) |>
  group_modify( ~ fit_saturating(.x, "H2O_r")) |>
  ungroup()

n_failed_h2os = sum(!h2os_saturating$converged)
n_failed_h2or = sum(!h2or_saturating$converged)
if (n_failed_h2os > 0 || n_failed_h2or > 0) {
  warning(
    glue::glue(
      "{n_failed_h2os} curve(s) failed to converge for H2O_s and ",
      "{n_failed_h2or} curve(s) failed to converge for H2O_r. ",
      "Inspect objects/curve-vpd-sensitivity.rds."
    )
  )
}

curve_vpd_sensitivity = curve_vpd |>
  select(curve,
         leaf_type = curve_type,
         vpd_rate = sat_rate,
         vpd_r2 = r2_sat) |>
  left_join(
    h2os_saturating |>
      select(
        curve,
        h2os_rate = sat_rate,
        h2os_r2 = r2_sat,
        h2os_asym = sat_asym,
        h2os_r0 = sat_r0
      ),
    by = "curve"
  ) |>
  left_join(
    h2or_saturating |>
      select(
        curve,
        h2or_rate = sat_rate,
        h2or_r2 = r2_sat,
        h2or_asym = sat_asym,
        h2or_r0 = sat_r0
      ),
    by = "curve"
  ) |>
  mutate(
    # amplitude of change over the curve (mmol/mol); direction is implicit
    # in the sign of (asym - r0)
    h2os_amplitude = abs(h2os_asym - h2os_r0),
    h2or_amplitude = abs(h2or_asym - h2or_r0)
  )

write_rds(curve_vpd_sensitivity, "objects/curve-vpd-sensitivity.rds")

sensitivity_cor = cor(curve_vpd_sensitivity |> select(vpd_rate, h2os_rate, h2or_rate),
                      use = "pairwise.complete.obs")

p_h2os_vs_vpd = ggplot(curve_vpd_sensitivity, aes(vpd_rate, h2os_rate)) +
  geom_point(alpha = 0.4) +
  geom_abline(
    slope = 1,
    intercept = 0,
    linetype = "dashed",
    color = "grey50"
  ) +
  labs(x = "VPD saturation rate (1/s)", y = expression(italic(W)[s] ~ "saturation rate (1/s)"))

p_h2or_vs_vpd = ggplot(curve_vpd_sensitivity, aes(vpd_rate, h2or_rate)) +
  geom_point(alpha = 0.4) +
  labs(x = "VPD saturation rate (1/s)", y = expression(italic(W)[r] ~ "saturation rate (1/s)"))

p_vpd_sensitivity = (p_h2os_vs_vpd + p_h2or_vs_vpd) +
  plot_annotation(tag_levels = "a")

ggsave("figures/vpd-h2o-rate.pdf",
       width = 8,
       height = 4)
