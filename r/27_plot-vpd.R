# WORKING HERE - breakin up into multiple smaller scripts and revising contect from PostAI

# The same saturating model is also fit to the two raw humidity channels the VPD 
# calculation is built from: H2O_s (sample cell, i.e. the air actually exchanged 
# with the leaf) and H2O_r (reference cell, i.e. incoming air upstream of the 
# leaf). We expect H2O_s to behave like VPD (both reflect leaf transpiration 
# declining as stomata close), and H2O_r to be comparatively flat (it should 
# mostly reflect the controlled/scrubbed incoming air, not the leaf response). 
# This is a check, not a covariate source -- see bottom of script.

source("r/header.R")

curve_vpd = read_rds("objects/curve-vpd-summary.rds")
pairs_data = curve_vpd |>
  select(median_vpd, final_vpd, vpd_rate = sat_rate)

# Pairwise correlations (Pearson)
vpd_cor = cor(pairs_data, use = "pairwise.complete.obs")
print(vpd_cor)

# Scatterplot matrix, colored by leaf type since Table S3 highlighted
# leaf-type differences in realized VPD
plot_data = curve_vpd |>
  mutate(
    leaf_type = factor(
      case_when(
        curve_type == "1-sided RH" ~ "pseudohypo",
        curve_type == "2-sided RH" ~ "amphi"
      ),
      levels = c("amphi", "pseudohypo")
    )
  )

p_median_final = ggplot(plot_data, aes(median_vpd, final_vpd, color = leaf_type)) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "lm", se = FALSE) +
  scale_color_manual(values = c(amphi = col_amphi, pseudohypo = col_pseudohypo)) +
  labs(x = "Median VPD (kPa)", y = "Final VPD (kPa)", color = "Leaf type")

p_median_rate = ggplot(plot_data, aes(median_vpd, sat_rate, color = leaf_type)) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "lm", se = FALSE) +
  scale_color_manual(values = c(amphi = col_amphi, pseudohypo = col_pseudohypo)) +
  labs(x = "Median VPD (kPa)", y = "VPD saturation rate (1/s)", color = "Leaf type")

p_final_rate = ggplot(plot_data, aes(final_vpd, sat_rate, color = leaf_type)) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "lm", se = FALSE) +
  scale_color_manual(values = c(amphi = col_amphi, pseudohypo = col_pseudohypo)) +
  labs(x = "Final VPD (kPa)", y = "VPD saturation rate (1/s)", color = "Leaf type")

p_vpd_pairs = (p_median_final + p_median_rate + p_final_rate) +
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = "A") &
  theme(legend.position = "bottom")

print(p_vpd_pairs)

write_rds(p_vpd_pairs, "objects/plot-vpd-pairs.rds")

# Diagnostic: distribution of saturating-fit R2 across curves, to document
# fit quality alongside the reviewer response.
p_r2_sat = ggplot(curve_vpd, aes(r2_sat)) +
  geom_histogram(bins = 30) +
  labs(
    x = expression(R^2 ~ "of saturating VPD fit"),
    y = "Number of curves"
  )

print(p_r2_sat)

write_rds(p_r2_sat, "objects/plot-vpd-sat-r2.rds")

# --- Sensitivity/validation check: raw H2O_r and H2O_s trajectories -------
#
# VPDleaf is calculated from H2O_s (see `vpd_traj` above); this section
# checks that (a) the VPD saturation rate is not an artifact of the
# saturating-fit procedure by confirming it closely tracks the saturation
# rate fit directly to H2O_s, and (b) H2O_r (reference/incoming air, largely
# outside the leaf's influence) is comparatively flat, as expected if it is
# not part of the leaf-driven signal.

h2os_saturating = rh_curves |>
  group_by(curve) |>
  group_modify(~ fit_saturating(.x, "H2O_s")) |>
  ungroup()

h2or_saturating = rh_curves |>
  group_by(curve) |>
  group_modify(~ fit_saturating(.x, "H2O_r")) |>
  ungroup()

n_failed_h2os = sum(!h2os_saturating$converged)
n_failed_h2or = sum(!h2or_saturating$converged)
if (n_failed_h2os > 0 || n_failed_h2or > 0) {
  warning(glue::glue(
    "{n_failed_h2os} curve(s) failed to converge for H2O_s and ",
    "{n_failed_h2or} curve(s) failed to converge for H2O_r. ",
    "Inspect objects/curve-vpd-sensitivity.rds."
  ))
}

curve_vpd_sensitivity = curve_vpd |>
  select(curve, leaf_type = curve_type, vpd_rate = sat_rate, vpd_r2 = r2_sat) |>
  left_join(
    h2os_saturating |>
      select(curve, h2os_rate = sat_rate, h2os_r2 = r2_sat,
             h2os_asym = sat_asym, h2os_r0 = sat_r0),
    by = "curve"
  ) |>
  left_join(
    h2or_saturating |>
      select(curve, h2or_rate = sat_rate, h2or_r2 = r2_sat,
             h2or_asym = sat_asym, h2or_r0 = sat_r0),
    by = "curve"
  ) |>
  mutate(
    # amplitude of change over the curve (mmol/mol); direction is implicit
    # in the sign of (asym - r0)
    h2os_amplitude = abs(h2os_asym - h2os_r0),
    h2or_amplitude = abs(h2or_asym - h2or_r0)
  )

write_rds(curve_vpd_sensitivity, "objects/curve-vpd-sensitivity.rds")

sensitivity_cor = cor(
  curve_vpd_sensitivity |> select(vpd_rate, h2os_rate, h2or_rate),
  use = "pairwise.complete.obs"
)
print(sensitivity_cor)

message(glue::glue(
  "Median H2O_s amplitude: {round(median(curve_vpd_sensitivity$h2os_amplitude, na.rm = TRUE), 2)} mmol/mol; ",
  "median H2O_r amplitude: {round(median(curve_vpd_sensitivity$h2or_amplitude, na.rm = TRUE), 2)} mmol/mol."
))

p_h2os_vs_vpd = ggplot(curve_vpd_sensitivity, aes(vpd_rate, h2os_rate)) +
  geom_point(alpha = 0.4) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
  labs(
    x = "VPD saturation rate (1/s)",
    y = expression(H[2] * O[s] ~ "saturation rate (1/s)")
  )

p_h2or_vs_vpd = ggplot(curve_vpd_sensitivity, aes(vpd_rate, h2or_rate)) +
  geom_point(alpha = 0.4) +
  labs(
    x = "VPD saturation rate (1/s)",
    y = expression(H[2] * O[r] ~ "saturation rate (1/s)")
  )

p_h2os_amplitude = ggplot(curve_vpd_sensitivity, aes(h2os_amplitude)) +
  geom_histogram(bins = 40) +
  scale_x_log10() +
  labs(x = expression(H[2] * O[s] ~ "amplitude |Asym - R0| (mmol/mol)"), y = "Number of curves")

p_h2or_amplitude = ggplot(curve_vpd_sensitivity, aes(h2or_amplitude)) +
  geom_histogram(bins = 40) +
  scale_x_log10() +
  labs(x = expression(H[2] * O[r] ~ "amplitude |Asym - R0| (mmol/mol)"), y = "Number of curves")

p_vpd_sensitivity = (p_h2os_vs_vpd + p_h2or_vs_vpd) /
  (p_h2os_amplitude + p_h2or_amplitude) +
  plot_annotation(tag_levels = "A")

print(p_vpd_sensitivity)

write_rds(p_vpd_sensitivity, "objects/plot-vpd-sensitivity.rds")

# Focused 2-panel version of just the rate-vs-rate comparison (no amplitude
# histograms), for use as a standalone figure in the response to reviewers
# (Comment R1.1): (A) VPD saturation rate vs. H2O_s saturation rate -- these
# track closely because both reflect gsw-driven change in the sample cell,
# and VPD is calculated directly from H2O_s; (B) VPD saturation rate vs.
# H2O_r saturation rate -- H2O_r is the incoming/reference airstream,
# upstream of the leaf, and is comparatively flat and uncorrelated with the
# VPD rate.
p_vpd_h2o_rate = (p_h2os_vs_vpd + p_h2or_vs_vpd) +
  plot_annotation(tag_levels = "A")

print(p_vpd_h2o_rate)

write_rds(p_vpd_h2o_rate, "objects/plot-vpd-h2o-rate.rds")

# --- Relationship between tau and VPD saturation rate ---------------------
#
# The VPD saturation rate k is strongly correlated with the fitted kinetics
# time-constant tau (see previous section: k tracks the same gsw-driven
# signal as H2O_s, not the largely independent H2O_r), which is why k is not
# used as a covariate for tau (that would be close to circular). This
# section quantifies and plots that relationship directly, joining the
# curve-level tau estimates (data/joined-summary.rds) to the VPD saturation
# rate (curve_vpd, by curve/id).

tau_vpd_rate = read_rds("data/joined-summary.rds") |>
  select(id, logtau_mean, leaf_type) |>
  left_join(
    curve_vpd |> select(curve, sat_rate, r2_sat, exclude_sat_rate),
    by = c("id" = "curve")
  )

tau_rate_cor = cor(
  tau_vpd_rate$logtau_mean,
  tau_vpd_rate$sat_rate,
  use = "pairwise.complete.obs"
)
message(glue::glue("Correlation between log(tau) and VPD saturation rate: {round(tau_rate_cor, 2)}"))

p_tau_vs_rate = ggplot(tau_vpd_rate, aes(sat_rate, logtau_mean, color = leaf_type)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "lm", se = FALSE) +
  scale_color_manual(values = c(amphi = col_amphi, pseudohypo = col_pseudohypo)) +
  labs(
    x = "VPD saturation rate, k (1/s)",
    y = expression(log ~ tau),
    color = "Leaf type"
  )

print(p_tau_vs_rate)

write_rds(p_tau_vs_rate, "objects/plot-tau-vpd-rate.rds")
