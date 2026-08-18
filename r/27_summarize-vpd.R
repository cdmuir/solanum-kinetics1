# --- Summarize the VPD trajectory of each curve with a saturating exponential ---
#
# VPD(t) = Asym + (R0 - Asym) * exp(-rate * t_sec) is used instead of a linear
# slope: raw trajectories rise then plateau, and the saturating model fits far
# better (median R2 = 0.99 vs. 0.89 for a straight line, across all curves).
# nlsLM() (Levenberg-Marquardt) with data-driven starting values is used
# because SSasymp()'s default self-start algorithm failed to converge for
# 10/2156 curves; the starting values below rescue all of them.

source("r/header.R")

rh_curves = read_rds("data/rh_curves.rds")

vpd_traj = rh_curves |>
  mutate(
    SVPleaf = li6800_svp(Tleaf),
    VPDleaf = SVPleaf - H2O_s * (Pa + ΔPcham) / 1000
  )

vpd_saturating = vpd_traj |>
  group_by(curve) |>
  group_modify(~ fit_saturating(.x, "VPDleaf")) |>
  ungroup()

# All curves should converge; flag any that don't or fit poorly so they can
# be inspected before use in the model.
n_failed = sum(!vpd_saturating$converged)
n_poor_fit = sum(vpd_saturating$r2_sat < 0.8, na.rm = TRUE)
if (n_failed > 0 || n_poor_fit > 0) {
  warning(glue::glue(
    "{n_failed} curve(s) failed to converge and {n_poor_fit} curve(s) have ",
    "R2 < 0.8 for the saturating VPD fit. Inspect objects/curve-vpd-summary.rds."
  ))
}

# Two curves have implausibly large sat_rate estimates (~0.04-0.05 /s, a
# ~2x jump from the next-highest values, ~0.02 /s) despite acceptable R2
# (0.84-0.86). Both are unusually short curves (12-13 observations), where
# the rate constant is poorly identified. 

# Eight curves have extremely small sat_rate estimates (<0.0001 /s) despite 
# acceptable R2 (0.81-0.997) in 7/8 cases. Excluded here, alongside curves
# with r2_sat < 0.8.

sat_rate_outliers = c("LA0407-J_pseudohypo_150", "LA0436-A_pseudohypo_150",
                      "LA1364-K_amphi_150", "LA1782-T_pseudohypo_150",
                      "LA2133-J_amphi_2000", "LA2951-R_amphi_2000",
                      "LA2951-T_amphi_150", "LA2951-W_pseudohypo_150",
                      "LA2951-Z_amphi_150", "LA2964-U_pseudohypo_2000")

curve_vpd = vpd_traj |>
  summarize(
    median_vpd = median(VPDleaf),
    # last observation in the trimmed curve (largest t_sec)
    final_vpd = VPDleaf[which.max(t_sec)],
    n_obs = n(),
    duration_sec = max(t_sec),
    .by = c(curve, acc_id, light_intensity, light_treatment, curve_type)
  ) |>
  left_join(vpd_saturating, by = "curve") |>
  mutate(
    exclude_sat_rate = r2_sat < 0.8 | curve %in% sat_rate_outliers,
    sat_rate = if_else(exclude_sat_rate, NA_real_, sat_rate),
    sat_rate_se = if_else(exclude_sat_rate, NA_real_, sat_rate_se),
    lighttreatment = recode_lighttreatment(light_treatment),
    lightintensity = recode_lightintensity(light_intensity),
    leaftype = recode_leaftype(curve_type)
  ) |>
  rename(accid = acc_id)

write_rds(curve_vpd, "objects/curve-vpd-summary.rds")
