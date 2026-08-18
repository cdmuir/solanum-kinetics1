# --- Summarize realized RH and VPD by treatment combination for a supporting table ---

source("r/header.R")

rh_curves = read_rds("data/rh_curves.rds")

curve_vpd = rh_curves |>
  mutate(SVPleaf = li6800_svp(Tleaf),
         VPDleaf = (SVPleaf - H2O_s * (Pa + ΔPcham) / 1000)) |>
  summarize(
    across(c(RH, VPDleaf), list(initial = first, median = median, final = last),
           .names = "{.fn}_{.col}"),
    .by = c(acc_id, light_intensity, light_treatment, curve_type)
  ) |>
  mutate(
    accid = acc_id,
    lighttreatment = recode_lighttreatment(light_treatment),
    lightintensity = recode_lightintensity(light_intensity),
    leaftype = recode_leaftype(curve_type)
  ) 
  
curve_vpd |>
  summarize(
    across(contains("RH"), mean),
    across(contains("VPD"), mean),
    .by = c(lightintensity, lighttreatment, leaftype)
  ) |>
  mutate(
    `Growth\nlight intensity` = lighttreatment,
    `Measurement\nlight intensity` = lightintensity,
    `Leaf type` = leaftype,
    RH_initial = prettyNum(100 * initial_RH, digits = 3),
    RH_median = prettyNum(100 * median_RH, digits = 3),
    RH_final = prettyNum(100 * final_RH, digits = 3),
    VPD_initial = formatC(initial_VPDleaf, format = "f", digits = 2),
    VPD_median = formatC(median_VPDleaf, format = "f", digits = 2),
    VPD_final = formatC(final_VPDleaf, format = "f", digits = 2)
  ) |>
  select(
    `Growth\nlight intensity`,
    `Measurement\nlight intensity`,
    `Leaf type`,
    RH_initial, RH_median, RH_final,
    VPD_initial, VPD_median, VPD_final
  ) |>
  arrange(`Growth\nlight intensity`,
          `Measurement\nlight intensity`,
          desc(`Leaf type`)) |>
  write_rds("objects/tbl-vpd.rds")

write_rds(curve_vpd, "objects/curve-vpd.rds")
