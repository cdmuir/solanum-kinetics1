# Summarize RH and VPD treatments in a supporting table

source("r/header.R")

rh_curves = read_rds("data/rh_curves.rds")

curve_vpd = rh_curves |>
  mutate(SVPleaf = li6800_svp(Tleaf),
         VPDleaf = (SVPleaf - H2O_s * (Pa + ΔPcham) / 1000)) |>
  summarize(
    RH = median(RH),
    VPDleaf = median(VPDleaf),
    .by = c(acc_id, light_intensity, light_treatment, curve_type)
  ) |>
  mutate(
    accid = acc_id,
    lighttreatment = factor(case_when(
      light_treatment == "high" ~ "sun",
      light_treatment == "low" ~ "shade"
    ), levels = c("shade", "sun")),
    lightintensity = factor(case_when(
      light_intensity == "150" ~ "low",
      light_intensity == "2000" ~ "high"
    ), levels = c("low", "high")) |>
      factor(levels = c("low", "high")),
    curvetype = factor(case_when(
      curve_type == "1-sided RH" ~ "pseudohypo",
      curve_type == "2-sided RH" ~ "amphi"
    ), levels = c("amphi", "pseudohypo"))
  ) 
  
curve_vpd |>
  summarize(
    # min_RH = min(RH),
    # max_RH = max(RH),
    # min_VPDleaf = min(VPDleaf),
    # max_VPDleaf = max(VPDleaf),
    RH = mean(RH),
    VPDleaf = mean(VPDleaf),
    .by = c(lightintensity, lighttreatment, curvetype)
  ) |>
  mutate(
    `Growth\nlight intensity` = lighttreatment,
    `Measurement\nlight intensity` = lightintensity,
    `Curve type` = curvetype,
    `RH (\\%)` = prettyNum(100 * RH, digits = 3),
    `VPD (kPa)` = formatC(VPDleaf, format = "f", digits = 2)
  ) |>
  select(
    `Growth\nlight intensity`,
    `Measurement\nlight intensity`,
    `Curve type`,
    `RH (\\%)`,
    `VPD (kPa)`
  ) |>
  arrange(`Growth\nlight intensity`,
          `Measurement\nlight intensity`,
          desc(`Curve type`)) |>
  write_rds("objects/tbl-vpd.rds")

write_rds(curve_vpd, "objects/curve-vpd.rds")
