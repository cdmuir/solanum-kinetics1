# --- Make the per-curve table of kinetic and anatomical parameter estimates ---

source("r/header.R")

selected_model = read_rds("objects/selected_model.rds")

curve_vpd = read_rds("objects/curve-vpd.rds") |>
  select(
    accid,
    lighttreatment,
    lightintensity,
    leaftype,
    initial_RH,
    median_RH,
    final_RH,
    initial_VPD = initial_VPDleaf,
    median_VPD = median_VPDleaf,
    final_VPD = final_VPDleaf
  )

tbl_estimates = selected_model$data |>
  select(
    accid,
    lighttreatment,
    lightintensity,
    leaftype,
    loggcl,
    loggi,
    loggmax,
    loglambdamean,
    loglambdasd,
    logtaumean,
    logtausd
  ) |>
  left_join(curve_vpd, by = join_by(accid, lighttreatment, lightintensity, leaftype)) |>
  separate_wider_delim(
    cols = accid,
    delim = "-",
    names = c("accession", "id")
  )

write_csv(tbl_estimates, "tables/tbl-estimates-curve.csv")

# Write tbl-estimates-curve-dictionary.csv

# ---- build dictionary ----
dict = tibble(
  variable = names(tbl_estimates),
  data_type = map_chr(tbl_estimates, type_label),
  acceptable_values = map_chr(tbl_estimates, acceptable_values),
) |>
  mutate(acceptable_values = case_when(
    variable == "id" ~ "A-Z",
    TRUE ~ acceptable_values
  )) |>
  full_join(tibble(
    variable = c(
      "accession",
      "id",
      "lighttreatment",
      "lightintensity",
      "leaftype",
      "loggcl",
      "loggi",
      "loggmax",
      "loglambdamean",
      "loglambdasd",
      "logtaumean",
      "logtausd",
      "initial_RH",
      "median_RH",
      "final_RH",
      "initial_VPD",
      "median_VPD",
      "final_VPD"
    ),
    description = c(
      "TGRC accession",
      "Unique identifier for each accession",
      "Growth light intensity treatment",
      "Measurement light intensity treatment",
      "Leaf type treatment",
      "Estimate of log-transformed guard cell length (um)",
      "Estimate of log-transformed initial stomatal conductance (mol m^-2 s^-1)",
      "Estimate of log-transformed anatomical maximum conductance (mol m^-2 s^-1)",
      "Point estimate of logit-transformed lag-time parameter",
      "Standard error of logit-transformed lag-time parameter",
      "Point estimate of log-transformed time-constant parameter",
      "Standard error of log-transformed time-constant parameter",
      "Relative humidity (proportion) at the start of the fitted interval",
      "Median relative humidity (proportion) during the fitted interval",
      "Relative humidity (proportion) at the end of the fitted interval",
      "Leaf-to-air vapor pressure deficit (kPa) at the start of the fitted interval",
      "Median leaf-to-air vapor pressure deficit (kPa) during the fitted interval",
      "Leaf-to-air vapor pressure deficit (kPa) at the end of the fitted interval"
    ),
  ), by = join_by(variable)) |>
  mutate(
    acceptable_values = na_if(acceptable_values, ""),
    acceptable_values = replace_na(acceptable_values, "")
  )

# ---- write to CSV ----
write_csv(dict, "tables/tbl-estimates-curve-dictionary.csv")
