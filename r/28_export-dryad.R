# Export text-based data sets to deposit on Dryad
source("r/header.R")

fit = read_rds("objects/best_model.rds")

fit$data |>
  select(
    accid,
    lighttreatment,
    lightintensity,
    curve_type,
    loggcl,
    logitfgmax,
    loglambdamean,
    loglambdasd,
    logtaumean,
    logtausd
  ) |>
  write_csv("dryad/main.csv")

# Export README once data sets are ready