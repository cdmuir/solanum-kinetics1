# Join stomatal kinetic data with stomatal anatomical data
source("r/header.R")

# NEED TO ADD STEP FILTERING OUT NON-CONVERGED MODELS
tmp = full_join(
  read_rds("objects/brms-summary.rds") |>
    filter(variable == "b_logtau_Intercept") |>
    select(mean, sd, id) |>
    mutate(
      accession = str_extract(id, "^(LA[0-9]{4}A{0,1}|nelsonii|sandwicense)"),
      replicate = str_extract(id, "-([A-Z][AB]{0,1})_", group = 1),
      curve_type = str_extract(id, "amphi|pseudohypo"),
      light_intensity = str_extract(id, "150|2000")
    ),
  read_rds("data/stomata.rds"),
  by = join_by(accession, replicate)
)

tmp |>
  filter(is.na(mean))

tmp |>
  filter(is.na(upper_stomatal_density_mm2))
