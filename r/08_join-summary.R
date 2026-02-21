# Join stomatal kinetic data with stomatal anatomical data
source("r/header.R")

pars_summary = read_rds("objects/pars-summary.rds") |>
  filter(
    variable %in% c("ginit", "gfinal", "b_logtau_Intercept", "b_loglambda_Intercept"),
    rhat < convergence_criteria$rhat_max,
    ess_bulk > convergence_criteria$ess_min
  ) |>
  mutate(
    variable = case_when(
      variable == "gfinal" ~ "gfinal",
      variable == "ginit" ~ "ginit",
      variable == "b_logtau_Intercept" ~ "logtau",
      variable == "b_loglambda_Intercept" ~ "loglambda"
    )
  ) |>
  select(variable, mean, sd, id) |>
  pivot_wider(
    names_from = variable,
    values_from = c(mean, sd),
    names_glue = "{variable}_{.value}"
  ) |>
  mutate(
    accession = str_extract(id, "^(LA[0-9]{4}A{0,1}|nelsonii|sandwicense)"),
    replicate = str_extract(id, "-([A-Z][AB]{0,1})_", group = 1),
    curve_type = str_extract(id, "amphi|pseudohypo"),
    light_intensity = str_extract(id, "150|2000")
  )

joined_data = read_rds("data/joined-data.rds") |>
  summarize(
    guard_cell_length_um = first(guard_cell_length_um),
    gmax = first(gmax),
    .by = c(acc, id, curve_type)
  ) |>
  rename(accession = acc, replicate = id) |>
  mutate(
    curve_type = case_when(
      curve_type == "1-sided RH" ~ "amphi",
      curve_type == "2-sided RH" ~ "pseudohypo"
    )
  )

left_join(pars_summary,
          joined_data,
          by = join_by(accession, replicate, curve_type)) |>
  # missing kinetic data or any data contributed to gmax_ratio
  filter(!is.na(logtau_mean), !is.na(gmax)) |>
  left_join(
    read_rds("data/plant_info.rds") |>
      select(acc_id, accession, replicate, light_treatment),
    by = join_by(accession, replicate)
  ) |>
  mutate(
    light_treatment = light_treatment |>
      factor(levels = c("low", "high")) |>
      fct_recode(shade = "low", sun = "high"),
    light_intensity = light_intensity |>
      factor(levels = c("150", "2000")) |>
      fct_recode(low = "150", high = "2000"),
    f_gmax = ginit_mean / gmax
  ) |>
  write_rds("data/joined-summary.rds")
