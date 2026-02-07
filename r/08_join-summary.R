# NEED TO FIX CALCULATION OF gsmax (divide by 1000 and make different for lower versus total)
# NEED TO DECIDE WHETHER THIS HAPPENS HERE OR IN 01_join-data.R
# Join stomatal kinetic data with stomatal anatomical data
source("r/header.R")

full_join(
  read_rds("objects/pars-summary.rds") |>
    filter(
      variable %in% c("ginit", "b_logtau_Intercept", "b_loglambda_Intercept"),
      rhat < convergence_criteria$rhat_max,
      ess_bulk > convergence_criteria$ess_min
    ) |>
    mutate(
      variable = case_when(
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
    ),
  read_rds("data/stomata.rds"),
  by = join_by(accession, replicate)
) |>
  # missing kinetic data or any data contributed to gmax_ratio
  filter(!is.na(logtau_mean), !is.na(gmax_ratio)) |>
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
    guard_cell_length_um = lower_guard_cell_length_um * (1 - stomatal_ratio) + upper_guard_cell_length_um * stomatal_ratio,
    total_gmax = lower_gmax + upper_gmax,
    f_gmax = ginit_mean / total_gmax
  ) |>
  write_rds("data/joined-summary.rds")

foo |>
  filter(acc_id == "LA1777-B") |>
  select(id, lower_gmax, upper_gmax)
