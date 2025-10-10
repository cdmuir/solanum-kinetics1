# Join stomatal kinetic data with stomatal anatomical data
source("r/header.R")

full_join(
  read_rds("objects/brms-summary.rds") |>
    filter(
      variable %in% c("b_logtau_Intercept", "log_gi"),
      rhat < convergence_criteria$rhat_max,
      ess_bulk > convergence_criteria$ess_min
    ) |>
    mutate(
      variable = case_when(
        variable == "b_logtau_Intercept" ~ "log_tau",
        variable == "log_gi" ~ "log_gi"
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
  filter(!is.na(log_tau_mean), !is.na(gmax_ratio)) |>
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
    f_gmax = exp(log_gi_mean) / total_gmax
  ) |>
  write_rds("data/joined-data.rds")
