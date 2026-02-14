# Placeholder for now with some prelim code
source("r/header.R")

fit_amphi_high = read_rds("objects/fit_amphi_high.rds")
fit_amphi_low = read_rds("objects/fit_amphi_low.rds")

# Variance decomposition
fit_amphi_low |>
  as_draws_df() |>
  select(starts_with("."), starts_with("sd_"), starts_with("sigma")) |>
  rename_with(.fn = \(.x) {str_replace(.x, "sigma_", "sd_resid__")}, .cols = starts_with("sigma_")) |>
  rename_with(.fn = \(.x) {str_remove(.x, "_Intercept")}, .cols = ends_with("_Intercept")) |>
  pivot_longer(
    cols = -starts_with("."),
    names_sep = "__",
    names_to = c("component", "trait"),
    values_to = "sd"
  ) |>
  mutate(component = str_remove(component, "sd_")) |>
  mutate(var = sd ^ 2, .keep = "unused") |>
  pivot_wider(values_from = var, names_from = component) |>
  mutate(total_var = phy + accession + resid) |>
  mutate(across(accession:resid, ~ .x / total_var)) |>
  split(~ trait) |>
  map(summarize_draws) |>
  map(filter, variable != "trait") |>
  imap_dfr(\(.x, .y) {
    .x |> mutate(trait = .y)
  }) |>
  select(variable, median, q5, q95, trait)
