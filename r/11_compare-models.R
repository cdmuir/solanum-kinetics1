# --- Compare candidate models using LOOIC ------------------------------------

source("r/header.R")

plan(multisession, workers = 19)

fits = read_rds("objects/df_forms.rds") |>
  mutate(
    fit = future_map(file, read_rds),
    model = file |>
      str_replace("objects/fit", "model") |>
      str_remove(".rds")
  ) |>
  filter(!is.na(fit)) |>
  mutate(loo = map(fit, \(.x) .x$criteria$loo))

converged = fits$fit |>
  future_map_lgl(check_convergence, convergence_criteria, .progress = TRUE)

looic_table = fits$loo |>
  set_names(fits$model) |>
  loo_compare() 

map2_dfr(fits$fit, fits$model, \(.fit, .name) {
  tibble(par = .fit |>
           as_draws_df() |>
           dplyr::select(contains("b_")) |>
           colnames()) |>
    mutate(par = str_remove(par, "b_") |>
             str_replace("curve_type", "curvetype")) |>
    separate_wider_delim(par,
                         delim = "_",
                         names = c("response", "explanatory")) |>
    filter(
      str_detect(response, "^log(lambda|tau)mean$"),
      explanatory %in% c("loggcl", "loggi", "loggmax")
    ) |>
    mutate(
      response = fct_recode(response, `$\\lambda$` = "loglambdamean", `$\\tau$` = "logtaumean"),
      explanatory = fct_recode(explanatory, `\\gcl` = "loggcl", `\\gi` = "loggi", `\\gmax` = "loggmax"),
      model = .name
    )
  
}) |>
  mutate(value = TRUE) |>
  pivot_wider(id_cols = c(response, model), names_from = explanatory) |>
  pivot_wider(
    id_cols = model,
    names_from = response,
    values_from = c(`\\gcl`, `\\gi`, `\\gmax`),
    names_glue = "{response}__{.value}"
  ) |>
  full_join(tibble(
    model = looic_table$model,
    `$\\Delta \\mathrm{LOOIC}$` = -2 * looic_table[, "elpd_diff"],
    SE = 2 * looic_table[, "se_diff"]
  ),
  by = join_by(model)) |>
  mutate(across(where(is_logical), \(.x) ifelse(.x, "\\cmark", "")),
         plausible = abs(`$\\Delta \\mathrm{LOOIC}$`) <= 2 * SE
         ) |>
  arrange(`$\\Delta \\mathrm{LOOIC}$`) |>
  mutate(
    `$\\Delta \\mathrm{LOOIC}$` = formatC(
      abs(`$\\Delta \\mathrm{LOOIC}$`),
      format = "f",
      digits = 2
    ),
    SE = formatC(SE, format = "f", digits = 2),
    mutate(across(everything(), \(.x) replace_na(.x, "")))
  ) |>
  write_rds("objects/tbl-comparison.rds")
