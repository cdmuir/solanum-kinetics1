# Compare models using LOOIC
source("r/header.R")

plan(multisession, workers = 19)

fits = read_rds("objects/fits.rds") |>
  mutate(loo = map(fit, \(.x) .x$criteria$loo))

converged = fits$fit |>
  future_map_lgl(check_convergence, convergence_criteria)

assert_true(all(converged))

looic_table = fits$loo |>
  set_names(paste0("model", seq_along(fits$loo))) |>
  loo_compare()

fits = fits |>
  mutate(model = paste0("model", row_number()))

map2_dfr(fits$fit, fits$model, \(.fit, .name) {
  tibble(par = .fit |>
           as_draws_df() |>
           select(contains("b_")) |>
           colnames()) |>
    mutate(par = str_remove(par, "b_") |>
             str_replace("curve_type", "curvetype")) |>
    separate_wider_delim(par,
                         delim = "_",
                         names = c("response", "explanatory")) |>
    filter(
      str_detect(response, "^log(lambda|tau)mean$"),
      explanatory %in% c("logitfgmax", "loggcl")
    ) |>
    mutate(
      response = fct_recode(response, `$\\lambda$` = "loglambdamean", `$\\tau$` = "logtaumean"),
      explanatory = fct_recode(explanatory, `\\gcl` = "loggcl", `\\fgmax` = "logitfgmax"),
      model = .name
    )
  
}) |>
  mutate(value = TRUE) |>
  pivot_wider(id_cols = c(response, model), names_from = explanatory) |>
  pivot_wider(
    id_cols = model,
    names_from = response,
    values_from = c(`\\gcl`, `\\fgmax`),
    names_glue = "{response}__{.value}"
  ) |>
  full_join(tibble(
    model = rownames(looic_table),
    `$\\Delta \\text{LOOIC}$` = -2 * looic_table[, "elpd_diff"],
    SE = 2 * looic_table[, "se_diff"]
  ),
  by = "model") |>
  mutate(across(where(is_logical), \(.x) ifelse(.x, "\\cmark", ""))) |>
  arrange(`$\\Delta \\text{LOOIC}$`) |>
  mutate(
    `$\\Delta \\text{LOOIC}$` = formatC(
      abs(`$\\Delta \\text{LOOIC}$`),
      format = "f",
      digits = 2
    ),
    SE = formatC(SE, format = "f", digits = 2),
    mutate(across(everything(), \(.x) replace_na(.x, "")))
  ) |>
  write_rds("objects/tbl-comparison.rds")


# Write best model
write_rds(fits$fit[[as.numeric(str_extract(rownames(looic_table)[1], "\\d+"))]], "objects/best_model.rds")
