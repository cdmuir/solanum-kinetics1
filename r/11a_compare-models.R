# Compare models using LOOIC
source("r/header.R")

# selected_model = "model6" # see note at bottom

plan(multisession, workers = 9)

# TODO: REDO AND MAKE COMPACT ONCE ALL FITS ARE DONE
fits = read_rds("objects/df_forms.rds") |>
  mutate(
    fit = map(file, \(.x) {
      if (file.exists(.x)) {
        read_rds(.x)
      } else {
        NA
      }
    }),
    model = file |>
      str_replace("objects/fit", "model") |>
      str_remove(".rds")
  ) |>
  filter(!is.na(fit)) |>
  mutate(loo = map(fit, \(.x) .x$criteria$loo))

converged = df_forms$fit |>
  future_map_lgl(check_convergence, convergence_criteria)

assert_true(all(converged))

looic_table = fits$loo |>
  set_names(fits$model) |>
  loo_compare() 

tmp = map2_dfr(fits$fit, fits$model, \(.fit, .name) {
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
         plausible = abs(`$\\Delta \\mathrm{LOOIC}$`) <= 2 * SE,
         # selected = model == selected_model
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


# Write "best" model
# Note: I reran the models several times and the order of the top four models
# changed, consistent with the difference in LOOIC being caused by sampling 
# variability. After reviewing model estimates, I determined that model 6 is the
# clearest to interpret.

# Previous version (lowest LOOIC)
# write_rds(fits$fit[[as.numeric(str_extract(rownames(looic_table)[1], "\\d+"))]], "objects/best_model.rds")

# Current version (model 6)
write_rds(fits$fit[[as.numeric(str_remove(selected_model, "model"))]], "objects/best_model.rds")
