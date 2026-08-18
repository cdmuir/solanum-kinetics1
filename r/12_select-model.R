# --- Select the model used for parameter estimation and inference ------------
#
# Criteria: (1) the model is among the plausible set based on LOOIC, and (2) it
# has the highest proportion of significant focal parameters based on 95% CIs.

source("r/header.R")

# Plausible models
plausible_models = read_rds("objects/tbl-comparison.rds") |>
  filter(plausible) |>
  pull(model)

tbl_comparison = read_rds("objects/tbl-comparison.rds")

# Focal variables
# Phylo effect of gcl, gi, gmax on tau/lambda
# Fixed effect of gcl, gi, gmax on tau/lambda

focal_variables = crossing(
  type = c("b", "cor_phy"),
  response = c("logtaumean", "loglambdamean"),
  nesting(
    explanatory = c("loggcl", "loggi", "loggmax"),
    sign_expected = c(1, 1, -1)
  )
) |>
  mutate(
    response = case_when(
      type == "b" ~ response,
      type == "cor_phy" ~ str_c(response, "Intercept", sep = "_")
    ),
    explanatory = case_when(
      type == "b" ~ explanatory,
      type == "cor_phy" ~ str_c(explanatory, "Intercept", sep = "_")
    ),
    sep = case_when(type == "b" ~ "_", type == "cor_phy" ~ "__"),
    var = glue("{type}{sep}{response}{sep}{explanatory}")
  )

# Get focal parameter estimates and 95% CIs from each plausible model
estimates_by_model1 = map_dfr(plausible_models, \(.m) {
  fit = glue("objects/fits/fit_{n}.rds", n = str_extract(.m, "[0-9]+")) |>
    read_rds() |>
    as_draws_df()
  
  vars_in_model = colnames(fit)[colnames(fit) %in% focal_variables$var]
  
  fit |>
    dplyr::select(vars_in_model) |>
    summarize_draws(estimate = median, ~ quantile2(.x, probs = c(0.025, 0.975))) |>
    mutate(model = .m)
})

estimates_by_model2 = full_join(
  estimates_by_model1,
  expand(estimates_by_model1, variable, model),
  by = join_by(variable, model)
) |>
  mutate(
    model = factor(model, levels = plausible_models),
    overlaps_zero = sign(`q2.5`) != sign(`q97.5`),
    sign_observed = ifelse(overlaps_zero, 0, sign(estimate))
  ) |>
  full_join(focal_variables, by = join_by("variable" == "var")) |>
  mutate(prediction_correct = sign_observed == sign_expected)

# Select model
selected_model = estimates_by_model2 |>
  filter(!is.na(estimate)) |>
  summarize(
    n_focal_par = n(),
    n_signif = sum(!overlaps_zero),
    n_correct = sum(prediction_correct),
    .by = model
  ) |>
  # Sort by proportion significant, tie break by n_focal_par
  arrange(desc(n_signif / n_focal_par), n_focal_par) |>
  pull(model) |>
  first() |>
  as.character()

# Write estimates_by_model2 for plotting
write_rds(estimates_by_model2, "objects/estimates_by_model.rds")

# Update tbl_comparison
tbl_comparison |>
  mutate(selected = model == selected_model) |>
  write_rds("objects/tbl-comparison.rds")

# Write selected model
selected_model |>
  str_replace("^models/fit_", "objects/fits/fit_") |>
  str_c(".rds") |>
  file.copy("objects/selected_model.rds", overwrite = TRUE)
