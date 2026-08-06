# Select "best_model" for parameter estimation and inference
source("r/header.R")

tmp = read_rds("objects/fit_01.rds") |>
  loo_moment_match()

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
  fit = glue("objects/fit_{n}.rds", n = str_extract(.m, "[0-9]+")) |>
    read_rds() |>
    as_draws_df()
  
  vars_in_model = colnames(fit)[colnames(fit) %in% focal_variables$var]
  
  fit |>
    select(vars_in_model) |>
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
  )

tmp = estimates_by_model2 |>
  full_join(focal_variables, by = join_by("variable" == "var")) |>
  mutate(prediction_correct = sign_observed == sign_expected) |>
  filter(!is.na(estimate)) 

tmp |>
  filter(model == "model_14") |> 
  select(response, explanatory, estimate, `q2.5`, `q97.5`, prediction_correct)

tmp |>
  summarize(
    n_focal_par = n(),
    n_signif = sum(!overlaps_zero),
    n_correct = sum(prediction_correct),
    .by = model
  )

selected_model = ?
  
tbl_comparison |> 
  mutate(selected = model == selected_model) |>
  write_rds("objects/tbl-comparison.rds")


write_rds(fits$fit[[as.numeric(str_remove(selected_model, "model"))]], "objects/best_model.rds")




# Write "best" model
# Note: I reran the models several times and the order of the top models
# changed, consistent with the difference in LOOIC being caused by sampling 
# variability. After reviewing model estimates, I determined that model 6 is the
# clearest to interpret.

# Previous version (lowest LOOIC)
# write_rds(fits$fit[[as.numeric(str_extract(rownames(looic_table)[1], "\\d+"))]], "objects/best_model.rds")

# Current version (model 6)
write_rds(fits$fit[[as.numeric(str_remove(selected_model, "model"))]], "objects/best_model.rds")


