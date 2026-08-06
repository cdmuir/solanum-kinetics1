# Compare estimates of focal parameters across plausible models
source("r/header.R")

# Plausible models
plausible_models = read_rds("objects/tbl-comparison.rds") |>
  filter(plausible) |>
  pull(model)

# Focal variables
# Phylo effect of gcl, gi, gmax on tau/lambda
# Fixed effect of gcl, gi, gmax on tau/lambda

focal_variables = crossing(
  type = c("b", "cor_phy"),
  response = c("logtaumean", "loglambdamean"),
  explanatory = c("loggcl", "loggi", "loggmax")
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
  mutate(model = factor(model, levels = model))
