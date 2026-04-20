# Compare models using LOOIC
source("r/header.R")

#### ADD TO HEADER/FUNCTIONS IF RETAINING
library(reformulas)

get_terms_vec = function(fit, resp) {
  
  f = fit$formula$forms[[resp]]
  
  lhs = formula(f)[[2]]
  rhs = formula(f)[[3]]
  
  rhs_fixed = nobars(rhs)
  
  f_fixed = formula(paste(deparse(lhs), "~", paste(deparse(rhs_fixed), collapse = " ")))
  
  terms_vec = attr(terms(f_fixed), "term.labels")
  terms_vec
}
#####

plan(multisession, workers = 9)

file_names = list.files("objects/fits", full.names = TRUE)

loos = file_names |>
  future_map(\(.x) {
    fit = read_rds(.x) |>
      add_criterion(criterion = "loo")
    fit$criteria$loo
  }, .progress = TRUE)

looic_table = loos[-3744] |>
  set_names(file_names[-3744]) |>
  loo_compare() |>
  as.data.frame() |>
  rownames_to_column("file") |>
  filter(abs(elpd_diff) <= 2 * se_diff) 

write_rds(looic_table, "objects/looic_table-alt.rds")

df_terms_tau = file_names |>
  future_map_dfr(\(.x) {
    fit = read_rds(.x)
    terms_vec = get_terms_vec(fit, resp = "logtaumean")
    tibble(file = .x, terms = terms_vec)
  }, .progress = TRUE) |>
  mutate(value = TRUE) |>
  pivot_wider(names_from = terms, values_from = value, values_fill = FALSE)

df_terms_lambda = file_names |>
  future_map_dfr(\(.x) {
    fit = read_rds(.x)
    terms_vec = get_terms_vec(fit, resp = "loglambdamean")
    tibble(file = .x, terms = terms_vec)
  }, .progress = TRUE) |>
  mutate(value = TRUE) |>
  pivot_wider(names_from = terms, values_from = value, values_fill = FALSE)

write_rds(df_terms_tau, "objects/df_terms_tau.rds")
write_rds(df_terms_lambda, "objects/df_terms_lambda.rds")

# looic_table = read_rds("objects/looic_table-alt.rds")
# df_terms_tau = read_rds("objects/df_terms_tau.rds")
# df_terms_lambda = read_rds("objects/df_terms_lambda.rds")

looic_table |>
  left_join(df_terms_tau, by = "file") 

looic_table |>
  left_join(df_terms_lambda, by = "file")



  