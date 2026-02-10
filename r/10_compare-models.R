# Compare models using LOOIC
source("r/header.R")

plan(multisession, workers = 16)

fits_amphi = read_rds("objects/fits_amphi.rds") |>
  mutate(loo = map(fit, \(.x) .x$criteria$loo))

converged = fits_amphi$fit |>
  future_map_lgl(check_convergence, convergence_criteria)

assert_true(all(converged))

looic_table = fits_amphi$loo |>
  set_names(paste0("model", seq_along(fits_amphi$loo))) |>
  loo_compare()

write_rds(
  fits_amphi$fit[[as.numeric(str_extract(rownames(looic_table)[1], "\\d+"))]],
  "objects/best_amphi_model.rds"
)
