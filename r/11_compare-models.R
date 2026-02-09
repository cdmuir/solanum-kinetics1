# i did this at work. might already have version going at home? merge/delete
source("r/header.R")
fits_amphi = read_rds("objects/fits_amphi.rds") |>
  mutate(loo = map(fit, \(.x) .x$criteria$loo))

fits_amphi$loo |>
  set_names(paste0("model", seq_along(fits_amphi$loo))) |>
loo_compare()
fits_amphi$fit[[2]]
