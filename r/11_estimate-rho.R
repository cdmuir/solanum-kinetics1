# Estimate evolutionary correlation between stomatal kinetics, guard cell size, and fraction gmax
source("r/header.R")

fit_amphi_high = read_rds("objects/fit_amphi_high.rds")
fit_amphi_low = read_rds("objects/fit_amphi_low.rds")

# Partial correlation
rho_high = summarize_parcor(fit_amphi_high)
rho_low = summarize_parcor(fit_amphi_low)

rho_high |>
  select(component, trait1, trait2, median, q5, q95)

rho_low |>
  select(component, trait1, trait2, median, q5, q95)
