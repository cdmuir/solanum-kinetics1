# Test whether VPD is associated with tau
source("r/header.R")

joined_summary = read_rds("data/joined-summary.rds") |>
  prepare_tau_anatomy_data(logtau_threshold)

curve_vpd = read_rds("objects/curve-vpd.rds")

joined_summary1 = joined_summary |>
  left_join(curve_vpd,
            by = join_by(curvetype, lightintensity, accid, lighttreatment)) |>
  filter(!is.na(VPDleaf))

phy = read_rds("data/phylogeny.rds")
A = vcv(phy, corr = TRUE)
thin = 1

# Define formula
form_vpd = bf(
  logtaumean | se(logtausd, sigma = TRUE) ~
    lighttreatment +
    lightintensity +
    curvetype +
    logitfgmax + 
    VPDleaf +
    (1 | accid) +
    (1 | a | accession) +
    (1 | b | gr(phy, cov = A))
)

fit_vpd = brm(
  formula = form_vpd,
  data = joined_summary1 |>
    mutate(phy = accession),
  data2 = list(A = A),
  cores = 4,
  chains = 4,
  iter = thin * 2e3,
  thin = thin,
  refresh = thin * 1e2,
  control = list(adapt_delta = 0.9),
  backend = "cmdstanr",
  family = student(),
  seed = 30420
)

write_rds(fit_vpd, "objects/fit-vpd.rds")

fit_vpd = read_rds("objects/fit-vpd.rds")
check_convergence(fit_vpd, convergence_criteria)