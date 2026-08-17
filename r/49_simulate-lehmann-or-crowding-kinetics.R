# --- Repeat the density-dependent crowding test from r/47-48, but with a
# GEOMETRY-DERIVED crowding strength (informed by Lehmann & Or 2015,
# gleaf_op_lehmann_or() in r/functions.R) rather than the ad hoc
# kappa0 * (d/d_ref)^1 form used in r/47.
#
# Lehmann & Or's mechanism ties crowding-induced diffusive resistance to
# the ratio of pore size to inter-stomatal spacing. Here that ratio,
# chi(alpha) = 2 * alpha * (amax/2) / (1/sqrt(d)), is computed directly
# from each curve's own anatomy (density AND pore size, not density
# alone), and multiplied by a single free calibration constant kappa0
# (kappa0 = 0 recovers gleaf_op() exactly for every curve). See r/48's
# side-check output for the typical range of amax/spacing ratios in this
# dataset (median ~0.10, up to ~0.20 at full aperture) which is used here
# to choose a plausible kappa0 grid.
#
# As in r/45-48, a single shared true kinetic rate constant k (identical
# closure rate for every curve) is imposed; any resulting covariation of
# tau_hat/lambda_hat with gi/gmax is therefore attributable entirely to
# the geometry-derived crowding nonlinearity.

source("r/header.R")

baseline_gec <- 0.02
baseline_xprime <- 6.8e-8

# kappa0 grid: since the geometric ratio r_max/lambda_s is already ~0.04-
# 0.2 across curves (see r/47/48 companion check below), kappa0 on the
# order of 1-50 is needed to span from negligible to strong crowding
# (kappa0 * chi(alpha = 1) ranging from ~0.1 to ~4 across this dataset).
kappa0_grid <- c(0, 5, 15, 50)

tbl_implied_aperture <- read_rds("objects/tbl-implied-aperture.rds")
pars_summary <- read_rds("objects/pars-summary.rds")
joined_data <- read_rds("data/joined-data.rds")

median_tau <- pars_summary |>
  filter(variable == "b_logtau_Intercept") |>
  summarise(median_tau = median(exp(mean))) |>
  pull(median_tau)

k_shared <- 1 / median_tau

anat <- tbl_implied_aperture |>
  filter(gec == baseline_gec, x_prime == baseline_xprime) |>
  select(id, d_m2, s_m2, gec, x_prime, amax_um,
         gi = ginit_mean, gf = gfinal_mean, gmax) |>
  mutate(amax_m = amax_um * 1e-6) |>
  left_join(
    pars_summary |> filter(variable == "sigma") |> select(id, sigma = mean),
    by = "id"
  ) |>
  filter(!is.na(sigma))

# Side check: typical pore-radius-to-spacing ratio at full aperture
# (chi(alpha = 1) / kappa0), for reference when interpreting kappa0.
r_over_spacing <- anat |>
  mutate(ratio = amax_m / 2 * sqrt(d_m2)) |>
  pull(ratio)
cli::cli_alert_info(
  "amax/2 * sqrt(d) (pore radius / spacing) ranges {round(min(r_over_spacing), 3)}-{round(max(r_over_spacing), 3)}, median {round(median(r_over_spacing), 3)}"
)

sim_pars <- expand_grid(anat, kappa0 = kappa0_grid) |>
  mutate(
    alpha_gi = pmap_dbl(list(gi, d_m2, s_m2, gec, x_prime, amax_m, kappa0), solve_alpha_lehmann_or),
    alpha_gf = pmap_dbl(list(gf, d_m2, s_m2, gec, x_prime, amax_m, kappa0), solve_alpha_lehmann_or)
  ) |>
  filter(!is.na(alpha_gi), !is.na(alpha_gf))

time_design <- joined_data |> select(curve, t_sec)

sim_one_lehmannor_curve <- function(id_val, kappa0_val, sim_pars, time_design, k) {
  p <- sim_pars[sim_pars$id == id_val & sim_pars$kappa0 == kappa0_val, ]
  ts <- time_design$t_sec[time_design$curve == id_val]
  alpha_t <- p$alpha_gf + (p$alpha_gi - p$alpha_gf) * exp(-k * ts)
  gsw_true <- gleaf_op_lehmann_or(alpha_t, p$d_m2, p$s_m2, p$gec, p$x_prime, p$amax_m, p$kappa0)
  tibble(
    id = id_val,
    kappa0 = kappa0_val,
    t_sec = ts,
    alpha = alpha_t,
    gsw_true = gsw_true,
    gsw_sim = gsw_true + rnorm(length(ts), 0, p$sigma)
  )
}

set.seed(148203)
plan(multisession, workers = 9)

grid <- sim_pars |> distinct(id, kappa0)

sim_lehmannor_data <- future_map2_dfr(
  grid$id, grid$kappa0,
  sim_one_lehmannor_curve,
  sim_pars = sim_pars,
  time_design = time_design,
  k = k_shared,
  .progress = TRUE,
  .options = furrr_options(seed = TRUE)
)

write_rds(sim_lehmannor_data, "objects/sim-lehmannor-crowding-data.rds")
write_rds(
  list(sim_pars = sim_pars, k_shared = k_shared,
       gec = baseline_gec, x_prime = baseline_xprime, kappa0_grid = kappa0_grid),
  "objects/sim-lehmannor-crowding-pars.rds"
)
