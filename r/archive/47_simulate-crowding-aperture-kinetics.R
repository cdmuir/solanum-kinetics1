# --- Test a specific mechanistic hypothesis: could density-dependent
# stomatal "crowding" (interference among neighboring diffusion shells,
# stronger at high stomatal density -- e.g. Franks & Beerling 2009, de Boer
# et al. 2016) decouple fitted tau_hat from gmax while preserving its
# correlation with gi?
#
# Rationale: r/45-46 showed that a single GLOBAL nonlinearity parameter
# kappa (same for every curve) mainly amplifies the tau_hat/lambda_hat
# correlation with gi and f_gmax, without decoupling tau_hat from gmax.
# The user's hypothesis is different: if crowding gets systematically
# WORSE at higher stomatal density (not a fixed global departure from the
# geometric model), then high-gmax curves (which in this dataset are
# high-gmax largely BECAUSE they have high density; see below) would have
# a flatter gs(alpha) relationship right where they operate, damping the
# apparent kinetics there -- so tau_hat could stay correlated with gi
# (which sets the operating point on that curve's OWN gs(alpha) curve) but
# become decoupled from -- or even weakly negatively correlated with --
# gmax (which here is a stand-in for density itself).
#
# Model: same gleaf_op_nonlinear() as r/45, but kappa is now a per-curve
# function of that curve's stomatal density:
#   kappa_i = kappa0 * (d_i / d_ref)^p
# d_ref = median density in the dataset (so kappa0 is the crowding
# strength for a "typical" leaf); p = 1 (crowding severity scales linearly
# with density; not measured, illustrative). kappa0 = 0 recovers the
# pore-geometry-only model (r/41) for every curve.
#
# As in r/45, a single shared true kinetic rate constant k (biologically
# identical closure rate for every curve) is imposed; any resulting
# covariation of tau_hat/lambda_hat with gi/gmax is therefore attributable
# entirely to the density-dependent nonlinearity, not to a true difference
# in intrinsic closure rate.

source("r/header.R")

baseline_gec <- 0.02
baseline_xprime <- 6.8e-8

# Crowding-strength grid at a "typical" (median-density) leaf; kappa0 = 0
# reproduces r/41 (no crowding) for every curve regardless of density.
kappa0_grid <- c(0, 1, 3, 10)
p_crowding  <- 1

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
  select(id, d_m2, s_m2, gec, x_prime,
         gi = ginit_mean, gf = gfinal_mean, gmax, amax_um) |>
  left_join(
    pars_summary |> filter(variable == "sigma") |> select(id, sigma = mean),
    by = "id"
  ) |>
  filter(!is.na(sigma))

d_ref <- median(anat$d_m2)

sim_pars <- expand_grid(anat, kappa0 = kappa0_grid) |>
  mutate(
    kappa = kappa0 * (d_m2 / d_ref)^p_crowding,
    alpha_gi = pmap_dbl(list(gi, d_m2, s_m2, gec, x_prime, kappa), solve_alpha_nonlinear),
    alpha_gf = pmap_dbl(list(gf, d_m2, s_m2, gec, x_prime, kappa), solve_alpha_nonlinear)
  ) |>
  filter(!is.na(alpha_gi), !is.na(alpha_gf))

time_design <- joined_data |> select(curve, t_sec)

sim_one_crowding_curve <- function(id_val, kappa0_val, sim_pars, time_design, k) {
  p <- sim_pars[sim_pars$id == id_val & sim_pars$kappa0 == kappa0_val, ]
  ts <- time_design$t_sec[time_design$curve == id_val]
  alpha_t <- p$alpha_gf + (p$alpha_gi - p$alpha_gf) * exp(-k * ts)
  gsw_true <- gleaf_op_nonlinear(alpha_t, p$d_m2, p$s_m2, p$gec, p$x_prime, p$kappa)
  tibble(
    id = id_val,
    kappa0 = kappa0_val,
    t_sec = ts,
    alpha = alpha_t,
    gsw_true = gsw_true,
    gsw_sim = gsw_true + rnorm(length(ts), 0, p$sigma)
  )
}

set.seed(302947)
plan(multisession, workers = 9)

grid <- sim_pars |> distinct(id, kappa0)

sim_crowding_data <- future_map2_dfr(
  grid$id, grid$kappa0,
  sim_one_crowding_curve,
  sim_pars = sim_pars,
  time_design = time_design,
  k = k_shared,
  .progress = TRUE,
  .options = furrr_options(seed = TRUE)
)

write_rds(sim_crowding_data, "objects/sim-crowding-aperture-data.rds")
write_rds(
  list(sim_pars = sim_pars, k_shared = k_shared,
       gec = baseline_gec, x_prime = baseline_xprime,
       kappa0_grid = kappa0_grid, p_crowding = p_crowding, d_ref = d_ref),
  "objects/sim-crowding-aperture-pars.rds"
)
