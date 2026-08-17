# --- Simulate gsw(t) from the same anatomy-independent aperture kinetics
# model as r/41_simulate-aperture-kinetics.R, but mapping aperture to
# conductance through gleaf_op_nonlinear() (r/functions.R) instead of the
# pore-geometry-only gleaf_op(). gleaf_op_nonlinear() matches the Ochoa et
# al. (2024) pore-diffusion model exactly at small fractional aperture, but
# increasingly saturates below it as aperture grows, with the degree of
# saturation controlled by a free parameter kappa (kappa = 0 recovers
# gleaf_op() exactly).
#
# Question: if every curve shares the SAME underlying pore-closure rate
# constant k (no true dependence of aperture kinetics on gi, gmax, or
# f_gmax), does an increasingly nonlinear (saturating) aperture-conductance
# relationship change whether/how strongly the fitted Weibull tau_hat and
# lambda_hat covary with gi and gmax? This is the same question as r/41-42,
# but for a phenomenological nonlinearity in the aperture-conductance
# mapping rather than a physically motivated series resistance (contrast
# with r/44_simulate-gias-series.R).
#
# Design mirrors r/41 exactly except:
#   - alpha_gi, alpha_gf are recalculated per kappa via
#     solve_alpha_nonlinear() (not reused from r/40's tbl_implied_aperture,
#     which was calibrated under the pore-geometry-only model, kappa = 0).
#   - Simulated gsw(t) uses gleaf_op_nonlinear() with that curve's kappa.
#   - gec and x_prime are fixed at the same baseline combination as r/41
#     (gec = 0.02 mol m^-2 s^-1, x_prime = 6.8e-8 m).
#   - k_shared is identical to r/41 (1 / median real tau); it does not
#     depend on kappa, so any change in tau_hat/lambda_hat covariation
#     across the kappa grid is attributable to the nonlinearity alone.

source("r/header.R")

baseline_gec <- 0.02
baseline_xprime <- 6.8e-8

# kappa grid: 0 = pore-geometry-only baseline (identical to r/41); larger
# kappa = increasingly strong saturation of gs relative to alpha at large
# aperture (alpha_eff = 1 / (1 + kappa) at full aperture, e.g. kappa = 10
# caps the effective pore-diffusion aperture at ~9% of true aperture even
# when the pore is fully open). Not measured; illustrative only.
kappa_grid <- c(0, 1, 3, 10)

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

sim_pars <- expand_grid(anat, kappa = kappa_grid) |>
  mutate(
    alpha_gi = pmap_dbl(list(gi, d_m2, s_m2, gec, x_prime, kappa), solve_alpha_nonlinear),
    alpha_gf = pmap_dbl(list(gf, d_m2, s_m2, gec, x_prime, kappa), solve_alpha_nonlinear)
  ) |>
  filter(!is.na(alpha_gi), !is.na(alpha_gf))

time_design <- joined_data |> select(curve, t_sec)

sim_one_nonlinear_curve <- function(id_val, kappa_val, sim_pars, time_design, k) {
  p <- sim_pars[sim_pars$id == id_val & sim_pars$kappa == kappa_val, ]
  ts <- time_design$t_sec[time_design$curve == id_val]
  alpha_t <- p$alpha_gf + (p$alpha_gi - p$alpha_gf) * exp(-k * ts)
  gsw_true <- gleaf_op_nonlinear(alpha_t, p$d_m2, p$s_m2, p$gec, p$x_prime, p$kappa)
  tibble(
    id = id_val,
    kappa = kappa_val,
    t_sec = ts,
    alpha = alpha_t,
    gsw_true = gsw_true,
    gsw_sim = gsw_true + rnorm(length(ts), 0, p$sigma)
  )
}

set.seed(781245)
plan(multisession, workers = 9)

grid <- sim_pars |> distinct(id, kappa)

sim_nonlinear_data <- future_map2_dfr(
  grid$id, grid$kappa,
  sim_one_nonlinear_curve,
  sim_pars = sim_pars,
  time_design = time_design,
  k = k_shared,
  .progress = TRUE,
  .options = furrr_options(seed = TRUE)
)

write_rds(sim_nonlinear_data, "objects/sim-nonlinear-aperture-data.rds")
write_rds(
  list(sim_pars = sim_pars, k_shared = k_shared,
       gec = baseline_gec, x_prime = baseline_xprime, kappa_grid = kappa_grid),
  "objects/sim-nonlinear-aperture-pars.rds"
)
