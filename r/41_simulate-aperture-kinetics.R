# --- Simulate gsw(t) from a simple, anatomy-independent aperture kinetics
# model, to ask: if every curve shared the SAME underlying pore-closure
# rate constant k (i.e. aperture kinetics have no true dependence on gi,
# gmax, or f_gmax), would fitting the phenomenological Weibull model to
# the resulting gsw(t) still produce fitted tau/lambda that covary with
# gi and gmax -- purely as an artifact of the nonlinear, saturating
# relationship between fractional aperture and conductance (Ochoa et al.
# 2024 Eq. 3)? This is a mechanistic complement to the curve-fitting null
# simulation in Notes S2 (r/34_simulate-null.R): that script asks whether
# estimating gi and tau from the same nonlinear fit can spuriously
# generate a gi-tau association; this script asks whether the anatomical
# aperture-to-conductance mapping itself, on its own, would predict
# tau/lambda to covary with gi and gmax even absent any true difference
# in intrinsic closure rate.
#
# Model:
#   - Aperture kinetics: da/dt = k * (a_final - a), i.e. a(t) = a_final +
#     (a_initial - a_final) * exp(-k * t). Since alpha = a / p is just a
#     rescaling of a, this is equivalent to d(alpha)/dt = k * (alpha_final
#     - alpha), so we work directly in alpha (fractional aperture) space.
#   - "All stomata on a leaf behave homogeneously": every pore has the
#     same alpha(t) at every instant, so leaf-level conductance is simply
#     gleaf_op(alpha(t), ...) (Ochoa et al. 2024 Eq. 3; r/functions.R).
#   - alpha_initial and alpha_final are calibrated (via r/40_calc-implied-
#     aperture.R) so that gleaf_op(alpha_initial) = observed gi and
#     gleaf_op(alpha_final) = observed gf for that curve exactly. This
#     preserves each curve's real gi, gf, and anatomy (d, s), while
#     imposing a single common kinetic rate constant k across all curves.
#   - k is fixed at 1 / median(real tau), so the simulated timescale is
#     broadly comparable to the real data; this is an arbitrary choice of
#     absolute timescale; testing the same qualitative question with a
#     different k should give qualitatively similar covariation patterns.
#   - Noise: i.i.d. Gaussian with each curve's own real residual SD
#     (`sigma`, from the Bayesian curve fits; objects/pars-summary.rds),
#     matching real per-curve noise magnitude and hence realistic
#     estimability of tau/lambda on refitting.
#   - gec and x_prime are fixed at one representative combination from the
#     grid explored in r/40_calc-implied-aperture.R (gec = 0.02 mol m^-2
#     s^-1, x_prime = 6.8e-8 m); since implied aperture was found to be
#     fairly insensitive to gec and only moderately sensitive to x_prime,
#     this is a reasonable default, but results should be checked against
#     other combinations if this analysis is reported.

source("r/header.R")

baseline_gec <- 0.02
baseline_xprime <- 6.8e-8

tbl_implied_aperture <- read_rds("objects/tbl-implied-aperture.rds")
pars_summary <- read_rds("objects/pars-summary.rds")
joined_data <- read_rds("data/joined-data.rds")

# Common kinetic rate constant k, shared across all curves by assumption
median_tau <- pars_summary |>
  filter(variable == "b_logtau_Intercept") |>
  summarise(median_tau = median(exp(mean))) |>
  pull(median_tau)

k_shared <- 1 / median_tau

sim_pars <- tbl_implied_aperture |>
  filter(gec == baseline_gec, x_prime == baseline_xprime) |>
  select(id, alpha_gi, alpha_gf, d_m2, s_m2, gec, x_prime,
         gi = ginit_mean, gf = gfinal_mean, gmax, amax_um) |>
  left_join(
    pars_summary |> filter(variable == "sigma") |> select(id, sigma = mean),
    by = "id"
  ) |>
  filter(!is.na(sigma))

time_design <- joined_data |> select(curve, t_sec)

sim_one_aperture_curve <- function(id_val, sim_pars, time_design, k) {
  p <- sim_pars[sim_pars$id == id_val, ]
  ts <- time_design$t_sec[time_design$curve == id_val]
  alpha_t <- p$alpha_gf + (p$alpha_gi - p$alpha_gf) * exp(-k * ts)
  gsw_true <- gleaf_op(alpha_t, p$d_m2, p$s_m2, p$gec, p$x_prime)
  tibble(
    id = id_val,
    t_sec = ts,
    alpha = alpha_t,
    gsw_true = gsw_true,
    gsw_sim = gsw_true + rnorm(length(ts), 0, p$sigma)
  )
}

set.seed(596184)
plan(multisession, workers = 9)

sim_aperture_data <- future_map_dfr(
  sim_pars$id,
  sim_one_aperture_curve,
  sim_pars = sim_pars,
  time_design = time_design,
  k = k_shared,
  .progress = TRUE,
  .options = furrr_options(seed = TRUE)
)

write_rds(sim_aperture_data, "objects/sim-aperture-data.rds")
write_rds(
  list(sim_pars = sim_pars, k_shared = k_shared,
       gec = baseline_gec, x_prime = baseline_xprime),
  "objects/sim-aperture-pars.rds"
)
