# --- Extend the aperture-kinetics model (r/43_fit-k-to-match-tau-lambda.R)
# with a mesophyll/intercellular-airspace diffusional resistance (g_ias) in
# series with the stomatal pore, following Kaiser (2009, Plant, Cell &
# Environment, 10.1111/j.1365-3040.2009.01990.x). Kaiser found that a
# pore-geometry-only model of aperture-conductance underestimates how
# steeply conductance rises with aperture, and that adding a diffusional
# resistance in series downstream of the pore substantially improved the
# fit to measured data.
#
# r/43 showed that first-order/exponential aperture kinetics, mapped
# through the pore-geometry-only model (gleaf_op()), can match any curve's
# real tau almost exactly (via k) but is essentially incapable of
# reproducing lambda != ~1 (i.e. cannot generate the accelerating,
# lambda > 1 dynamics seen in most real curves). This script asks: does
# adding a series g_ias resistance (gleaf_op_series(); r/functions.R)
# expand the range of lambda achievable by this otherwise identical
# mechanistic model?
#
# Because a stronger (lower) g_ias requires a larger fractional aperture to
# achieve the same observed gi/gf (some of the "conductance budget" is
# spent on the fixed series resistor), it also changes exactly WHERE on
# the nonlinear gleaf_op_series(alpha) curve the trajectory operates,
# which can introduce curvature (lambda != 1) that a k-only rescaling
# cannot. We repeat the per-curve k-optimization from r/43 across a small
# grid of assumed g_ias values (not directly measured), and compare
# against r/43's g_ias = Inf (i.e. no series resistance) results.

source("r/header.R")

baseline_gec <- 0.02
baseline_xprime <- 6.8e-8

# Illustrative g_ias grid (mol m^-2 s^-1); not directly measured. Chosen
# to span from a mild to a fairly strong limitation relative to the
# observed range of gi (~0.02-1.3 mol m^-2 s^-1) in this dataset.
gias_grid <- c(2, 1, 0.5)

tbl_implied_aperture <- read_rds("objects/tbl-implied-aperture.rds")
pars_summary <- read_rds("objects/pars-summary.rds")
joined_data <- read_rds("data/joined-data.rds")
# g_ias = Inf baseline (pore-geometry-only model) already computed in r/43
tbl_k_match_baseline <- read_rds("objects/tbl-k-match-tau-lambda.rds") |>
  mutate(g_ias = Inf)

real_pars <- pars_summary |>
  filter(variable %in% c("b_logtau_Intercept", "b_loglambda_Intercept")) |>
  select(id, variable, mean) |>
  pivot_wider(names_from = variable, values_from = mean) |>
  transmute(id, tau_real = exp(b_logtau_Intercept), lambda_real = exp(b_loglambda_Intercept))

sim_pars <- tbl_implied_aperture |>
  filter(gec == baseline_gec, x_prime == baseline_xprime) |>
  select(id, d_m2, s_m2, gec, x_prime,
         gi = ginit_mean, gf = gfinal_mean, gmax, amax_um) |>
  left_join(real_pars, by = "id") |>
  filter(!is.na(tau_real), !is.na(lambda_real))

time_design <- joined_data |> select(curve, t_sec)

# Deterministic Weibull fit (no noise; used purely to probe what
# (tau_hat, lambda_hat) a given k produces for a curve's anatomy)
fit_weibull_deterministic <- function(t_sec, gsw, gf0, gi0) {
  tryCatch(
    nls(
      gsw ~ gf + (gi - gf) * exp(-(t_sec / tau) ^ lambda),
      data = data.frame(t_sec = t_sec, gsw = gsw),
      start = list(gf = gf0, gi = gi0, tau = 200, lambda = 1.2),
      control = nls.control(maxiter = 200, warnOnly = TRUE)
    ),
    error = function(e) NULL
  )
}

# tau_hat, lambda_hat produced by rate constant k = exp(logk), given the
# aperture boundary conditions already implied by g_ias for this curve
tau_lambda_at_k_series <- function(logk, alpha_i, alpha_f, p, ts, g_ias) {
  k <- exp(logk)
  alpha_t <- alpha_f + (alpha_i - alpha_f) * exp(-k * ts)
  gsw_true <- gleaf_op_series(alpha_t, p$d_m2, p$s_m2, p$gec, p$x_prime, g_ias)
  m <- fit_weibull_deterministic(ts, gsw_true, gf0 = p$gf, gi0 = p$gi)
  if (is.null(m) || !m$convInfo$isConv) return(c(tau = NA_real_, lambda = NA_real_))
  cf <- coef(m)
  c(tau = unname(cf["tau"]), lambda = unname(cf["lambda"]))
}

loss_fn_series <- function(logk, alpha_i, alpha_f, p, ts, g_ias) {
  tl <- tau_lambda_at_k_series(logk, alpha_i, alpha_f, p, ts, g_ias)
  if (any(is.na(tl))) return(1e6)
  (log(tl["tau"]) - log(p$tau_real))^2 + (log(tl["lambda"]) - log(p$lambda_real))^2
}

fit_k_one_curve_series <- function(id_val, sim_pars, time_design, g_ias) {
  p <- sim_pars[sim_pars$id == id_val, ]
  ts <- time_design$t_sec[time_design$curve == id_val]

  alpha_i <- solve_alpha_series(p$gi, p$d_m2, p$s_m2, p$gec, p$x_prime, g_ias)
  alpha_f <- solve_alpha_series(p$gf, p$d_m2, p$s_m2, p$gec, p$x_prime, g_ias)

  if (is.na(alpha_i) || is.na(alpha_f)) {
    return(tibble(id = id_val, g_ias = g_ias, k_hat = NA_real_, tau_hat = NA_real_,
                  lambda_hat = NA_real_, loss = NA_real_))
  }

  opt <- tryCatch(
    optimize(loss_fn_series, interval = log(c(1e-5, 0.5)),
             alpha_i = alpha_i, alpha_f = alpha_f, p = p, ts = ts, g_ias = g_ias,
             tol = 1e-6),
    error = function(e) NULL
  )
  if (is.null(opt)) {
    return(tibble(id = id_val, g_ias = g_ias, k_hat = NA_real_, tau_hat = NA_real_,
                  lambda_hat = NA_real_, loss = NA_real_))
  }
  tl <- tau_lambda_at_k_series(opt$minimum, alpha_i, alpha_f, p, ts, g_ias)
  tibble(id = id_val, g_ias = g_ias, k_hat = exp(opt$minimum),
         tau_hat = unname(tl["tau"]), lambda_hat = unname(tl["lambda"]),
         loss = opt$objective)
}

plan(multisession, workers = 9)

grid <- expand_grid(id = sim_pars$id, g_ias = gias_grid)

k_match_series_fits <- future_map2_dfr(
  grid$id, grid$g_ias,
  fit_k_one_curve_series,
  sim_pars = sim_pars,
  time_design = time_design,
  .progress = TRUE,
  .options = furrr_options(seed = TRUE)
)

tbl_k_match_gias_series <- k_match_series_fits |>
  left_join(sim_pars, by = "id") |>
  bind_rows(
    tbl_k_match_baseline |>
      select(id, g_ias, k_hat, tau_hat, lambda_hat, loss) |>
      left_join(sim_pars, by = "id")
  ) |>
  mutate(
    fgmax = gi / gmax,
    tau_logresid = log(tau_hat) - log(tau_real),
    lambda_logresid = log(lambda_hat) - log(lambda_real)
  )

write_rds(tbl_k_match_gias_series, "objects/tbl-k-match-gias-series.rds")
