# --- For each curve, find the single rate constant k (in da/dt = k *
# (a_final - a), i.e. first-order/exponential closure in fractional
# aperture) that makes the resulting Weibull fit (tau_hat, lambda_hat) as
# close as possible -- jointly, on the log scale -- to that curve's real
# fitted tau and lambda.
#
# This differs from r/41_simulate-aperture-kinetics.R, which imposed a
# single k shared across ALL curves to test whether anatomy alone (via
# the nonlinear aperture-conductance mapping) predicts covariation of
# tau_hat/lambda_hat with gi and gmax. Here we instead ask, per curve: is
# there ANY k for which this simple mechanistic model reproduces the
# observed (tau, lambda) pair? Since k only rescales time, and alpha(t) is
# evaluated at each curve's fixed real sampling times (not rescaled with
# k), tau_hat and lambda_hat are not fully separable functions of k in
# general -- but empirically (see below), lambda_hat here turns out to be
# very insensitive to k for a given curve's anatomy, while tau_hat can be
# tuned almost freely via k. This means this simple model can generally
# match a curve's real tau almost exactly (by choosing k appropriately),
# but often cannot also match its real lambda if lambda_real is far from
# the value intrinsic to that curve's anatomy under first-order aperture
# kinetics (empirically close to lambda ~ 1, i.e. no acceleration).
#
# For each curve we use noise-free (deterministic) simulated gsw(t), since
# the goal is to characterize what this mechanistic family can achieve in
# principle, not to re-test estimability under noise (see r/41-r/42 for
# that).

source("r/header.R")

baseline_gec <- 0.02
baseline_xprime <- 6.8e-8

tbl_implied_aperture <- read_rds("objects/tbl-implied-aperture.rds")
pars_summary <- read_rds("objects/pars-summary.rds")
joined_data <- read_rds("data/joined-data.rds")

real_pars <- pars_summary |>
  filter(variable %in% c("b_logtau_Intercept", "b_loglambda_Intercept")) |>
  select(id, variable, mean) |>
  pivot_wider(names_from = variable, values_from = mean) |>
  transmute(id, tau_real = exp(b_logtau_Intercept), lambda_real = exp(b_loglambda_Intercept))

sim_pars <- tbl_implied_aperture |>
  filter(gec == baseline_gec, x_prime == baseline_xprime) |>
  select(id, alpha_gi, alpha_gf, d_m2, s_m2, gec, x_prime,
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

# tau_hat, lambda_hat produced by rate constant k = exp(logk) for one curve
tau_lambda_at_k <- function(logk, p, ts) {
  k <- exp(logk)
  alpha_t <- p$alpha_gf + (p$alpha_gi - p$alpha_gf) * exp(-k * ts)
  gsw_true <- gleaf_op(alpha_t, p$d_m2, p$s_m2, p$gec, p$x_prime)
  m <- fit_weibull_deterministic(ts, gsw_true, gf0 = p$gf, gi0 = p$gi)
  if (is.null(m) || !m$convInfo$isConv) return(c(tau = NA_real_, lambda = NA_real_))
  cf <- coef(m)
  c(tau = unname(cf["tau"]), lambda = unname(cf["lambda"]))
}

# Joint squared log-error to (tau_real, lambda_real); large penalty on
# non-convergence so optimize() avoids degenerate k values.
loss_fn <- function(logk, p, ts) {
  tl <- tau_lambda_at_k(logk, p, ts)
  if (any(is.na(tl))) return(1e6)
  (log(tl["tau"]) - log(p$tau_real))^2 + (log(tl["lambda"]) - log(p$lambda_real))^2
}

fit_k_one_curve <- function(id_val, sim_pars, time_design) {
  p <- sim_pars[sim_pars$id == id_val, ]
  ts <- time_design$t_sec[time_design$curve == id_val]

  opt <- tryCatch(
    optimize(loss_fn, interval = log(c(1e-5, 0.5)), p = p, ts = ts, tol = 1e-6),
    error = function(e) NULL
  )
  if (is.null(opt)) {
    return(tibble(id = id_val, k_hat = NA_real_, tau_hat = NA_real_,
                  lambda_hat = NA_real_, loss = NA_real_))
  }
  k_hat <- exp(opt$minimum)
  tl <- tau_lambda_at_k(opt$minimum, p, ts)
  tibble(id = id_val, k_hat = k_hat, tau_hat = unname(tl["tau"]),
         lambda_hat = unname(tl["lambda"]), loss = opt$objective)
}

plan(multisession, workers = 9)

k_match_fits <- future_map_dfr(
  sim_pars$id,
  fit_k_one_curve,
  sim_pars = sim_pars,
  time_design = time_design,
  .progress = TRUE,
  .options = furrr_options(seed = TRUE)
)

tbl_k_match_tau_lambda <- k_match_fits |>
  left_join(sim_pars, by = "id") |>
  mutate(
    fgmax = gi / gmax,
    tau_logresid = log(tau_hat) - log(tau_real),
    lambda_logresid = log(lambda_hat) - log(lambda_real)
  )

write_rds(tbl_k_match_tau_lambda, "objects/tbl-k-match-tau-lambda.rds")
