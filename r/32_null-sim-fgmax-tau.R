# Reviewer comment R1.3: gi (and hence fgmax = gi/gmax) is estimated from the
# same nonlinear Weibull curve fit that also produces tau and lambda, raising
# the concern that parameter covariance within the fit alone could induce a
# spurious gi-tau (fgmax-tau) association even with no true biological
# relationship. This script runs a null simulation to test that directly:
#
#   1. Simulate curves where the *true* gi and gf are fixed at their real
#      median values (so true fgmax has zero variance -- by construction it
#      cannot be truly associated with tau), while tau and lambda retain
#      realistic variability (matched 1:1 to real curves), using each real
#      curve's own measurement time points and residual noise magnitude.
#   2. Re-fit the same functional form (nls, as a fast approximation to the
#      production Bayesian fit -- see r/32b_validate-nls-vs-bayes.R for a
#      validation of this approximation) to each simulated curve to recover
#      gi_hat, gf_hat, tau_hat, lambda_hat.
#   3. Compute fgmax_hat = gi_hat / gmax_fixed and test its correlation with
#      tau_hat (on the same logit/log scales used in the real brms model).
#   4. Repeat over many replicate simulated datasets to build a null
#      distribution of this correlation, for comparison against the real,
#      observed fgmax-tau association.

source("r/header.R")

n_replicates = 200
sim_seed = 4471

pars_summary = read_rds("objects/pars-summary.rds")
joined_data = read_rds("data/joined-data.rds")

# --- Real per-curve building blocks: tau, lambda, sigma (residual SD), gi, gf

real_pars = pars_summary |>
  filter(variable %in% c("b_logtau_Intercept", "b_loglambda_Intercept", "sigma", "ginit", "gfinal")) |>
  select(id, variable, mean) |>
  pivot_wider(names_from = variable, values_from = mean) |>
  transmute(
    id,
    tau = exp(b_logtau_Intercept),
    lambda = exp(b_loglambda_Intercept),
    sigma,
    gi = ginit,
    gf = gfinal
  )

time_design = joined_data |> select(curve, t_sec)

# True gi/gf are fixed at these real median values for every simulated curve.
median_gi = median(real_pars$gi, na.rm = TRUE)
median_gf = median(real_pars$gf, na.rm = TRUE)

# Anatomical gmax is not part of the curve-fitting procedure (it comes from
# stomatal density/geometry, not the Weibull fit), so its value cannot affect
# whether curve-fitting induces a gi-tau artifact. We fix it here at the same
# value as median_gi purely so fgmax_hat is on a realistic, interpretable
# 0-1-ish scale; the correlation between fgmax_hat and tau_hat is invariant
# to this choice (dividing gi_hat by any positive constant just rescales it).
gmax_fixed = median_gi

# --- Simulate one dataset: real time design/tau/lambda/sigma per curve,
# true gi/gf replaced with the fixed medians -----------------------------

sim_one_dataset = function(real_pars, time_design, median_gi, median_gf) {
  real_pars |>
    left_join(time_design, by = c("id" = "curve")) |>
    mutate(
      gsw_true = median_gf + (median_gi - median_gf) * exp(-(t_sec / tau) ^ lambda),
      gsw_sim = gsw_true + rnorm(n(), 0, sigma)
    )
}

# --- Fast re-fit of the same functional form via nls() -------------------

fit_nls_one = function(df) {
  gf0 = min(df$gsw_sim)
  gi0 = max(df$gsw_sim)
  m = tryCatch(
    nls(
      gsw_sim ~ gf + (gi - gf) * exp(-(t_sec / tau) ^ lambda),
      data = df,
      start = list(gf = gf0, gi = gi0, tau = 200, lambda = 1.2),
      control = nls.control(maxiter = 200, warnOnly = TRUE)
    ),
    error = function(e) NULL
  )
  if (is.null(m) || !isTRUE(summary(m)$convInfo$isConv)) {
    return(tibble(converged = FALSE, gi_hat = NA_real_, gf_hat = NA_real_,
                  tau_hat = NA_real_, lambda_hat = NA_real_))
  }
  cf = coef(m)
  tibble(
    converged = TRUE,
    gi_hat = cf[["gi"]],
    gf_hat = cf[["gf"]],
    tau_hat = cf[["tau"]],
    lambda_hat = cf[["lambda"]]
  )
}

# --- One full replicate: simulate + refit all curves + compute correlation

run_one_replicate = function(replicate_id, real_pars, time_design, median_gi, median_gf, gmax_fixed) {
  sim_data = sim_one_dataset(real_pars, time_design, median_gi, median_gf)

  fits = sim_data |>
    group_by(id) |>
    group_modify(~ fit_nls_one(.x)) |>
    ungroup()

  fits_valid = fits |>
    filter(converged, gi_hat > 0, gi_hat < 1, tau_hat > 0) |>
    mutate(
      fgmax_hat = gi_hat / gmax_fixed,
      logitfgmax_hat = qlogis(pmin(pmax(fgmax_hat, 1e-6), 1 - 1e-6)),
      logtau_hat = log(tau_hat)
    )

  ct = suppressWarnings(cor.test(fits_valid$logitfgmax_hat, fits_valid$logtau_hat))

  tibble(
    replicate = replicate_id,
    n_curves = nrow(fits_valid),
    n_converged = sum(fits$converged),
    n_total = nrow(fits),
    cor = unname(ct$estimate),
    cor_lower = ct$conf.int[1],
    cor_upper = ct$conf.int[2]
  )
}

# --- Run all replicates ----------------------------------------------------

set.seed(sim_seed)

null_sim_results = map_dfr(seq_len(n_replicates), \(b) {
  run_one_replicate(b, real_pars, time_design, median_gi, median_gf, gmax_fixed)
}, .progress = TRUE)

write_rds(null_sim_results, "objects/null-sim-fgmax-tau.rds")

# Save one example simulated dataset (true + estimated parameters) alongside
# its fits, for diagnostic/illustrative plotting.
set.seed(sim_seed)
example_sim_data = sim_one_dataset(real_pars, time_design, median_gi, median_gf)
example_fits = example_sim_data |>
  distinct(id, tau, lambda, sigma) |>
  left_join(
    example_sim_data |>
      group_by(id) |>
      group_modify(~ fit_nls_one(.x)) |>
      ungroup(),
    by = "id"
  )
write_rds(example_fits, "objects/null-sim-fgmax-tau-example.rds")

# --- Compare to the real, observed fgmax-tau association -----------------

joined_summary = read_rds("data/joined-summary.rds")
real_cor = joined_summary |>
  filter(!is.na(f_gmax), !is.na(logtau_mean)) |>
  transmute(
    logitfgmax = qlogis(pmin(pmax(f_gmax, 1e-6), 1 - 1e-6)),
    logtau_mean
  ) |>
  {\(.df) cor.test(.df$logitfgmax, .df$logtau_mean)}()

real_cor_summary = tibble(
  cor = unname(real_cor$estimate),
  cor_lower = real_cor$conf.int[1],
  cor_upper = real_cor$conf.int[2]
)
write_rds(real_cor_summary, "objects/real-cor-fgmax-tau.rds")

# Empirical p-value: proportion of null-replicate correlations at least as
# extreme (same sign, larger magnitude) as the observed real correlation.
empirical_p = mean(abs(null_sim_results$cor) >= abs(real_cor_summary$cor))
write_rds(empirical_p, "objects/null-sim-empirical-p.rds")

message(glue::glue(
  "Null simulation: {n_replicates} replicates, median correlation = ",
  "{formatC(median(null_sim_results$cor), format = 'f', digits = 3)} ",
  "(range {formatC(min(null_sim_results$cor), format = 'f', digits = 3)} to ",
  "{formatC(max(null_sim_results$cor), format = 'f', digits = 3)}). ",
  "Observed real correlation = {formatC(real_cor_summary$cor, format = 'f', digits = 3)}. ",
  "Empirical p = {formatC(empirical_p, format = 'f', digits = 4)}."
))

# --- Figure: null distribution vs. the observed real correlation ---------

p_null_sim = ggplot(null_sim_results, aes(cor)) +
  geom_histogram(bins = 30) +
  geom_vline(xintercept = real_cor_summary$cor, color = "firebrick", linewidth = 1) +
  annotate(
    "text",
    x = real_cor_summary$cor - 0.02,
    y = Inf,
    label = "observed\n(real data)",
    color = "firebrick",
    hjust = 1,
    vjust = 1.3
  ) +
  scale_x_continuous(limits = range(c(null_sim_results$cor, real_cor_summary$cor)) + c(-0.02, 0.02)) +
  labs(
    x = expression("Null-simulation correlation: logit(" * f[gmax] * ") vs. log(" * tau * ")"),
    y = "Number of replicates"
  )

ggsave("figures/null-sim-fgmax-tau.pdf", p_null_sim, width = 6, height = 4)
