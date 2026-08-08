# Reviewer comment R1.6: the fgmax path analysis (r/22_plot-mediation.R) is
# not a randomized-treatment causal mediation analysis, and its "sequential
# ignorability" assumption (no unmeasured confounder of the fgmax-tau
# relationship) may not hold. This script implements an Imai, Keele &
# Yamamoto (2010)-style sensitivity analysis for a linear mediation model,
# parameterized by rho = Cor(e_M, e_Y), the residual correlation between the
# mediator (fgmax) and outcome (tau) equations.
#
# Key insight: the selected joint multivariate model (objects/best_model.rds)
# already estimates this residual correlation freely
# (rescor__logtaumean__loggi, from set_rescor(TRUE) in r/10_fit-all.R).
# Because the fgmax equation's predictors (lighttreatment + lightintensity +
# leaftype + the same 3 group-level terms) are a strict subset of tau's
# (... + loggcl + loggi + the same group-level terms), the classic
# Zellner SUR-equivalence result means jointly estimating this correlation
# does not change the point estimate of b_logtaumean_loggi relative to
# a naive model that assumes rho = 0 (sequential ignorability). This lets us
# treat b_logtaumean_loggi as the "naive" beta2 in the sensitivity
# framework, and rescor__logtaumean__loggi as a genuine, data-driven
# estimate of rho -- not just a hypothetical value to sweep over.
#
# For a linear mediator model M = ... + e_M (Var = sigma_M^2) and outcome
# model Y = ... + beta2*M + e_Y (Var = sigma_Y^2) with Cor(e_M, e_Y) = rho,
# the bias-adjusted mediator effect is:
#   beta2(rho) = beta2_hat - rho * sigma_Y / sigma_M
# so ACME(rho) = gamma1 * beta2(rho), and the breakdown point (ACME = 0) is:
#   rho* = beta2_hat * sigma_M / sigma_Y
# (derived from the standard omitted-variable-bias formula for endogenous
# regressors under joint normality; verified against a small hand derivation
# -- see plan notes -- rather than a package implementation, since our model
# is not a simple two-equation OLS setup that `mediation::medsens()` expects.)

source("r/header.R")

selected_model = read_rds("objects/selected_model.rds")

draws = as_draws_df(selected_model, variable = c(
  "b_logtaumean_loggi",
  "b_loggi_lighttreatmentsun",
  "b_loggi_lightintensityhigh",
  "b_loggi_leaftypepseudohypo",
  "b_logtaumean_lighttreatmentsun",
  "b_logtaumean_lightintensityhigh",
  "b_logtaumean_leaftypepseudohypo",
  "sigma_loggi",
  "sigma_logtaumean",
  "rescor__logtaumean__loggi"
))

beta2 = draws$b_logtaumean_loggi
sigma_M = draws$sigma_loggi
sigma_Y = draws$sigma_logtaumean
rho_hat = draws$rescor__logtaumean__loggi

# Breakdown point: value of rho at which ACME = 0. Identical across
# treatments because ACME(rho) = gamma1 * beta2(rho), and gamma1 does not
# depend on rho, so the zero-crossing depends only on beta2(rho).
rho_star = beta2 * sigma_M / sigma_Y

rho_summary = tibble(
  parameter = c("rho_hat (estimated)", "rho_star (breakdown)"),
  estimate = c(median(rho_hat), median(rho_star)),
  lower = c(quantile(rho_hat, 0.025), quantile(rho_star, 0.025)),
  upper = c(quantile(rho_hat, 0.975), quantile(rho_star, 0.975))
)

# Posterior probability that the model's own estimated residual correlation
# is smaller in magnitude than the breakdown point, i.e. that taking rho_hat
# at face value as confounding-relevant correlation would not be enough to
# nullify the mediated effect.
prob_below_breakdown = mean(rho_hat < rho_star)

write_rds(
  list(rho_summary = rho_summary, prob_below_breakdown = prob_below_breakdown),
  "objects/mediation-sensitivity-rhostar.rds"
)

message(glue::glue(
  "rho_hat = {formatC(median(rho_hat), format = 'f', digits = 2)} ",
  "[{formatC(quantile(rho_hat, 0.025), format = 'f', digits = 2)}, ",
  "{formatC(quantile(rho_hat, 0.975), format = 'f', digits = 2)}]; ",
  "rho* = {formatC(median(rho_star), format = 'f', digits = 2)} ",
  "[{formatC(quantile(rho_star, 0.025), format = 'f', digits = 2)}, ",
  "{formatC(quantile(rho_star, 0.975), format = 'f', digits = 2)}]; ",
  "P(rho_hat < rho*) = {formatC(prob_below_breakdown, format = 'f', digits = 2)}"
))

# --- Sensitivity curves: ACME and proportion mediated vs. rho -------------

rho_grid = seq(-0.9, 0.9, by = 0.02)

treatments = list(
  `Growth light (sun)` = list(
    gamma1 = draws$b_loggi_lighttreatmentsun,
    direct = draws$b_logtaumean_lighttreatmentsun
  ),
  `Measurement light (high)` = list(
    gamma1 = draws$b_loggi_lightintensityhigh,
    direct = draws$b_logtaumean_lightintensityhigh
  ),
  `Leaf type (pseudohypo)` = list(
    gamma1 = draws$b_loggi_leaftypepseudohypo,
    direct = draws$b_logtaumean_leaftypepseudohypo
  )
)

sens_curve = map_dfr(names(treatments), function(tx) {
  gamma1 = treatments[[tx]]$gamma1
  direct = treatments[[tx]]$direct
  map_dfr(rho_grid, function(rho) {
    beta2_rho = beta2 - rho * sigma_Y / sigma_M
    acme = gamma1 * beta2_rho
    prop_mediated = acme / (direct + acme)
    tibble(
      treatment = tx,
      rho = rho,
      acme_median = median(acme),
      acme_lower = quantile(acme, 0.025),
      acme_upper = quantile(acme, 0.975),
      prop_median = median(prop_mediated),
      prop_lower = quantile(prop_mediated, 0.025),
      prop_upper = quantile(prop_mediated, 0.975)
    )
  })
})

write_rds(sens_curve, "objects/mediation-sensitivity-curve.rds")

# --- Figures ---------------------------------------------------------------

rho_star_ci = quantile(rho_star, c(0.025, 0.5, 0.975))
rho_hat_ci = quantile(rho_hat, c(0.025, 0.5, 0.975))

p_sens = ggplot(sens_curve, aes(rho, acme_median, color = treatment, fill = treatment)) +
  geom_ribbon(aes(ymin = acme_lower, ymax = acme_upper), alpha = 0.15, color = NA) +
  geom_line(linewidth = 0.9) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
  annotate("rect", xmin = rho_hat_ci[1], xmax = rho_hat_ci[3],
           ymin = -Inf, ymax = Inf, fill = "grey60", alpha = 0.3) +
  geom_vline(xintercept = rho_hat_ci[2], linetype = "dotted", color = "black") +
  geom_vline(xintercept = rho_star_ci[2], linetype = "solid", color = "black") +
  labs(
    x = expression("Assumed residual correlation, " * rho),
    y = "Indirect effect on log(tau)\n(through gi)",
    color = "Treatment", fill = "Treatment"
  ) +
  guides(color = guide_legend(nrow = 2), fill = guide_legend(nrow = 2)) +
  theme(legend.position = "bottom")

ggsave("figures/mediation-sensitivity.pdf", p_sens, width = 7, height = 5.5)

rho_compare = bind_rows(
  tibble(value = rho_hat, type = "Estimated (rescor in model)"),
  tibble(value = rho_star, type = "Breakdown point (rho*)")
)

p_rho_compare = ggplot(rho_compare, aes(value, fill = type)) +
  geom_density(alpha = 0.5, color = NA) +
  labs(x = expression(rho), y = "Posterior density", fill = NULL) +
  theme(legend.position = "bottom")

ggsave("figures/mediation-sensitivity-rho-compare.pdf", p_rho_compare, width = 6, height = 4)
