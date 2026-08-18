# --- Sensitivity analysis for the sequential-ignorability assumption in the gi-tau path ---
#
# Tests whether the gi -> tau path analysis is robust to an unmeasured
# confounder of the gi-tau relationship, following an Imai, Keele & Yamamoto
# (2010)-style sensitivity analysis for a linear mediation model, parameterized
# by rho = Cor(e_M, e_Y), the residual correlation between the mediator (gi)
# and outcome (tau) equations. The selected joint multivariate model
# (objects/selected_model.rds) already estimates this residual correlation
# (rescor__logtaumean__loggi, from set_rescor(TRUE) in r/10_fit-all.R).
#
# For a linear mediator model M = ... + e_M (Var = sigma_M^2) and outcome model
# Y = ... + beta2*M + e_Y (Var = sigma_Y^2) with Cor(e_M, e_Y) = rho, the bias-
# adjusted mediator effect is beta2(rho) = beta2_hat - rho * sigma_Y / sigma_M,
# so ACME(rho) = gamma1 * beta2(rho), and the breakdown point (ACME = 0) is
# rho* = beta2_hat * sigma_M / sigma_Y.

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
    x = "Assumed residual correlation, $\\rho$",
    y = "Indirect effect on $\\log \\left( \\tau \\right)$ through $\\log \\left( g_\\mathrm{i} \\right)$",
    color = "Treatment", fill = "Treatment"
  ) +
  guides(color = guide_legend(nrow = 2), fill = guide_legend(nrow = 2)) +
  theme(legend.position = "bottom")

tikz(
  "figures/mediation-sensitivity.tex",
  standAlone = TRUE,
  width = 7,
  height = 5.5
)
print(p_sens)
dev.off()

system("cd figures; pdflatex mediation-sensitivity.tex; rm mediation-sensitivity.aux mediation-sensitivity.log")
