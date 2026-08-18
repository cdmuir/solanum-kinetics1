# --- Compare posterior estimates with and without the VPD covariate ----------
#
# Compares key posterior estimates between the original model
# (objects/selected_model.rds) and the model refit with final VPD as a curve-
# level covariate of tau and lambda (objects/selected_model_vpd.rds; see
# r/30_refit-vpd.R), to check whether the gi -> tau effect and the gcl-tau
# phylogenetic correlation are robust to accounting for realized VPD exposure
# (Notes S1).

source("r/header.R")

selected_model = read_rds("objects/selected_model.rds")
selected_model_vpd = read_rds("objects/selected_model_vpd.rds")

pars_of_interest = c(
  "b_logtaumean_loggi",
  "cor_phy__logtaumean_Intercept__loggcl_Intercept"
)

par_labels = c(
  b_logtaumean_loggi = "log(italic(g)[i])~on~log(tau)",
  cor_phy__logtaumean_Intercept__loggcl_Intercept = "phylogenetic~corr(log(italic(l)[gc])*','~log(tau))"
)

comp_pars = bind_rows(
  extract_pars(selected_model, "Original"),
  extract_pars(selected_model_vpd, "VPD-adjusted")
) |>
  mutate(
    parameter = factor(
      par_labels[parameter],
      levels = par_labels
    ),
    model = factor(model, levels = c("Original", "VPD-adjusted"))
  )

write_rds(comp_pars, "objects/tbl-comparison-vpd.rds")

# The fixed effect of final_vpd itself on log(tau) (only defined in the
# VPD-adjusted model, so it does not fit the two-model comparison above).
final_vpd_effect = as_draws_df(selected_model_vpd, variable = "b_logtaumean_final_vpd") |>
  as_tibble() |>
  summarize(
    estimate = mean(b_logtaumean_final_vpd),
    lower = quantile(b_logtaumean_final_vpd, 0.025),
    upper = quantile(b_logtaumean_final_vpd, 0.975)
  )

write_rds(final_vpd_effect, "objects/tbl-final-vpd-effect.rds")

p_comparison = ggplot(comp_pars, aes(x = model, y = estimate, color = model)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_pointrange(aes(ymin = lower, ymax = upper), linewidth = 0.9) +
  facet_wrap(~ parameter, scales = "free_y", labeller = label_parsed) +
  scale_color_manual(values = c(Original = "grey30", `VPD-adjusted` = col_amphi)) +
  labs(x = NULL, y = "Posterior estimate (95% CI)", color = NULL) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    strip.text = element_text(size = 10)
  )

ggsave("figures/vpd-comparison.pdf", p_comparison, width = 9, height = 4)
