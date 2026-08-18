# --- Compare the null distribution to the real, observed gi-tau association ---

source("r/header.R")

null_sim_results = read_rds("objects/null-sim-fgmax-tau.rds")
joined_summary = read_rds("data/joined-summary.rds")
real_cor = joined_summary |>
  filter(!is.na(ginit_mean), !is.na(logtau_mean)) |>
  {\(.df) cor.test(.df$ginit_mean, .df$logtau_mean)}()

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
    x = expression("Null-simulation correlation: log(" * italic(g)[i] * ") vs. log(" * tau * ")"),
    y = "Frequency"
  )

ggsave("figures/null-sim-fgmax-tau.pdf", p_null_sim, width = 4.5, height = 4)
