# Plot effects of fgmax on tau and lambda - WORK IN PROGRESS
source("r/header.R")

fit = read_rds("objects/best_model.rds")
dat = fit$data

df_new = dat |>
  summarize(
    min_logitfgmax = min(logitfgmax),
    max_logitfgmax = max(logitfgmax),
    .by = c(lighttreatment, lightintensity, leaftype)
  ) |>
  reframe(
    logitfgmax = seq(min_logitfgmax, max_logitfgmax, length.out = 100),
    .by = c(lighttreatment, lightintensity, leaftype)
  ) |>
  mutate(
    logtausd = 0,
    loglambdasd = 0,
    variable = glue("...{row_number()}")
  )

df_pred_tau = posterior_epred(fit,
                              newdata = df_new,
                              re_formula = NA,
                              resp = "logtaumean") |>
  as_draws_df() |>
  summarize_draws(median, quantile2, .args = list(probs = c(0.025, 0.975))) |>
  full_join(df_new, by = join_by(variable)) |>
  rename(logtaumean = median)

df_pred_lambda = posterior_epred(fit,
                                 newdata = df_new,
                                 re_formula = NA,
                                 resp = "loglambdamean") |>
  as_draws_df() |>
  summarize_draws(median, quantile2, .args = list(probs = c(0.025, 0.975))) |>
  full_join(df_new, by = join_by(variable)) |>
  rename(loglambdamean = median)

# fgmax vs. tau
gp_tau = ggplot(dat, aes(plogis(logitfgmax), exp(logtaumean), color = leaftype)) +
  geom_point(alpha = 0.5) +
  geom_ribbon(
    data = df_pred_tau,
    aes(
      x = plogis(logitfgmax),
      ymin = exp(`q2.5`),
      ymax = exp(`q97.5`),
      fill = leaftype,
      color = NULL
    ),
    alpha = 0.3
  ) +
  geom_line(data = df_pred_tau, aes(x = plogis(logitfgmax))) +
  facet_grid(lightintensity ~ lighttreatment) +
  scale_x_continuous(transform = "logit", breaks = c(0.025, 0.1, 0.4)) +
  scale_y_log10(breaks = c(50, 100, 200, 400)) +
  scale_fill_manual(values = c(col_amphi, col_pseudohypo)) +
  scale_color_manual(values = c(col_amphi, col_pseudohypo)) +
  labs(
    x = "$f_\\mathrm{gmax}$ (logit-scale)",
    y = "$\\tau$ (s, log-scale)",
    color = "Leaf type:",
    fill = "Leaf type:"
  ) +
  theme(legend.position = "bottom")

# fgmax vs. lambda
gp_lambda = ggplot(dat, aes(plogis(logitfgmax), exp(loglambdamean), color = leaftype)) +
  geom_point(alpha = 0.5) +
  geom_ribbon(
    data = df_pred_lambda,
    aes(
      x = plogis(logitfgmax),
      ymin = exp(`q2.5`),
      ymax = exp(`q97.5`),
      fill = leaftype,
      color = NULL
    ),
    alpha = 0.3
  ) +
  geom_line(data = df_pred_lambda, aes(x = plogis(logitfgmax))) +
  facet_grid(lightintensity ~ lighttreatment) +
  scale_x_continuous(transform = "logit", breaks = c(0.025, 0.1, 0.4)) +
  scale_y_log10(breaks = c(1, 2)) +
  scale_fill_manual(values = c(col_amphi, col_pseudohypo)) +
  scale_color_manual(values = c(col_amphi, col_pseudohypo)) +
  labs(
    x = "$f_\\mathrm{gmax}$ (logit-scale)",
    y = "$\\lambda$ (unitless, log-scale)",
    color = "Leaf type:",
    fill = "Leaf type:"
  ) +
  theme(legend.position = "bottom")

# Annotate and write
gp1 = annotate_figure(
  gp_tau,
  top = ggpubr::text_grob("        Growth light intensity"),
  right = ggpubr::text_grob("Measurement light intensity        ", rot = -90)
)

# ggsave(
#   filename = "figures/fgmax-tau.pdf",
#   width = 5,
#   height = 5
# )

options(
  tikzLatexPackages = c(
    getOption("tikzLatexPackages"),
    "\\usepackage{siunitx}"
  )
)

tikz(
  "figures/fgmax-tau.tex",
  standAlone = TRUE,
  width = 5,
  height = 5
)
print(gp1)
dev.off()

system("cd figures; pdflatex fgmax-tau.tex; rm fgmax-tau.aux fgmax-tau.log")

gp2 = annotate_figure(
  gp_lambda,
  top = ggpubr::text_grob("        Growth light intensity"),
  right = ggpubr::text_grob("Measurement light intensity        ", rot = -90)
)

# ggsave(
#   plot = gp_lambda,
#   filename = "figures/fgmax-lambda.pdf",
#   width = 5,
#   height = 5
# )

tikz(
  "figures/fgmax-lambda.tex",
  standAlone = TRUE,
  width = 5,
  height = 5
)
print(gp2)
dev.off()

system("cd figures; pdflatex fgmax-lambda.tex; rm fgmax-lambda.aux fgmax-lambda.log")
