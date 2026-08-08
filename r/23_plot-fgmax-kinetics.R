# Plot effects of gi on tau and lambda
source("r/header.R")

selected_model = read_rds("objects/selected_model.rds")
dat = selected_model$data

df_new = dat |>
  summarize(
    min_loggi = min(loggi),
    max_loggi = max(loggi),
    loggmax = median(loggmax),
    .by = c(lighttreatment, lightintensity, leaftype)
  ) |>
  reframe(
    loggi = seq(min_loggi, max_loggi, length.out = 100),
    loggmax = first(loggmax),
    .by = c(lighttreatment, lightintensity, leaftype)
  ) |>
  mutate(
    logtausd = 0,
    loglambdasd = 0,
    variable = glue("...{row_number()}")
  )

df_pred_tau = posterior_epred(selected_model,
                              newdata = df_new,
                              re_formula = NA,
                              resp = "logtaumean") |>
  as_draws_df() |>
  summarize_draws(median, quantile2, .args = list(probs = c(0.025, 0.975))) |>
  full_join(df_new, by = join_by(variable)) |>
  rename(logtaumean = median)

df_pred_lambda = posterior_epred(selected_model,
                                 newdata = df_new,
                                 re_formula = NA,
                                 resp = "loglambdamean") |>
  as_draws_df() |>
  summarize_draws(median, quantile2, .args = list(probs = c(0.025, 0.975))) |>
  full_join(df_new, by = join_by(variable)) |>
  rename(loglambdamean = median)

# gi vs. tau
gp_tau = ggplot(dat, aes(exp(loggi), exp(logtaumean), color = leaftype)) +
  geom_point(alpha = 0.5) +
  geom_ribbon(
    data = df_pred_tau,
    aes(
      x = exp(loggi),
      ymin = exp(`q2.5`),
      ymax = exp(`q97.5`),
      fill = leaftype,
      color = NULL
    ),
    alpha = 0.3
  ) +
  geom_line(data = df_pred_tau, aes(x = exp(loggi))) +
  facet_grid(lightintensity ~ lighttreatment) +
  scale_x_log10(breaks = c(0.01, 0.1, 1)) +
  scale_y_log10(breaks = c(50, 100, 200, 400)) +
  scale_fill_manual(values = c(col_amphi, col_pseudohypo)) +
  scale_color_manual(values = c(col_amphi, col_pseudohypo)) +
  labs(
    x = "$g_\\mathrm{i}$ ($\\SI{}{\\mole\\per\\meter\\squared\\per\\second}$, log-scale)",
    y = "$\\tau$ (s, log-scale)",
    color = "Leaf type:",
    fill = "Leaf type:"
  ) +
  theme(legend.position = "bottom")

# gi vs. lambda
gp_lambda = ggplot(dat, aes(exp(loggi), exp(loglambdamean), color = leaftype)) +
  geom_point(alpha = 0.5) +
  geom_ribbon(
    data = df_pred_lambda,
    aes(
      x = exp(loggi),
      ymin = exp(`q2.5`),
      ymax = exp(`q97.5`),
      fill = leaftype,
      color = NULL
    ),
    alpha = 0.3
  ) +
  geom_line(data = df_pred_lambda, aes(x = exp(loggi))) +
  facet_grid(lightintensity ~ lighttreatment) +
  scale_x_log10(breaks = c(0.01, 0.1, 1)) +
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

options(
  tikzLatexPackages = c(
    getOption("tikzLatexPackages"),
    "\\usepackage{siunitx}"
  )
)

tikz(
  "figures/gi-tau.tex",
  standAlone = TRUE,
  width = 5,
  height = 5
)
print(gp1)
dev.off()

system("cd figures; pdflatex gi-tau.tex; rm gi-tau.aux gi-tau.log")

gp2 = annotate_figure(
  gp_lambda,
  top = ggpubr::text_grob("        Growth light intensity"),
  right = ggpubr::text_grob("Measurement light intensity        ", rot = -90)
)

tikz(
  "figures/gi-lambda.tex",
  standAlone = TRUE,
  width = 5,
  height = 5
)
print(gp2)
dev.off()

system("cd figures; pdflatex gi-lambda.tex; rm gi-lambda.aux gi-lambda.log")

