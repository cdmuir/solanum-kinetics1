# --- Plot effects of fgmax components (gi and gmax) on tau and lambda --------

source("r/header.R")

selected_model = read_rds("objects/selected_model.rds")
dat = selected_model$data

# gi varies; gmax fixed
df_gi = dat |>
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

# gmax varies; gi fixed
df_gmax = dat |>
  summarize(
    min_loggmax = min(loggmax),
    max_loggmax = max(loggmax),
    loggi = median(loggi),
    .by = c(lighttreatment, lightintensity, leaftype)
  ) |>
  reframe(
    loggmax = seq(min_loggmax, max_loggmax, length.out = 100),
    loggi = first(loggi),
    .by = c(lighttreatment, lightintensity, leaftype)
  ) |>
  mutate(
    logtausd = 0,
    loglambdasd = 0,
    variable = glue("...{row_number()}")
  )

# tau as gi and gmax vary
df_gi_tau = posterior_epred(
  selected_model,
  newdata = df_gi,
  re_formula = NA,
  resp = "logtaumean"
) |>
  as_draws_df() |>
  summarize_draws(median, quantile2, .args = list(probs = c(0.025, 0.975))) |>
  full_join(df_gi, by = join_by(variable)) |>
  rename(logtaumean = median)

df_gmax_tau = posterior_epred(
  selected_model,
  newdata = df_gmax,
  re_formula = NA,
  resp = "logtaumean"
) |>
  as_draws_df() |>
  summarize_draws(median, quantile2, .args = list(probs = c(0.025, 0.975))) |>
  full_join(df_gmax, by = join_by(variable)) |>
  rename(logtaumean = median)

# lambda as gi and gmax vary
df_gi_lambda = posterior_epred(
  selected_model,
  newdata = df_gi,
  re_formula = NA,
  resp = "loglambdamean"
) |>
  as_draws_df() |>
  summarize_draws(median, quantile2, .args = list(probs = c(0.025, 0.975))) |>
  full_join(df_gi, by = join_by(variable)) |>
  rename(loglambdamean = median)

df_gmax_lambda = posterior_epred(
  selected_model,
  newdata = df_gmax,
  re_formula = NA,
  resp = "loglambdamean"
) |>
  as_draws_df() |>
  summarize_draws(median, quantile2, .args = list(probs = c(0.025, 0.975))) |>
  full_join(df_gmax, by = join_by(variable)) |>
  rename(loglambdamean = median)

# gi vs. tau
gi_tau = ggplot(dat, aes(exp(loggi), exp(logtaumean), color = leaftype)) +
  geom_point(alpha = 0.5) +
  geom_ribbon(
    data = df_gi_tau,
    aes(
      x = exp(loggi),
      ymin = exp(`q2.5`),
      ymax = exp(`q97.5`),
      fill = leaftype,
      color = NULL
    ),
    alpha = 0.3
  ) +
  geom_line(data = df_gi_tau, aes(x = exp(loggi))) +
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

# gmax vs. tau
gmax_tau = ggplot(dat, aes(exp(loggmax), exp(logtaumean), color = leaftype)) +
  geom_point(alpha = 0.5) +
  # geom_ribbon(
  #   data = df_gmax_tau,
  #   aes(
  #     x = exp(loggmax),
  #     ymin = exp(`q2.5`),
  #     ymax = exp(`q97.5`),
  #     fill = leaftype,
  #     color = NULL
  #   ),
  #   alpha = 0.3
  # ) +
  # geom_line(data = df_gmax_tau, aes(x = exp(loggmax))) +
  facet_grid(lightintensity ~ lighttreatment) +
  scale_x_log10(breaks = c(0.01, 0.1, 1)) +
  scale_y_log10(breaks = c(50, 100, 200, 400)) +
  scale_fill_manual(values = c(col_amphi, col_pseudohypo)) +
  scale_color_manual(values = c(col_amphi, col_pseudohypo)) +
  labs(
    x = "$g_\\mathrm{max}$ ($\\SI{}{\\mole\\per\\meter\\squared\\per\\second}$, log-scale)",
    y = "$\\tau$ (s, log-scale)",
    color = "Leaf type:",
    fill = "Leaf type:"
  ) +
  theme(legend.position = "bottom")

# gi vs. lambda
gi_lambda = ggplot(dat, aes(exp(loggi), exp(loglambdamean), color = leaftype)) +
  geom_point(alpha = 0.5) +
  geom_ribbon(
    data = df_gi_lambda,
    aes(
      x = exp(loggi),
      ymin = exp(`q2.5`),
      ymax = exp(`q97.5`),
      fill = leaftype,
      color = NULL
    ),
    alpha = 0.3
  ) +
  geom_line(data = df_gi_lambda, aes(x = exp(loggi))) +
  facet_grid(lightintensity ~ lighttreatment) +
  scale_x_log10(breaks = c(0.01, 0.1, 1)) +
  scale_y_log10(breaks = c(1, 2)) +
  scale_fill_manual(values = c(col_amphi, col_pseudohypo)) +
  scale_color_manual(values = c(col_amphi, col_pseudohypo)) +
  labs(
    x = "$g_\\mathrm{i}$ ($\\SI{}{\\mole\\per\\meter\\squared\\per\\second}$, log-scale)",
    y = "$\\lambda$ (unitless, log-scale)",
    color = "Leaf type:",
    fill = "Leaf type:"
  ) +
  theme(legend.position = "bottom")

# gmax vs. lambda
gmax_lambda = ggplot(dat, aes(exp(loggmax), exp(loglambdamean), color = leaftype)) +
  geom_point(alpha = 0.5) +
  # geom_ribbon(
  #   data = df_gmax_lambda,
  #   aes(
  #     x = exp(loggmax),
  #     ymin = exp(`q2.5`),
  #     ymax = exp(`q97.5`),
  #     fill = leaftype,
  #     color = NULL
  #   ),
  #   alpha = 0.3
  # ) +
  # geom_line(data = df_gmax_lambda, aes(x = exp(loggmax))) +
  facet_grid(lightintensity ~ lighttreatment) +
  scale_x_log10(breaks = c(0.01, 0.1, 1)) +
  scale_y_log10(breaks = c(1, 2)) +
  scale_fill_manual(values = c(col_amphi, col_pseudohypo)) +
  scale_color_manual(values = c(col_amphi, col_pseudohypo)) +
  labs(
    x = "$g_\\mathrm{max}$ ($\\SI{}{\\mole\\per\\meter\\squared\\per\\second}$, log-scale)",
    y = "$\\lambda$ (unitless, log-scale)",
    color = "Leaf type:",
    fill = "Leaf type:"
  ) +
  theme(legend.position = "bottom")

# Annotate and write
gp1 = annotate_figure(
  gi_tau,
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
  gmax_tau,
  top = ggpubr::text_grob("        Growth light intensity"),
  right = ggpubr::text_grob("Measurement light intensity        ", rot = -90)
)

tikz(
  "figures/gmax-tau.tex",
  standAlone = TRUE,
  width = 5,
  height = 5
)
print(gp2)
dev.off()

system("cd figures; pdflatex gmax-tau.tex; rm gmax-tau.aux gmax-tau.log")

gp3 = annotate_figure(
  gi_lambda,
  top = ggpubr::text_grob("        Growth light intensity"),
  right = ggpubr::text_grob("Measurement light intensity        ", rot = -90)
)

tikz(
  "figures/gi-lambda.tex",
  standAlone = TRUE,
  width = 5,
  height = 5
)
print(gp3)
dev.off()

system("cd figures; pdflatex gi-lambda.tex; rm gi-lambda.aux gi-lambda.log")

gp4 = annotate_figure(
  gmax_lambda,
  top = ggpubr::text_grob("        Growth light intensity"),
  right = ggpubr::text_grob("Measurement light intensity        ", rot = -90)
)

tikz(
  "figures/gmax-lambda.tex",
  standAlone = TRUE,
  width = 5,
  height = 5
)
print(gp4)
dev.off()

system("cd figures; pdflatex gmax-lambda.tex; rm gmax-lambda.aux gmax-lambda.log")
