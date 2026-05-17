# Simplified fgmax vs. tau plot: amphi, shade treatment, high light intensity
source("r/header.R")

annotate <- ggplot2::annotate

fit <- read_rds("objects/best_model.rds")
dat <- fit$data

# Summarize mean fgmax and tau by accession (on transformed scale, then back-transform)
df_acc <- dat |>
  filter(leaftype == "amphi", lighttreatment == "shade", lightintensity == "high") |>
  summarize(
    fgmax = plogis(mean(logitfgmax)),
    tau   = exp(mean(logtaumean)),
    .by   = accession
  )

ggplot(df_acc, aes(fgmax, tau)) +
  geom_smooth(method = "lm", color = col_amphi, fill = col_amphi, alpha = 0.2) +
  geom_point(color = col_amphi) +
  scale_x_continuous(
    transform = "logit",
    breaks    = c(0.025, 0.05, 0.1, 0.2, 0.4),
    labels    = scales::label_number(accuracy = 0.001)
  ) +
  scale_y_log10(breaks = c(50, 100, 200, 400)) +
  labs(
    x = expression(f[gmax]~"(logit scale)"),
    y = expression(tau~"(s, log scale)")
  ) +
  theme_classic(base_size = 12)
