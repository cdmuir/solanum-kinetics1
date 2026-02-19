# Plot accession-level gcl against tau (a bit messy, but getting close)
source("r/header.R")

fit = read_rds("objects/best_model.rds")

df_new = crossing(
  phy = unique(fit$data$phy),
  curve_type = unique(fit$data$curve_type),
  lightintensity = unique(fit$data$lightintensity),
  lighttreatment = unique(fit$data$lighttreatment),
  logtausd = 0
) |>
  mutate(accession = phy,
         variable = paste0("...", row_number())) |> full_join(
           fit$data |>
             summarize(
               logfgmax = median(logfgmax),
               .by = c(curve_type, lightintensity, lighttreatment)
             ),
           by = join_by(curve_type, lightintensity, lighttreatment)
         )

df_pred_loggcl = posterior_epred(fit, newdata = df_new, re_formula = ~ phy, resp = "loggcl") |>
  as_draws_df() |>
  summarize_draws() |>
  full_join(df_new, by = "variable") |>
  select(accession, curve_type, lightintensity, lighttreatment, loggcl = median)

df_pred_logtau = posterior_epred(fit, newdata = df_new, re_formula = ~ phy, resp = "logtaumean") |>
  as_draws_df() |>
  summarize_draws() |>
  full_join(df_new, by = "variable") |>
  select(accession, curve_type, lightintensity, lighttreatment, logtau = median)

Sigma_median = fit |>
  as_draws_df() |>
  select(
    starts_with("."),
    sd_phy__loggcl_Intercept,
    sd_phy__logtaumean_Intercept,
    cor_phy__logtaumean_Intercept__loggcl_Intercept
  ) |>
  transmute(
    vx = sd_phy__loggcl_Intercept^2,
    vy = sd_phy__logtaumean_Intercept^2,
    cxy = cor_phy__logtaumean_Intercept__loggcl_Intercept * sd_phy__loggcl_Intercept * sd_phy__logtaumean_Intercept
  ) %>%
  summarise(across(everything(), median)) %>%
  {
    matrix(
      c(.$vx, .$cxy, .$cxy, .$vy),
      nrow = 2,
      byrow = TRUE,
      dimnames = list(c("X", "Y"), c("X", "Y"))
    )
  }

df_acc = fit$data |>
  summarize(
    loggcl = median(loggcl),
    logtaumean = median(logtaumean), 
    .by = c(phy, curve_type, lightintensity, lighttreatment)) 

df_mu = df_acc |>
  summarize(
    loggcl = mean(loggcl),
    logtaumean = mean(logtaumean), 
    .by = c(curve_type, lightintensity, lighttreatment))

df_ellipse = df_mu |>
  mutate(ell = map2(loggcl, logtaumean, \(.x, .y) {
    ell = ellipse_points(mu = c(.x, .y), Sigma_median, level = 0.95)
  })) |>
  unnest(cols = ell) |>
  select(-loggcl, -logtaumean) |>
  rename(loggcl = x, logtaumean = y)

p = ggplot(df_acc, aes(exp(loggcl), exp(logtaumean), color = curve_type)) +
  geom_polygon(
    data = df_ellipse,
    aes(exp(loggcl), exp(logtaumean), fill = curve_type),
    color = "black",
    alpha = 0.25,
    inherit.aes = FALSE
  ) +
  geom_point() +
  facet_grid(lightintensity ~ lighttreatment) +
  scale_x_log10() +
  scale_y_log10(breaks = c(100, 200, 400)) +
  scale_fill_manual(values = c("steelblue", "tomato")) +
  scale_color_manual(values = c("steelblue", "tomato")) +
  labs(
    x = expression(Guard ~ cell ~ length ~ (paste(mu, "m"))),
    y = expression(tau ~ (s)),
    color = "Curve type",
    fill = "Curve type"
  ) +
  theme(legend.position = "bottom")

ggpubr::annotate_figure(
  p,
  top = ggpubr::text_grob("        Growth light intensity"),
  right = ggpubr::text_grob("Measurement light intensity        ", rot = -90)
)

ggsave("figures/gcl-tau.pdf", width = 5, height = 5)
