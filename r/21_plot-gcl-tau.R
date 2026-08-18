# --- Plot accession-level guard cell length against tau ----------------------

source("r/header.R")

selected_model = read_rds("objects/selected_model.rds")

df_new = crossing(
  phy = unique(selected_model$data$phy),
  leaftype = unique(selected_model$data$leaftype),
  lightintensity = unique(selected_model$data$lightintensity),
  lighttreatment = unique(selected_model$data$lighttreatment),
  logtausd = 0
) |>
  mutate(accession = phy,
         variable = paste0("...", row_number())) |> full_join(
           selected_model$data |>
             mutate(logfgmax = loggi - loggmax) |>
             summarize(
               logfgmax = median(logfgmax),
               loggi = median(loggi),
               loggmax = median(loggmax),
               .by = c(leaftype, lightintensity, lighttreatment)
             ),
           by = join_by(leaftype, lightintensity, lighttreatment)
         )

df_pred_loggcl = posterior_epred(selected_model, newdata = df_new, re_formula = ~ phy, resp = "loggcl") |>
  as_draws_df() |>
  summarize_draws() |>
  full_join(df_new, by = "variable") |>
  select(accession, leaftype, lightintensity, lighttreatment, loggcl = median)

df_pred_logtau = posterior_epred(selected_model, newdata = df_new, re_formula = ~ phy, resp = "logtaumean") |>
  as_draws_df() |>
  summarize_draws() |>
  full_join(df_new, by = "variable") |>
  select(accession, leaftype, lightintensity, lighttreatment, logtau = median)

Sigma_median = selected_model |>
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

df_acc = selected_model$data |>
  summarize(
    loggcl = median(loggcl),
    logtaumean = median(logtaumean), 
    .by = c(phy, leaftype, lightintensity, lighttreatment)) 

df_mu = df_acc |>
  summarize(
    loggcl = mean(loggcl),
    logtaumean = mean(logtaumean), 
    .by = c(leaftype, lightintensity, lighttreatment))

df_ellipse = df_mu |>
  mutate(ell = map2(loggcl, logtaumean, \(.x, .y) {
    ell = ellipse_points(mu = c(.x, .y), Sigma_median, level = 0.95)
  })) |>
  unnest(cols = ell) |>
  select(-loggcl, -logtaumean) |>
  rename(loggcl = x, logtaumean = y)

p = ggplot(df_acc, aes(exp(loggcl), exp(logtaumean), color = leaftype)) +
  geom_polygon(
    data = df_ellipse,
    aes(exp(loggcl), exp(logtaumean), fill = leaftype),
    color = "black",
    alpha = 0.25,
    inherit.aes = FALSE
  ) +
  geom_point() +
  facet_grid(lightintensity ~ lighttreatment) +
  scale_x_log10() +
  scale_y_log10(breaks = c(100, 200, 400)) +
  scale_fill_manual(values = c(col_amphi, col_pseudohypo)) +
  scale_color_manual(values = c(col_amphi, col_pseudohypo)) +
  labs(
    x = "Guard cell length ($\\si{\\micro\\meter}$, log-scale)",
    y = "$\\tau$ (s, log-scale)",
    color = "Leaf type:",
    fill = "Leaf type:"
  ) +
  theme(legend.position = "bottom")

gp1 = annotate_figure(
  p,
  top = ggpubr::text_grob("       Growth light intensity"),
  right = ggpubr::text_grob("Measurement light intensity        ", rot = -90)
)

options(
  tikzLatexPackages = c(
    getOption("tikzLatexPackages"),
    "\\usepackage{siunitx}"
  )
)

tikz(
  "figures/gcl-tau.tex",
  standAlone = TRUE,
  width = 5,
  height = 5
)
print(gp1)
dev.off()

system("cd figures; pdflatex gcl-tau.tex; rm gcl-tau.aux gcl-tau.log")
