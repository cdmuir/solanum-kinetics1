# Plot accession-level gcl against tau (a bit messy, but getting close)
source("r/header.R")

fit_amphi = read_rds("objects/best_amphi_model.rds")

df_new = crossing(
  phy = unique(fit_amphi$data$phy),
  lightintensity = unique(fit_amphi$data$lightintensity),
  lighttreatment = unique(fit_amphi$data$lighttreatment),
  logtausd = 0
) |>
  mutate(accession = phy,
         variable = paste0("...", row_number())) |> full_join(
           fit_amphi$data |>
             summarize(
               logfgmax = median(logfgmax),
               .by = c(lightintensity, lighttreatment)
             ),
           by = join_by(lightintensity, lighttreatment)
         )

df_pred_loggcl = posterior_epred(fit_amphi, newdata = df_new, re_formula = ~ phy, resp = "loggcl") |>
  as_draws_df() |>
  summarize_draws() |>
  full_join(df_new, by = "variable") |>
  select(accession, lightintensity, lighttreatment, loggcl = median)

df_pred_logtau = posterior_epred(fit_amphi, newdata = df_new, resp = "logtaumean") |>
  as_draws_df() |>
  summarize_draws() |>
  full_join(df_new, by = "variable") |>
  select(accession, lightintensity, lighttreatment, logtau = median)

Sigma_median = fit_amphi |>
  as_draws_df() |>
  select(starts_with("."), sd_phy__loggcl_Intercept, sd_phy__logtaumean_Intercept, cor_phy__logtaumean_Intercept__loggcl_Intercept) |>
transmute(vx = sd_phy__loggcl_Intercept^2,
            vy = sd_phy__logtaumean_Intercept^2,
            cxy = cor_phy__logtaumean_Intercept__loggcl_Intercept * sd_phy__loggcl_Intercept * sd_phy__logtaumean_Intercept) %>%
  summarise(across(everything(), median)) %>%
  { matrix(c(.$vx, .$cxy, .$cxy, .$vy),
           nrow = 2, byrow = TRUE,
           dimnames = list(c("X","Y"), c("X","Y"))) }

df_acc = fit_amphi$data |>
  summarize(
    loggcl = median(loggcl),
    logtaumean = median(logtaumean), 
    .by = c(phy, lightintensity, lighttreatment)) 

df_mu = df_acc |>
  summarize(
    loggcl = mean(loggcl),
    logtaumean = mean(logtaumean), 
    .by = c(lightintensity, lighttreatment))

df_ellipse = df_mu |>
  mutate(ell = map2(loggcl, logtaumean, \(.x, .y) {
    ell = ellipse_points(mu = c(.x, .y), Sigma_median, level = 0.95)
  })) |>
  unnest(cols = ell) |>
  select(-loggcl, -logtaumean) |>
  rename(loggcl = x, logtaumean = y)


ggplot(df_acc, aes(loggcl, logtaumean)) +
  geom_path(data = df_ellipse, aes(loggcl, logtaumean), inherit.aes = FALSE) +
  geom_point() +
  facet_grid(lightintensity ~ lighttreatment)


# stuff from chatgpt on plotting ellipse below
ellipse_points <- function(mu, Sigma, level = 0.95, n = 200) {
  stopifnot(length(mu) == 2, all(dim(Sigma) == c(2, 2)))
  r <- sqrt(qchisq(level, df = 2))           # radius for chosen level
  theta <- seq(0, 2*pi, length.out = n)
  
  # unit circle
  circle <- rbind(cos(theta), sin(theta))
  
  # transform circle -> ellipse: mu + A %*% circle, where A A^T = Sigma
  A <- chol(Sigma)                           # upper-triangular
  pts <- t(circle) %*% A                     # (n x 2)
  tibble(x = mu[1] + pts[,1], y = mu[2] + pts[,2])
}

mu <- c(0, 0)
Sigma <- matrix(c(4, 1.5,
                  1.5, 2), nrow = 2, byrow = TRUE)

ell <- ellipse_points(mu, Sigma, level = 0.95)

ggplot(df, aes(x, y)) +
  geom_point() +
  geom_path(data = ell, aes(x, y), inherit.aes = FALSE)
