# some notes, etc. relating to inertia model that might be useful for later

# gmin = minimum stomatal conductance
# gmax = anatomical maximum conductance
# ginit = initial conductance
# gstar = target conductance

f1 = function(t_sec, gmax, gmin, ginit, gstar, k, tau) {
  gmax + (gmin - gmax) * (((gstar - gmax) / (gmin - gmax)) ^ (1 / k) + (((ginit - gmax) / (gmin - gmax)) ^ (1 / k) - ((gstar - gmax) / (gmin - gmax)) ^ (1 / k)) * exp(-t_sec / tau)) ^ k}

# Formula with k, however fit is basically insensitive to k
form3 = gsw ~ gmax + (gmin - gmax) * (((gstar - gmax) / (gmin - gmax)) ^ (1 / k) + (((ginit - gmax) / (gmin - gmax)) ^ (1 / k) - ((gstar - gmax) / (gmin - gmax)) ^ (1 / k)) * exp(-t_sec / exp(logtau))) ^ k 

# try reparameterizing as inverse k - much better!
form3 = gsw ~ gmax + (gmin - gmax) * (((gstar - gmax) / (gmin - gmax)) ^ ik + (((ginit - gmax) / (gmin - gmax)) ^ ik - ((gstar - gmax) / (gmin - gmax)) ^ ik) * exp(-t_sec / exp(logtau))) ^ (1 / ik) 

# this derivation shows that when k = 1, model reduces to the vico model
# hence, when dynamics are insensitive to k, model will tend to estimate k near 1 and gop doesn't matter
gmax = 2.1; gmin = 0.01; ginit = 0.5; gstar = 0.1; k = 1; logtau = log(100); t_sec = 10
gmax + (gmin - gmax) * (((gstar - gmax) / (gmin - gmax)) ^ (1 / k) + (((ginit - gmax) / (gmin - gmax)) ^ (1 / k) - ((gstar - gmax) / (gmin - gmax)) ^ (1 / k)) * exp(-t_sec / exp(logtau))) ^ k 
gmax + (gmin - gmax) * (((gstar - gmax) / (gmin - gmax)) + (((ginit - gmax) / (gmin - gmax)) - ((gstar - gmax) / (gmin - gmax))) * exp(-t_sec / exp(logtau)))
gmax + (gmin - gmax) * ((gstar - gmax) / (gmin - gmax) + ((ginit - gmax) / (gmin - gmax) - (gstar - gmax) / (gmin - gmax)) * exp(-t_sec / exp(logtau)))
gmax + ((gstar - gmax) + ((ginit - gmax) - (gstar - gmax)) * exp(-t_sec / exp(logtau)))
gmax + (gstar - gmax + (ginit - gmax - gstar + gmax) * exp(-t_sec / exp(logtau)))
gmax + (gstar - gmax + (ginit - gstar) * exp(-t_sec / exp(logtau)))
gstar + (ginit - gstar) * exp(-t_sec / exp(logtau))

# this better shows sensitivity to k on inverse scale
tibble(
  k = 1 / 10 ^ seq(0, 3, 0.1),
  gsw = f1(10, 5, 0.02, 0.65, 0.27, k = k, tau = 228)
) |>
  ggplot(aes(1 /k, gsw)) +
  geom_line()

# example fits
# sample(unique(joined_data$curve), 1)
df1 = filter(joined_data, curve == "LA2951-M_amphi_2000") # good example of stretched curve
df1 = filter(joined_data, curve == "LA1316-T_amphi_150")

# try reparameterizing as inverse k - much better!
form3 = gsw ~ gmax + (gmin - gmax) * (((gstar - gmax) / (gmin - gmax)) ^ ik + (((ginit - gmax) / (gmin - gmax)) ^ ik - ((gstar - gmax) / (gmin - gmax)) ^ ik) * exp(-t_sec / exp(logtau))) ^ (1 / ik) 

bform3 = bf(form3, gmin ~ 1, gstar ~ 1, ginit ~ 1, ik ~ 1, logtau ~ 1, nl = TRUE)
pri = get_interia_prior(df1)

fit = fit_rh1(
  formula = bform3,
  data = df1,
  prior = pri,
  thin = 2,
  adapt_delta = 0.8,
  seed = 360036340
)

df2 = df1 |>
  select(-gsw) |>
  mutate(
    gsw = predict(fit)[,"Estimate"]
  )
ggplot(df1, aes(t_sec, gsw)) +
  geom_point() +
  geom_line(data = df2)


# g as a function of P
gP <- function(P, gmax, gmin, k) {
  gmax + (-gmax + gmin)*exp(k*log(P + 1))
}

tibble(
  P = seq(0, 6, 0.1),
  g = gP(P, gmax = 5, gmin = 0.01, k = -1.5)
) |>
  ggplot(aes(P, g)) +
  geom_line()


source("r/header.R")
form_inertia = gsw ~ gmax * (1 -
                               ((1 - gstar / gmax)^ik +
                                  ((1 - ginit / gmax)^ik -
                                     (1 - gstar / gmax)^ik) * exp(-t_sec / exp(logtau)))^(1 / ik))

bform_inertia = bf(form_inertia, gmax ~ 1, gstar ~ 1, ginit ~ 1, ik ~ 1, logtau ~ 1, nl = TRUE)

joined_data = read_rds("data/joined-data.rds")
df1 = filter(joined_data, curve == "LA2951-M_amphi_2000") # good example of stretched curve

pri = c(
  set_prior(
    "normal(5, 1)",
    nlpar = "gmax",
    lb = max(df1$gsw)
  ),
  set_prior(
    glue("normal({mu}, 1)", mu = min(df1$gsw)),
    nlpar = "gstar",
    lb = min(df1$gsw) * 0.95,
    ub = min(df1$gsw) * 1.05
  ),
  set_prior(
    glue("normal({mu}, 1)", mu = max(df1$gsw)),
    nlpar = "ginit",
    lb = max(df1$gsw) * 0.9,
    ub = first(df1$gmax)
  ),
  prior(
    normal(0, 100),
    nlpar = "ik",
    ub = 0
  ),
  prior(
    normal(log(300), 1),
    nlpar = "logtau",
    lb = log(10),
    ub = log(10000)
  )
)

fit = fit_rh1(
  formula = bform_inertia,
  data = df1,
  prior = pri,
  thin = 2,
  adapt_delta = 0.8,
  seed = 360036340
)
