# VERY MUCH IN FLUX - take away so far is that fit is much better on inverse of k because of numerical issue
# I THINK THE THING TO TEST IS WHETHER EFFECT OF GOP / GMAX ON TAU GOES AWAY IN THIS MODEL COMPARED TO THE CDWEIBULL ONE
# Fit model 3 (power function) to each curve separately
# this is not quite correct since it's using raw g, not g normalized distance to gmax
source("r/header.R")

sk_dir = c("objects/sk-curves3/")
if (!dir.exists(sk_dir)) {
  dir.create(sk_dir)
}

joined_data = read_rds("data/joined-data.rds")

plan(multisession, workers = 9)
# gmin = minimum stomatal conductance
# gmax = anatomical maximum conductance
# ginit = initial conductance
# gstar = target conductance

f1 = function(t_sec, gmax, gmin, ginit, gstar, k, tau) {
  gmax + (gmin - gmax) * (((gstar - gmax) / (gmin - gmax)) ^ (1 / k) + (((ginit - gmax) / (gmin - gmax)) ^ (1 / k) - ((gstar - gmax) / (gmin - gmax)) ^ (1 / k)) * exp(-t_sec / tau)) ^ k}

# sample(unique(joined_data$curve), 1)
df1 = filter(joined_data, curve == "LA2951-M_amphi_2000") # good example of stretched curve
df1 = filter(joined_data, curve == "LA1316-T_amphi_150")
df2 = df1 |>
  select(-gsw) |>
  mutate(
    gsw = f1(
      t_sec = t_sec,
      gmax = gmax,
      gmin = 0.02,
      ginit = 0.1,
      gstar = 0.05,
      k = 1,
      tau = 100
    )
  )

ggplot(df1, aes(t_sec, gsw)) +
  geom_point() +
  geom_line(data = df2)

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



bform3 = bf(form3, gmin ~ 1, gstar ~ 1, ginit ~ 1, ik ~ 1, logtau ~ 1, nl = TRUE)
# bform3 = bf(form3, gmin ~ 1, gstar ~ 1, ginit ~ 1, logtau ~ 1, nl = TRUE)

pri = c(
  prior(
    normal(0.05, 1),
    nlpar = "gmin",
    lb = 0,
    ub = 0.1
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
    ub = first(df1$gmax) #max(df1$gsw) * 1.1
  ),
  prior(
    normal(0, 100),
    nlpar = "ik"
  ),
  prior(
    normal(log(300), 1),
    nlpar = "logtau",
    lb = log(10),
    ub = log(10000)
  )
)

fit = fit_rh1(
  formula = bform3,
  data = df1,
  # data = mutate(df1, k = 1),
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


rh_curves |>
  split(~ curve) |>
  magrittr::extract(1:9) |>
  future_iwalk(\(df, curve_id) {
    file = paste0(sk_dir, curve_id, ".rds")
    
    fit = fit_rh1(
      formula = bform3,
      data = df,
      prior = pri,
      thin = 2,
      adapt_delta = 0.8,
      seed = 360036340 + df$ci[1]
    )
    
    write_rds(fit, file)
    
  }, .progress = TRUE, .options = furrr_options(seed = TRUE))

