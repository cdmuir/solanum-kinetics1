# trying out alternative way to fit inertia model by setting k constant within an accession, but allowing other parameters to vary among curves.
source("r/header.R")

plant_info = read_rds("data/plant_info.rds") |>
  select(acc_id, light_treatment)

joined_data = read_rds("data/joined-data.rds") |>
  filter(acc == "LA2172", str_detect(curve, "amphi")) |>
  unite("acc_id", acc, id, sep = "-") |>
  mutate(light = str_extract(curve, "150|2000")) |>
  left_join(plant_info, by = "acc_id")

joined_data

form = gsw ~ exp(loggmax) * (1 - (ustar ^ ik + (uinit ^ ik - ustar ^ ik) * exp(-t_sec / exp(logtau))) ^ (1 / ik))

bform = bf(form, loggmax ~ -1 + curve, ustar ~ -1 + curve, uinit ~ -1 + curve, ik ~ light_treatment, logtau ~ curve, nl = TRUE)

pri = c(
  prior(
    uniform(0, 1),
    nlpar = "ustar",
    class = "b",
    lb = 0,
    ub = 1
  ),
  prior(
    uniform(0, 1),
    nlpar = "uinit",
    class = "b",
    lb = 0,
    ub = 1
  ),
  prior(
    normal(log(300), 1),
    nlpar = "logtau",
    class = "b",
    lb = log(10),
    ub = log(10000)
  )
)

# get_prior(bform, data = joined_data)
thin = 1
fit =   brm(
  formula = bform,
  data = joined_data,
  prior = pri,
  iter = thin * 2000,
  thin = thin,
  chains = 1,
  cores = 1,
  backend = "cmdstanr",
  control = list(adapt_delta = 0.8),
  seed = 20251230
)

bayes_R2(fit)

df1 = fit$data |>
  mutate(y_hat = predict(fit)[,"Estimate"])
  

ggplot(df1, aes(t_sec, gsw)) +
  facet_wrap(~curve) +
  geom_point() +
  geom_line(aes(y = y_hat), color = "tomato")

fit |>
  as_draws_df()
