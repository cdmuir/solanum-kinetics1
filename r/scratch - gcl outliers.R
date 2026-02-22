# quick look at gcl outliers on 2/21. There are 3 leaves I should go back and look at: LA0407-P (this one def does not look like other LA0407), LA4116-I, LA0429-F
source("r/header.R")

# fit when these individuals were in htere
fit = read_rds("objects/best_model.rds")

hist(fit$data$loggcl)

fit1 = lm(loggcl ~ lighttreatment + curve_type + phy, data = fit$data)
anova(fit1)
qqnorm(resid(fit1))
fit$data$resid = resid(fit1)

fit$data |>
  select(phy, accid, resid) |>
  arrange(desc(abs(resid))) |> as_tibble() |>
  print(n = 50)

1 LA0407  LA0407-P   0.537
2 LA0407  LA0407-P   0.537
3 LA4116  LA4116-I   0.531
4 LA4116  LA4116-I   0.531
5 LA0429  LA0429-F   0.501
6 LA0429  LA0429-F   0.501
7 LA0407  LA0407-P   0.492
8 LA0407  LA0407-P   0.492
9 LA4116  LA4116-I   0.478
10 LA4116  LA4116-I   0.478
11 LA0429  LA0429-F   0.452
12 LA0429  LA0429-F   0.452

fit$data |>
  filter(phy == "LA0407", curve_type == "amphi", lightintensity == "high") |>
  select(accid, loggcl)
