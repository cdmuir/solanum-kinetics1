# figure showing negative relationship between estimates of phylogenetic and individual level relationships between tau and gcl
fits = read_rds("objects/fits.rds")
fits$fit[[5]] |>
  as_draws_df() |>
  select(cor_phy__logtaumean_Intercept__loggcl_Intercept, b_logtaumean_loggcl) |>
  ggplot(aes(x = cor_phy__logtaumean_Intercept__loggcl_Intercept, y = b_logtaumean_loggcl)) +
  geom_point() 
  