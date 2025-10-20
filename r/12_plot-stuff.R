source("r/header.R")

fit = read_rds("objects/fit.rds")

# Predictions for log_gcl against log_tau
df_new1 = fit$data |>
  summarize(
    min_log_gcl = min(log_gcl),
    max_log_gcl = max(log_gcl),
    log_invpi = mean(log_invpi),
    .by = c(accession, light_intensity, light_treatment)
  ) |>
  crossing(i = seq(0, 1, length.out = 20)) |>
  mutate(log_gcl = min_log_gcl + i * (max_log_gcl - min_log_gcl),
         log_tau_sd = 0) |>
  select(-i, -min_log_gcl, -max_log_gcl)

df_pred1 = posterior_epred(fit, newdata = df_new1, re_formula = ~ accesssion) |>
  as_draws_df() |>
  summarise_draws() |>
  bind_cols(df_new1) |>
  rename(log_tau_mean = mean)

gp = map(unique(fit$data$accession), \(.acc) {
  
  df0 = fit$data |>
    filter(accession == .acc)
  df1 = df_pred1 |>
    filter(accession == .acc)
  
  ggplot(df1, aes(log_gcl, log_tau_mean, color = light_treatment, fill = light_treatment)) +
    geom_ribbon(mapping = aes(ymin = q5, ymax = q95), alpha = 0.5) +
    geom_line(linewidth = 2) +
    geom_point(data = df0) + 
    facet_wrap(~light_intensity)
})

pdf("figures/accession-tau.pdf", width = 8, height = 6)
for (i in seq_along(gp)) {
  print(gp[[i]] + ggtitle(names(gp)[i]))
}
dev.off()

# testing
.acc = "LA1777"
df0 = fit$data |>
  filter(accession == .acc)
df1 = df_pred1 |>
  filter(accession == .acc)

ggplot(df1, aes(log_gcl, log_tau_mean, color = light_treatment, fill = light_treatment)) +
  geom_ribbon(mapping = aes(ymin = q5, ymax = q95), alpha = 0.5) +
  geom_line(linewidth = 2) +
  geom_point(data = df0) + 
  facet_wrap(~light_intensity)
})


ggplot(.x, aes(log_invpi, log_tau_mean, color = light_treatment)) +
  geom_point() +
  facet_wrap(~light_intensity)


# Among species correlation between log_gcl and log_tau
df_new1 = fit$data |>
  summarize(
    log_gcl = mean(log_gcl), 
    log_invpi = mean(log_invpi),
    .by = c(accession, light_treatment)
  ) |>
  crossing(
    log_tau_sd = 0,
    light_intensity = unique(fit$data$light_intensity)
  )
  

df_pred1 = posterior_epred(fit, newdata = df_new1, re_formula = NA) |>
  as_draws_df() |>
  summarise_draws() |>
  bind_cols(df_new1)

# How are gcl and invpi related?
ggplot(filter(df_pred1, light_intensity == "high"), aes(log_gcl, log_invpi, color = light_treatment)) +
  geom_point()

# "raw" accession-level correlation between gcl and tau
ggplot(df_pred1, aes(log_gcl, mean, ymin = q5, ymax = q95, color = light_treatment)) +
  facet_wrap(~light_intensity) +
  geom_pointrange()

# "raw" accession-level correlation between invpi and tau
ggplot(df_pred1, aes(log_invpi, mean, ymin = q5, ymax = q95, color = light_treatment)) +
  facet_wrap(~light_intensity) +
  geom_pointrange()
