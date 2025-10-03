source("r/header.R")

plant_info = read_rds("data/plant-info.rds")
stomata = read_rds("data/stomata.rds") |>
  mutate(acc_id = paste(accession, replicate, sep = "-"),
         log_gcl = log(stomatal_ratio * upper_guard_cell_length_um + (1 - stomatal_ratio) * lower_guard_cell_length_um))
rh_curves = read_rds("data/trimmed_rh_curves.rds") |>
  summarize(max_gsw = max(gsw, na.rm = TRUE),
            .by = c(acc_id, light_treatment, light_intensity, curve_type)) |>
  mutate(curve_type = case_when(
    curve_type == "1-sided RH" ~ "pseudohypo",
    curve_type == "2-sided RH" ~ "amphi"
  ))

# nls ----
nls_summary = read_rds("objects/nls_summary.rds") |>
  separate_wider_delim(curve,
                       "_",
                       names = c("acc_id", "curve_type", "light_intensity"),
                       cols_remove = FALSE) |>
  mutate(acc = str_extract(acc_id, "^(LA[0-9]{4}A*|nelsonii|sandwicense)"))

nls_summary |>
  filter(curve_type == "amphi", light_intensity == "2000", parameter == "tau", model == "cdweibull", !is.na(value)) |>
  summarize(tau = median(value), .by = acc) |>
  arrange(tau) |>
  print(n = Inf)
# 
## CDWeibull usually preferred
nls_summary |>
  summarize(
    aic = first(aic),
    .by = c(curve, model)
  ) |>
  pivot_wider(names_from = model, values_from = aic) |>
  mutate(delta_aic = cdweibull - vico) |>
  filter(!is.na(delta_aic)) |>
  arrange(desc(delta_aic)) 

df_nls = nls_summary |>
  pivot_wider(names_from = "parameter") |>
  filter(!is.na(aic)) |>
  left_join(dplyr::select(plant_info, accession, acc_id, light_treatment), by = join_by(acc_id)) |>
  left_join(stomata, by = join_by(acc_id, accession)) |>
  left_join(rh_curves, by = join_by(acc_id, light_treatment, light_intensity, curve_type)) |>
  mutate(log_tau = log(tau)) |>
  filter(model == "cdweibull", log_tau < 7)

ggplot(df_nls, aes(log_gcl, tau, color = light_treatment)) +
  facet_grid(light_intensity ~ curve_type, scales = "free") +
  geom_point(alpha = 0.5) +
  scale_x_log10() +
  scale_y_log10()

fit_nls = brm(
  bf(
    log_tau ~ log_gcl * light_treatment * light_intensity + (1 |
                                                               acc_id) + (1 + log_gcl |
                                                                            accession)
  ) +
    bf(log_gcl ~ light_treatment + (1 | accession)) +
    set_rescor(FALSE),
  data = filter(df_nls, curve_type == "amphi"),
  chains = 1,
  backend = "cmdstanr"
)

conditional_effects(fit_nls)


df_new = crossing(
  light_treatment = c("low", "high"),
  light_intensity = c("150", "2000"),
  curve_type = c("pseudohypo", "amphi"),
  accession = unique(df_nls$accession),
)

ranef(fit_nls)$accession[,,"Intercept"] |>
  as_tibble() |>
  arrange(Estimate) |>
  mutate(x = row_number()) |>
  ggplot(aes(x, Estimate, ymin = `Q2.5`, ymax = `Q97.5`)) +
  geom_pointinterval()

filter(df_nls, accession == "LA1777") |>
  ggplot(aes(log_gcl, log_tau, color = light_intensity)) +
  facet_wrap(light_treatment ~ curve_type, scales = "free") +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  scale_x_log10() +
  scale_y_log10()

df_nls |>
  summarize(
    log_tau = mean(log_tau, na.rm = TRUE),
    log_gcl = mean(log_gcl, na.rm = TRUE),
    .by = c(accession, light_treatment, light_intensity, curve_type)
  ) |>
  filter(curve_type == "amphi") |>
  ggplot(aes(log_gcl, log_tau, color = light_treatment)) +
  facet_wrap(light_treatment ~ light_intensity) +
  geom_point()

df_nls |>
  filter(curve_type == "amphi") |>
  ggplot(aes((gf + dg) / ((lower_gmax + upper_gmax) / 1000), log_tau)) +
  geom_point() 

df_nls |>
  filter(curve_type == "amphi") |>
  ggplot(aes(max_gsw, (gf + dg))) +
  geom_point()

df_nls |>
  filter(curve_type == "amphi", accession == "LA2933") |>
  ggplot(aes((gf + dg) / total_stomatal_density_mm2, log_tau, color = light_intensity)) +
  geom_point() 

df_nls |>
  filter(curve_type == "amphi", accession == "LA2172") |>
  ggplot(aes((gf + dg) / (lower_gmax + upper_gmax), log_tau, color = light_intensity)) +
  geom_point() +
  scale_x_log10() 

# brms ---
sk_curves = list.files("objects/sk-curves")

sk_pars = sk_curves |>
  map_dfr(\(.f) {
    .m = read_rds(paste0("objects/sk-curves/", .f))
    
    fixef(.m) |>
      as_tibble(rownames = "parameter") |>
      pivot_longer(cols = -parameter,
                   names_to = "stat",
                   values_to = "value") |>
      mutate(curve = str_remove(.f, "\\.rds$")) |>
      separate_wider_delim(curve,
                           "_",
                           names = c("acc_id", "curve_type", "light_intensity"))
  }, .progress = TRUE)
  
sk_pars = sk_pars |>
  pivot_wider(names_from = "stat") |>
  mutate(parameter = str_remove(parameter, "_Intercept"))

sk_pars |> 
  filter(
    parameter != "dg" | Estimate < 1,
    parameter != "dg" | Estimate > -1,
    parameter != "gf" | Estimate > 0,
    parameter != "lambda" | Estimate < 1500,
    parameter != "tau" | Estimate < 1000,
  ) |>
  ggplot(aes(Estimate, color = curve_type)) +
  facet_wrap(light_intensity ~ parameter, scales = "free") +
  geom_density()

sk_pars |>
  filter(parameter == "tau") |>
  mutate(Estimate = log(Estimate)) |>
  arrange(desc(Estimate)) |>
  print(n = 100)
  pull(Estimate) |> hist()

  sk_pars |>
    select(parameter:Estimate) |>
    pivot_wider(names_from = "parameter",
                values_from = "Estimate") |>
    ggplot(aes(dg, lambda)) +
    geom_point() +
    xlim(-1, 1) +
    ylim(-3, 3)
  