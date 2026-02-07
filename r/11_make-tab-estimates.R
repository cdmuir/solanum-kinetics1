# Make table of parameter estimates for each accession
source("r/header.R")

phy = read_rds("data/phylogeny.rds")
A = vcv(phy, corr = TRUE)

fit_amphi_high = read_rds("objects/fit_amphi_high.rds")
fit_amphi_low = read_rds("objects/fit_amphi_low.rds")
fit_pseudohypo_high = read_rds("objects/fit_pseudohypo_high.rds")
fit_pseudohypo_low = read_rds("objects/fit_pseudohypo_low.rds")

assert_true(check_convergence(fit_amphi_high, convergence_criteria))
assert_true(check_convergence(fit_amphi_low, convergence_criteria))
assert_true(check_convergence(fit_pseudohypo_high, convergence_criteria))
assert_true(check_convergence(fit_pseudohypo_low, convergence_criteria))

tab_estimates = crossing(
  leaf_type = c("amphi", "pseudohypo"),
  light_intensity = c("high", "low")
) |>
  mutate(tab = map2(leaf_type, light_intensity, \(.x, .y) {
    fit = get(glue("fit_{.x}_{.y}"))
    
    sample_sizes = fit$data |>
      summarize(n_plant = n(), .by = c(lighttreatment, phy))
    
    df_new = crossing(
      lighttreatment = unique(fit$data$lighttreatment),
      phy = unique(fit$data$phy),
      logtausd = 0
    ) |>
      mutate(variable = paste0("...", row_number()))
    
    df_pred_tau = posterior_epred(fit, newdata = df_new, resp = "logtaumean") |>
      as_draws_df() |>
      summarise_draws() |>
      mutate(
        tau_estimate = exp(median),
        tau_lowerCI = exp(q5),
        tau_upperCI = exp(q95),
      ) |>
      select(variable, starts_with("tau_"))
    
    df_pred_gcl = posterior_epred(fit, newdata = df_new, resp = "loggcl") |>
      as_draws_df() |>
      summarise_draws() |>
      mutate(
        gcl_estimate = exp(median),
        gcl_lowerCI = exp(q5),
        gcl_upperCI = exp(q95),
      ) |>
      select(variable, starts_with("gcl_"))
    
    df_pred_fgmax = posterior_epred(fit, newdata = df_new, resp = "logfgmax") |>
      as_draws_df() |>
      summarise_draws() |>
      mutate(
        fgmax_estimate = exp(median),
        fgmax_lowerCI = exp(q5),
        fgmax_upperCI = exp(q95),
      ) |>
      select(variable, starts_with("fgmax_"))
    
    df_pred = df_pred_tau |>
      full_join(df_pred_gcl, by = join_by(variable)) |>
      full_join(df_pred_fgmax, by = join_by(variable)) |>
      full_join(df_new, by = join_by(variable)) |>
      select(-variable, -logtausd)
    
    full_join(sample_sizes, df_pred, by = join_by(lighttreatment, phy)) |>
      select(accession = phy,
             light_treatment = lighttreatment,
             n_plant,
             where(is.numeric)) |>
      pivot_longer(
        cols = -c(accession, light_treatment, n_plant),
        names_to = c("trait", "statistic"),
        names_sep = "_",
        values_to = "value"
      ) |>
      pivot_wider(names_from = statistic, values_from = value)
    
  })) |>
  unnest(cols = tab) |>
  select(
    accession,
    light_treatment,
    leaf_type,
    light_intensity,
    n_plant,
    trait,
    estimate,
    lowerCI,
    upperCI
  ) |>
  arrange(trait, accession, light_treatment, leaf_type, light_intensity)

# Write tab-estimates.csv
write_csv(tab_estimates, "tables/tab-estimates.csv")

# Write tab-estimates-dictionary.csv

# ---- build dictionary ----
dict = tibble(
  variable = names(tab_estimates),
  data_type = map_chr(tab_estimates, type_label),
  acceptable_values = map_chr(tab_estimates, acceptable_values),
) |>
  full_join(tibble(
    variable = c(
      "accession",
      "light_treatment",
      "leaf_type",
      "light_intensity",
      "n_plant",
      "trait",
      "estimate",
      "lowerCI",
      "upperCI"
    ),
    description = c(
      "TGRC accession",
      "Growth light intensity treatment",
      "Leaf type treatment",
      "Measurement light intensity treatment",
      "Number of replicate plants",
      "Trait name",
      "Posterior estimate",
      "Lower 95% CI",
      "Upper 95% CI"
    ),
  ), by = join_by(variable)) |>
  mutate(
    acceptable_values = na_if(acceptable_values, ""),
    acceptable_values = replace_na(acceptable_values, "")
  )

# ---- write to CSV ----
write_csv(dict, "tables/tab-estimates-dictionary.csv")
