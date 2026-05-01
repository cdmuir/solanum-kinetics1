# Plot variance components for gcl
source("r/header.R")

fit_gcl = read_rds("objects/fit-gcl.rds")

assert_true(check_convergence(fit_gcl, convergence_criteria))

variance_parts = fit_gcl |>
  as_draws_df() |>
  select(starts_with("."), starts_with("sd_"), sigma) |>
  mutate(across(c(starts_with("sd_"), sigma), ~ .x^2)) |>
  rename_with(~ str_replace(.x, "sd_", "var_"), starts_with("sd_")) |>
  rename(
    var_accid = var_acc_id__Intercept,
    var_accession = var_accession__Intercept,
    var_resid = sigma
  ) |>
  mutate(
    var_phy = var_phy__Intercept + var_phy__light_treatmenthigh + var_phy__surfaceupper + `var_phy__light_treatmenthigh:surfaceupper`,
    var_total  = var_accid + var_accession + var_phy + var_resid,
    prop_accid    = var_accid    / var_total,
    prop_accession = var_accession / var_total,
    prop_phy       = var_phy       / var_total,
    prop_resid     = var_resid     / var_total
  )

summary_parts = variance_parts |>
  select(starts_with("."), starts_with("prop_")) |>
  pivot_longer(cols = starts_with("prop_"), names_to = "part", values_to = "prop") |>
  group_by(part) |>
  point_interval() |>
    mutate(part1 = case_when(
      part == "prop_resid" ~ "Within leaf + error",
      part == "prop_accid" ~ "Among individuals",
      part == "prop_accession" ~ "Among accessions\n(nonphylogenetic)",
      part == "prop_phy" ~ "Among accessions\n(phylogenetic)"
    ))

ggplot(summary_parts, aes(0, prop, fill = part1)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_viridis_d() +
  labs(y = "proportion of variance", fill = "variance component") +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.x = element_blank())
