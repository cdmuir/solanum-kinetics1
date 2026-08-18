# --- Make the table of all model parameter estimates and confidence intervals ---

source("r/header.R")

selected_model = read_rds("objects/selected_model.rds")
ci_level = 0.95
digits = 2 

df_summary1 = summarize_draws(selected_model,
                              estimate = median,
                              quantile2,
                              .args = list(probs = c((1 - ci_level) / 2, 1 - (1 - ci_level) / 2), names = FALSE)) |>
  mutate(
    across(where(is_double), ~ formatC(., format = "f", digits = digits)),
    `Estimate [95\\% CI]` = glue("{estimate} [{quantile2.1}, {quantile2.2}]")
    ) |>
  filter(
    !str_detect(variable, "^Intercept_"),
    !str_detect(variable, "^nu"),
    !str_detect(variable, "^r_"),
    !str_detect(variable, "^rescor__"),
    !str_detect(variable, "^sigma_"),
    variable != "lprior",
    variable != "lp__") |>
# Classify parameters
  mutate(
    Type = case_when(
      str_detect(variable, "^b_") ~ "fixed effects",
      str_detect(variable, "^sd_") ~ "random effect SDs",
      str_detect(variable, "^cor_") ~ "random effect correlations",
      TRUE ~ NA_character_,
    ),
    resp = trait_latex_from_paramname(variable)
  )

# Fixed effects
df_fixed = df_summary1 |>
  filter(Type == "fixed effects") |>
  mutate(
    explanatory = case_when(
      str_detect(variable, "_Intercept$") ~ "intercept (shade, low light, amphi)",
      str_detect(variable, "_lighttreatmentsun$") ~ "sun",
      str_detect(variable, "_lightintensityhigh$") ~ "high light",
      str_detect(variable, "_leaftypepseudohypo$") ~ "pseudohypo leaf type",
      str_detect(variable, "_loggcl$") ~ trait_latex_label("loggcl"),
      str_detect(variable, "_loggi$") ~ trait_latex_label("loggi"),
      str_detect(variable, "_loggmax$") ~ trait_latex_label("loggmax"),
      TRUE ~ NA_character_
    ),
    description = map2_chr(resp, explanatory, \(.r, .e) {
      if (str_detect(.e, "^intercept")) {
      glue("{.r} {.e}")
      } else {
        glue("effect of {.e} on {.r}")
      }
    })
  )

# Random effects
df_random_sd = df_summary1 |>
  filter(Type == "random effect SDs") |>
  mutate(
    group1 = str_extract(variable, "(?<=_)[^_]+(?=__)"),
    group2 = recode_variance_component(group1),
    description = glue("SD in {resp}")
  )

df_random_cor = df_summary1 |>
  filter(Type == "random effect correlations") |>
  mutate(
    map_dfr(variable, \(.x) {
      .x |>
        str_remove_all("_Intercept") |>
      str_split("__") |>
      extract2(1) |>
      as_tibble() |>
      mutate(name = c("group1", "resp1", "resp2")) |>
      pivot_wider()
    }),
    group2 = recode_variance_component(str_remove(group1, "^cor_")),
    resp1 = trait_latex_label(resp1),
    resp2 = trait_latex_label(resp2),
    description = glue("correlation between in {resp1} and {resp2}")
  )

write_rds(list(
  df_fixed = df_fixed,
  df_random_sd = df_random_sd,
  df_random_cor = df_random_cor
), "objects/tbl-fit-summary.rds")
