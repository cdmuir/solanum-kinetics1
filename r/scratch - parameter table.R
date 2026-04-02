# Make table of all model parameter estimates and CIs
source("r/header.R")

# working here
fit = read_rds("objects/best_model.rds")
ci_level = 0.95
digits = 2 

summ = summarise_draws(fit, estimate = median, quantile2, .args = list(probs = c((1 - ci_level) / 2, 1 - (1 - ci_level) / 2)))

# OLD FUNCTION

summ_df = as.data.frame(summ)
summ_df$param = rownames(summ_df)
colnames(summ_df) = c("Estimate", "Est.Error", "CI_low", "CI_high", "Parameter")

# Classify parameters
summ_df = summ_df |>
  mutate(
    Type = case_when(
      grepl("^b_", Parameter) &
        !grepl("sigma|Intercept", Parameter) ~ "fixed effects",
      grepl("^b_(aa|llma)_Intercept$", Parameter) ~ "fixed effects",
      grepl("^bsp_", Parameter) ~ "fixed effects",
      grepl("^b_", Parameter) &
        grepl("sigma", Parameter) ~ "distributional parameters",
      grepl("^sd_", Parameter) ~ "random effect SDs",
      TRUE ~ "Other",
    ),
    resp = case_when(
      str_detect(Parameter, "_aa_") ~ "aa",
      str_detect(Parameter, "_llma_") ~ "lma",
      TRUE ~ NA_character_
    )
  )

# Clean parameter names
summ_df <- summ_df %>%
  mutate(
    Parameter = case_when(
      # AA
      Parameter == "b_aa_Intercept" ~ "$\\mathrm{AA}$ intercept (shade, low light)",
      Parameter == "b_aa_light_treatmenthigh" ~ "effect of sun on $\\mathrm{AA}$",
      Parameter == "b_aa_light_intensity2000" ~ "effect of high light on $\\mathrm{AA}$",
      Parameter == "b_aa_light_treatmenthigh:light_intensity2000" ~ "effect of sun $\\times$ high light interaction on $\\mathrm{AA}$",
      Parameter == "bsp_aa_millma" ~ "effect of $\\log \\mathrm{LMA}$ on $\\mathrm{AA}$",
      Parameter == "b_sigma_aa_Intercept" ~ "$\\log \\sigma$ intercept (shade, low light)",
      Parameter == "b_sigma_aa_light_treatmenthigh" ~ "effect of sun on $\\log \\sigma$",
      Parameter == "b_sigma_aa_light_intensity2000" ~ "effect of high light on $\\log \\sigma$",
      Parameter == "sd_acc__aa_Intercept" ~ "$\\mathrm{AA}$ among populations",
      Parameter == "sd_acc__aa_light_treatmenthigh" ~ "effect of sun on $\\mathrm{AA}$ among populations",
      Parameter == "sd_acc__aa_light_intensity2000" ~ "effect of high light on $\\mathrm{AA}$ among populations",
      Parameter == "sd_acc__aa_light_treatmenthigh:light_intensity2000" ~ "effect of sun $\\times$ high light interaction on $\\mathrm{AA}$ among populations",
      
      # LMA
      Parameter == "b_llma_Intercept" ~ "$\\log \\mathrm{LMA}$ intercept (shade)",
      Parameter == "b_llma_light_treatmenthigh" ~ "effect of sun on $\\log \\mathrm{LMA}$",
      Parameter == "b_sigma_llma_Intercept" ~ "$\\log \\sigma$ intercept (shade)",
      Parameter == "b_sigma_llma_light_treatmenthigh" ~ "effect of sun on $\\log \\sigma$",
      Parameter == "sd_acc__llma_Intercept" ~ "$\\log \\mathrm{LMA}$ among populations",
      Parameter == "sd_acc__llma_light_treatmenthigh" ~ "effect of sun on $\\log \\mathrm{LMA}$ among populations"
    )
  )

# Select and format
summ_df |>
  filter(Type != "Other") |>
  mutate(across(c(Estimate, CI_low, CI_high), ~ round(., digits))) |>
  mutate(`Estimate [95% CI]` = paste0(Estimate, " [", CI_low, ", ", CI_high, "]")) |>
  select(resp, Type, Parameter, `Estimate [95% CI]`) |>
  arrange(resp, match(
    Type,
    c(
      "fixed effects",
      "distributional parameters",
      "random effect SDs"
    )
  )) |>
  as_tibble()
}

summary_tbl = make_brms_summary_table(fit_aa2)

tbl_aa_fixed = summary_tbl |> filter(Type == "fixed effects", resp == "aa") |> select(-Type)
tbl_aa_dist  = summary_tbl |> filter(Type == "distributional parameters", resp == "aa") |> select(-Type)
tbl_aa_rand  = summary_tbl |> filter(Type == "random effect SDs", resp == "aa") |> select(-Type)
tbl_lma_fixed = summary_tbl |> filter(Type == "fixed effects", resp == "lma") |> select(-Type)
tbl_lma_rand  = summary_tbl |> filter(Type == "random effect SDs", resp == "lma") |> select(-Type)

bind_rows(tbl_aa_fixed,
          tbl_aa_dist,
          tbl_aa_rand,
          tbl_lma_fixed,
          tbl_lma_rand) |>
  select(-resp) |>
  kable(
    col.names = c("Parameter", "Estimate [95\\% CI]"),
    booktabs = TRUE,
    format = "latex",
    escape = FALSE
  ) |>
  kable_styling(latex_options = c("striped")) |>
  pack_rows(
    "Response: $\\mathrm{AA}$",
    1,
    nrow(tbl_aa_fixed) + nrow(tbl_aa_dist) + nrow(tbl_aa_rand),
    escape = FALSE
  ) |>
  pack_rows(
    "~~Fixed effects",
    1,
    nrow(tbl_aa_fixed),
    escape = FALSE,
    bold = FALSE,
    italic = TRUE
  ) |>
  pack_rows(
    "~~Distributional parameters",
    nrow(tbl_aa_fixed) + 1,
    nrow(tbl_aa_fixed) + nrow(tbl_aa_dist),
    escape = FALSE,
    bold = FALSE,
    italic = TRUE
  ) |>
  pack_rows(
    "~~Random effect SDs",
    nrow(tbl_aa_fixed) + nrow(tbl_aa_dist) + 1,
    nrow(tbl_aa_fixed) + nrow(tbl_aa_dist) + nrow(tbl_aa_rand),
    escape = FALSE,
    bold = FALSE,
    italic = TRUE
  ) |>
  pack_rows(
    "Response: $\\log~\\\\mathrm{LMA}$",
    nrow(tbl_aa_fixed) + nrow(tbl_aa_dist) + nrow(tbl_aa_rand) + 1,
    nrow(tbl_aa_fixed) + nrow(tbl_aa_dist) + nrow(tbl_aa_rand) + nrow(tbl_lma_fixed) + nrow(tbl_lma_rand),
    escape = FALSE
  ) |>
  pack_rows(
    "~~Fixed effects",
    nrow(tbl_aa_fixed) + nrow(tbl_aa_dist) + nrow(tbl_aa_rand) + 1,
    nrow(tbl_aa_fixed) + nrow(tbl_aa_dist) + nrow(tbl_aa_rand) + nrow(tbl_lma_fixed),
    escape = FALSE,
    bold = FALSE,
    italic = TRUE
  ) |>
  pack_rows(
    "~~Random effect SDs",
    nrow(tbl_aa_fixed) + nrow(tbl_aa_dist) + nrow(tbl_aa_rand) + nrow(tbl_lma_fixed) + 1,
    nrow(tbl_aa_fixed) + nrow(tbl_aa_dist) + nrow(tbl_aa_rand) + nrow(tbl_lma_fixed) + nrow(tbl_lma_rand),
    escape = FALSE,
    bold = FALSE,
    italic = TRUE
  ) |>
  column_spec(1, width = "4.5in") |>
  column_spec(2, width = "2in")

```
