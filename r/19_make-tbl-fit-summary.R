# Make table of all model parameter estimates and CIs
source("r/header.R")

loggcl_latex = "$\\log \\left( l_\\text{gc} \\right)$"
logitfgmax_latex = "$\\text{logit} \\left( f_\\text{gmax} \\right)$"
loglambda_latex = "$\\log \\left( \\lambda \\right)$"
logtau_latex = "$\\log \\left( \\tau \\right)$"

# working here
fit = read_rds("objects/best_model.rds")
ci_level = 0.95
digits = 2 

df_summary1 = summarise_draws(fit,
                              estimate = median,
                              quantile2,
                              .args = list(probs = c((1 - ci_level) / 2, 1 - (1 - ci_level) / 2), names = FALSE)) |>
  mutate(
    across(where(is_double), ~ formatC(., format = "f", digits = digits)),
    `Estimate [95% CI]` = glue("{estimate} [{quantile2.1}, {quantile2.2}]")
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
    resp = case_when(
      str_detect(variable, "_loggcl_") ~ loggcl_latex,
      str_detect(variable, "_logitfgmax_") ~ logitfgmax_latex,
      str_detect(variable, "_loglambdamean_") ~ loglambda_latex,
      str_detect(variable, "_logtaumean_") ~ logtau_latex,
      TRUE ~ NA_character_
    )
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
      str_detect(variable, "_logitfgmax$") ~ "$\\text{logit} \\left( f_\\text{gmax} \\right)$",
      str_detect(variable, "_loggcl$") ~ "$\\log \\left( l_\\text{gc} \\right)$",
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
    group2 = case_when(
      group1 == "accession" ~ "population (nonphylogenetic)",
      group1 == "phy" ~ "phylogenetic",
      group1 == "accid" ~ "among-individual"
    ),
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
    group2 = case_when(
      group1 == "cor_accession" ~ "population (nonphylogenetic)",
      group1 == "cor_phy" ~ "phylogenetic",
      group1 == "cor_accid" ~ "among-individual"
    ),
    resp1 = case_when(
      resp1 == "loggcl" ~ loggcl_latex,
      resp1 == "logitfgmax" ~ logitfgmax_latex,
      resp1 == "loglambdamean" ~ loglambda_latex,
      resp1 == "logtaumean" ~ logtau_latex,
      TRUE ~ NA_character_
    ),
    resp2 = case_when(
      resp2 == "loggcl" ~ loggcl_latex,
      resp2 == "logitfgmax" ~ logitfgmax_latex,
      resp2 == "loglambdamean" ~ loglambda_latex,
      resp2 == "logtaumean" ~ logtau_latex,
      TRUE ~ NA_character_
    ),
    description = glue("correlation between in {resp1} and {resp2}")
  )

write_rds(list(
  df_fixed = df_fixed,
  df_random_sd = df_random_sd,
  df_random_cor = df_random_cor
), "objects/tbl-fit-summary.rds")

library(knitr)
library(kableExtra)
bind_rows(
  df_random_sd |>
    select(group2, Type, resp, `Estimate [95% CI]`, description),
  
  df_random_cor |>
    select(group2, Type, resp1, resp2, `Estimate [95% CI]`, description)
) |>
  arrange(group2, Type) |>
  select(description, `Estimate [95% CI]`) |>
  kable(
    col.names = c("Parameter description", "Estimate [95\\% CI]"),
    booktabs = TRUE,
    format = "latex",
    escape = FALSE
  ) |>
  kable_styling(latex_options = c("striped")) |>
  pack_rows(
    "among-individual random effects",
    1,
    3,
    escape = FALSE,
    bold = TRUE,
    italic = FALSE
  ) |>
  pack_rows(
    "phylogenetic random effects",
    4,
    13,
    escape = FALSE,
    bold = TRUE,
    italic = FALSE
  ) |>
  pack_rows(
    "population (nonphylogenetic) random effects",
    14,
    23,
    escape = FALSE,
    bold = TRUE,
    italic = FALSE
  ) |>
  column_spec(1, width = "4.5in") |>
  column_spec(2, width = "2in")
