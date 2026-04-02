# Make table of parameter estimates for each accession
source("r/header.R")

phy = read_rds("data/phylogeny.rds")
A = vcv(phy, corr = TRUE)

fit = read_rds("objects/best_model.rds")

# Sample sizes
sample_sizes = fit$data |>
  summarize(
    n_plant = n(),
    .by = c(phy, lighttreatment, lightintensity, leaftype)
  )

df_new = crossing(
  phy = unique(fit$data$phy),
  lighttreatment = unique(fit$data$lighttreatment),
  lightintensity = unique(fit$data$lightintensity),
  leaftype = unique(fit$data$leaftype),
  logtausd = 0,
  loglambdasd = 0
) |>
  mutate(accession = phy,
         variable = paste0("...", row_number()))

# Predictions for fgmax and gcl
df_pred_fgmax = get_posterior_epred(fit, df_new, "logitfgmax", "fgmax", inv = plogis)
df_pred_gcl = get_posterior_epred(fit, df_new, "loggcl", "gcl", inv = exp)

# Add predictions for fgmax and gcl to make predictions for tau and lambda
df_new1 = df_new |>
  left_join(df_pred_fgmax, by = join_by(variable)) |>
  left_join(df_pred_gcl, by = join_by(variable)) |>
  mutate(logitfgmax = qlogis(fgmax_estimate),
         loggcl = log(gcl_estimate))

df_pred_lambda = get_posterior_epred(fit, df_new1, "loglambdamean", "lambda", inv = exp)
df_pred_tau = get_posterior_epred(fit, df_new1, "logtaumean", "tau", inv = exp)

# Join and write
df_pred = df_pred_tau |>
  full_join(df_pred_lambda, by = join_by(variable)) |>
  full_join(df_pred_gcl, by = join_by(variable)) |>
  full_join(df_pred_fgmax, by = join_by(variable)) |>
  full_join(df_new, by = join_by(variable)) |>
  select(-variable, -logtausd, -loglambdasd)

tbl_estimates = full_join(sample_sizes,
                          df_pred,
                          by = join_by(phy, lighttreatment, lightintensity, leaftype)) |>
  select(accession,
         lighttreatment,
         lightintensity,
         leaftype,
         n_plant,
         where(is.numeric)) |>
  pivot_longer(
    cols = -c(accession, lighttreatment, lightintensity, leaftype, n_plant),
    names_to = c("trait", "statistic"),
    names_sep = "_",
    values_to = "value"
  ) |>
  pivot_wider(names_from = statistic, values_from = value)


# Write tbl-estimates-accession.csv
write_csv(tbl_estimates, "tables/tbl-estimates-accession.csv")

# Write tbl-estimates-accession-dictionary.csv

# ---- build dictionary ----
dict = tibble(
  variable = names(tbl_estimates),
  data_type = map_chr(tbl_estimates, type_label),
  acceptable_values = map_chr(tbl_estimates, acceptable_values),
) |>
  full_join(tibble(
    variable = c(
      "accession",
      "lighttreatment",
      "leaftype",
      "lightintensity",
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
write_csv(dict, "tables/tbl-estimates-accession-dictionary.csv")
