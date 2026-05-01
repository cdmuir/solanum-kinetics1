# Make table of parameter estimates for each accession
source("r/header.R")

fit = read_rds("objects/best_model.rds")

tbl_estimates = fit$data |>
  select(
    accid,
    lighttreatment,
    lightintensity,
    leaftype,
    loggcl,
    logitfgmax,
    loglambdamean,
    loglambdasd,
    logtaumean,
    logtausd
  ) |>
  separate_wider_delim(
    cols = accid,
    delim = "-",
    names = c("accession", "id")
  )

write_csv(tbl_estimates, "tables/tbl-estimates-curve.csv")

# Write tbl-estimates-curve-dictionary.csv

# ---- build dictionary ----
dict = tibble(
  variable = names(tbl_estimates),
  data_type = map_chr(tbl_estimates, type_label),
  acceptable_values = map_chr(tbl_estimates, acceptable_values),
) |>
  mutate(acceptable_values = case_when(
    variable == "id" ~ "A-Z",
    TRUE ~ acceptable_values
  )) |>
  full_join(tibble(
    variable = c(
      "accession",
      "id",
      "lighttreatment",
      "lightintensity",
      "leaftype",
      "loggcl",
      "logitfgmax",
      "loglambdamean",
      "loglambdasd",
      "logtaumean",
      "logtausd"
    ),
    description = c(
      "TGRC accession",
      "Unique identifier for each accession",
      "Growth light intensity treatment",
      "Measurement light intensity treatment",
      "Leaf type treatment",
      "Estimate of log-transformed guard cell length (um)",
      "Estimate of logit-transformed stomatal conductance as a fraction of maximum conductance",
      "Point estimate of logit-transformed lag-time parameter",
      "Standard error of logit-transformed lag-time parameter",
      "Point estimate of log-transformed time-constant parameter",
      "Standard error of log-transformed time-constant parameter"
    ),
  ), by = join_by(variable)) |>
  mutate(
    acceptable_values = na_if(acceptable_values, ""),
    acceptable_values = replace_na(acceptable_values, "")
  )

# ---- write to CSV ----
write_csv(dict, "tables/tbl-estimates-curve-dictionary.csv")
