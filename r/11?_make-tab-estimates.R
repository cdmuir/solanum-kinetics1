# Make table of parameter estimates for each accession
source("r/header.R")

phy = read_rds("data/phylogeny.rds")
A = vcv(phy, corr = TRUE)

fit_amphi_high = read_rds("objects/fit_amphi_high.rds")
fit_amphi_low = read_rds("objects/fit_amphi_low.rds")
fit_pseudohypo_high = read_rds("objects/fit_pseudohypo_high.rds")
fit_pseudohypo_low = read_rds("objects/fit_pseudohypo_low.rds")

# add assertion - waiting to figure out how many iterations are needed for convergence
# check_convergence(fit_amphi_high, convergence_criteria)
# check_convergence(fit_amphi_low, convergence_criteria)
# check_convergence(fit_pseudohypo_high, convergence_criteria)
# check_convergence(fit_pseudohypo_low, convergence_criteria)
fit = fit_pseudohypo_low
diag = c(
  fit |>
    as_draws_df() |>
    summarise_draws() |>
    as_tibble() |>
    filter(variable != "lprior") |>
    summarize(
      rhat_max = max(rhat, na.rm = TRUE),
      ess_min = min(ess_bulk, na.rm = TRUE)
    ) |>
    as.list(),
  n_divergent = nuts_params(fit) |>
    subset(Parameter == "divergent__") |>
    pull(Value) |>
    sum()
)
diag



df_new = crossing(
  lighttreatment = unique(fit_amphi_high$data$lighttreatment),
  phy = unique(fit_amphi_high$data$phy)
) |>
  mutate(variable = paste0("...", row_number()))

tmp = posterior_epred(fit_amphi_high, newdata = df_new, resp = "logfgmax") |>
  as_draws_df() |>
  summarise_draws() |>
  select(variable, log_fgmax = median, q5, q95) |>
  full_join(df_new, by = "variable")
dim(tmp)
tmp
