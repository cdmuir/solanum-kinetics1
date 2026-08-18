# --- Report Pareto k / PSIS-ESS / MCSE diagnostics for the plausible models ---
#
# Diagnostics for the plausible models identified in r/11_compare-models.R
#
# NOTE: the 64 mv-brms model fits live at "objects/fits/fit_NN.rds", as
# referenced by r/12_select-model.R and r/22_plot-collinear.R; this script
# reads from the same location.

source("r/header.R")

df_comparison = read_rds("objects/tbl-comparison.rds")

plausible_models = df_comparison |>
  filter(plausible) |>
  pull(model)

plausible_fits = tibble(model = plausible_models) |>
  mutate(
    file = str_replace(model, "^models/fit_", "objects/fits/fit_") |> paste0(".rds"),
    fit = map(file, read_rds),
    loo_obj = map(fit, \(.x) .x$criteria$loo)
  )

loo_obj = plausible_fits$loo_obj[[1]]

summarize_loo = function(loo_obj) {
  pk = pareto_k_table(loo_obj)
  n_curves = sum(pk[, "Count"])

  get_pct = function(row) {
    if (row %in% rownames(pk)) 100 * pk[row, "Proportion"] else 0
  }

  tibble(
    n_curves = n_curves,
    pct_good = get_pct("(-Inf, 0.7]"),
    pct_bad = get_pct("(0.7, 1]"),
    pct_very_bad = get_pct("(1, Inf)"),
    min_psis_ess = min(loo_obj$diagnostics$n_eff),
    # threshold = 1 (rather than loo's stricter sample-size-dependent default,
    # ~0.7): only Pareto k >= 1 ("very bad") make the MCSE estimate itself
    # undefined; k in (0.7, 1] ("bad") biases the underlying PSIS estimate
    # but the MCSE calculation remains defined.
    mcse_elpd_loo = mcse_loo(loo_obj, threshold = 1)
  )
}

tbl_loo_diagnostics = plausible_fits |>
  mutate(diag = map(loo_obj, summarize_loo)) |>
  select(model, diag) |>
  unnest(diag) |>
  mutate(model = str_remove(model, "^models/")) |>
  left_join(
    df_comparison |>
      mutate(model = str_remove(model, "^models/")) |>
      select(model, `$\\Delta \\mathrm{LOOIC}$`, SE, selected),
    by = "model"
  ) |>
  arrange(desc(selected), model) |>
  select(
    model,
    `$\\Delta \\mathrm{LOOIC}$`,
    SE,
    n_curves,
    pct_good,
    pct_bad,
    pct_very_bad,
    min_psis_ess,
    mcse_elpd_loo,
    selected
  )

write_rds(tbl_loo_diagnostics, "objects/tbl-loo-diagnostics.rds")
