# --- Simulate the null correlation between gi and log(tau) -------------------

source("r/header.R")

n_replicates = 1e3
sim_seed = 886858562

pars_summary = read_rds("objects/pars-summary.rds")
joined_data = read_rds("data/joined-data.rds")

# --- Real per-curve building blocks: tau, lambda, sigma (residual SD), gi, gf

real_pars = pars_summary |>
  filter(
    variable %in% c(
      "b_logtau_Intercept",
      "b_loglambda_Intercept",
      "sigma",
      "ginit",
      "gfinal"
    )
  ) |>
  select(id, variable, mean) |>
  pivot_wider(names_from = variable, values_from = mean) |>
  transmute(
    id,
    tau = exp(b_logtau_Intercept),
    lambda = exp(b_loglambda_Intercept),
    sigma,
    gi = ginit,
    gf = gfinal
  )

time_design = joined_data |> select(curve, t_sec)

# --- Run all replicates ----------------------------------------------------

set.seed(sim_seed)
plan(multisession, workers = 19)

null_sim_results = future_map_dfr(seq_len(n_replicates), \(b) {
  run_one_replicate(b, real_pars, time_design, quiet = TRUE)
}, .progress = TRUE)

write_rds(null_sim_results, "objects/null-sim-fgmax-tau.rds")

# Save one example simulated dataset (true + estimated parameters) alongside
# its fits, for diagnostic/illustrative plotting.
set.seed(sim_seed)
example_sim_data = sim_one_dataset(real_pars,
                                   time_design,
                                   "LA1364-X_pseudohypo_2000",
                                   "LA0716-T_pseudohypo_2000")
example_fits = bind_cols(example_sim_data, fit_nls_one(example_sim_data))
