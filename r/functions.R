# Functions for fitting curves using brms
fit_rh1 = function(formula,
                   data,
                   prior,
                   thin,
                   adapt_delta,
                   seed) {
  brm(
    formula = formula,
    data = data,
    prior = prior,
    iter = thin * 2000,
    thin = thin,
    chains = 1,
    cores = 1,
    backend = "cmdstanr",
    control = list(adapt_delta = adapt_delta),
    seed = seed
  )
  
}

check_convergence = function(fit, convergence_criteria) {
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
  
  if (diag$rhat_max < convergence_criteria$rhat_max &
      diag$ess_min > convergence_criteria$ess_min &
      diag$n_divergent <= convergence_criteria$n_divergent) {
    return(TRUE)
  } else {
    return(FALSE)
  }
  
}

refit_rh = function(path, convergence_criteria, max_thin = 10) {
  fit = read_rds(path)
  conv = check_convergence(fit, convergence_criteria)
  
  while (!conv & fit$fit@sim$thin <= max_thin) {
    fit = fit_rh1(
      formula = fit$formula,
      data = fit$data,
      prior = fit$prior,
      thin = fit$fit@sim$thin + 1,
      adapt_delta = min(0.8 * 1.1^fit$fit@sim$thin, 0.99),
      seed = 360036340 + fit$fit@sim$thin
    )
    conv = check_convergence(fit, convergence_criteria)
  }
  write_rds(fit, path)
}