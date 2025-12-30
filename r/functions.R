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

# Functions to calculate and summarize partial correlations from multivariate brms models

## Get correlation matrix from posterior
get_R = function(cor, p1, p2) {
  # cor is vector of correlations
  # p1 and p2 are vectors of matrix indices corresponding to cor
  n = length(cor)
  # rank ^2 - rank - 2 * n = 0
  rank = (1 + sqrt(1^2 - 4 * 1 * -2 * n)) / (2 * 1)
  R = diag(rank)
  for (i in seq_along(cor)) {
    R[p1[i], p2[i]] = cor[i]
    R[p2[i], p1[i]] = cor[i]
  }
  R
}

## Get partial correlations from off-diagonal elements of scaled inverse correlation matrix
get_parcor = function(Psi) {
  n = nrow(Psi)
  which(upper.tri(matrix(NA, n, n)), arr.ind = TRUE) |>
    as_tibble() |>
    mutate(parcor = map2_dbl(row, col, ~ -Psi[.x, .y] / sqrt(Psi[.x, .x] * Psi[.y, .y]))) |>
    rename(p1 = row, p2 = col)
}

## Summarize the partial correlations from a brms fit
summarize_parcor = function(fit) {
  
  resp_vars = formula(fit)$responses
  
  fit |>
    as_draws_df() |>
    select(starts_with("."),
           starts_with("cor_"),
           starts_with("rescor_")) |>
    rename_with(.fn = \(.x) {
      str_replace(.x, "rescor", "cor_resid")
    }, .cols = starts_with("rescor")) |>
    rename_with(.fn = \(.x) {
      str_remove_all(.x, "_Intercept")
    }, .cols = contains("_Intercept")) |>
    pivot_longer(
      cols = -starts_with("."),
      names_sep = "__",
      names_to = c("component", "trait1", "trait2"),
      values_to = "cor"
    ) |>
    mutate(
      component = str_remove(component, "cor_"),
      p1 = as.numeric(factor(trait1, levels = resp_vars)),
      p2 = as.numeric(factor(trait2, levels = resp_vars))
    ) |>
    summarize(R = list(get_R(cor, p1, p2)), .by = c(".draw", "component")) |>
    mutate(
      R_inv = map(R, ~ solve(.x)),
      D = map(R_inv, ~ diag(1 / sqrt(diag(
        .x
      )))),
      Psi = map2(R_inv, D, ~ .y %*% .x %*% .y)
    ) |>
    reframe(map_dfr(Psi, get_parcor, .progress = TRUE),
            .by = c(".draw", "component")) |>
    mutate(trait1 = resp_vars[p1], trait2 = resp_vars[p2]) |>
    select(-matches("^p[0-9]+$")) |>
    unite("pair", trait1, trait2, sep = "_") |>
    split(~ component + pair) |>
    map(summarize_draws) |>
    map(filter, variable == "parcor") |>
    imap_dfr(\(.x, .y) {
      comp_pair = str_split(.y, "\\.", n = 2)[[1]]
      comp = comp_pair[1]
      pair = comp_pair[2]
      .x |>
        mutate(component = comp, pair = pair)
    }) |>
    separate_wider_delim(pair, "_", names = c("trait1", "trait2"))
  
}

# Function to set priors for inertia model
get_interia_prior = function(.dat) {
  
  c(
    set_prior(
      "normal(0.01, 1)",
      nlpar = "gmin",
      lb = 0,
      ub = min(.dat$gsw)
    ),
    set_prior(
      glue("normal({mu}, 1)", mu = min(.dat$gsw)),
      nlpar = "gstar",
      lb = min(.dat$gsw) * 0.95,
      ub = min(.dat$gsw) * 1.05
    ),
    set_prior(
      glue("normal({mu}, 1)", mu = max(.dat$gsw)),
      nlpar = "ginit",
      lb = max(.dat$gsw) * 0.9,
      ub = first(.dat$gmax)
    ),
    prior(
      normal(0, 100),
      nlpar = "ik",
      ub = 0
    ),
    prior(
      normal(log(300), 1),
      nlpar = "logtau",
      lb = log(10),
      ub = log(10000)
    )
  )
  
}
