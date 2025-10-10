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

# Helper functions to make all possible model forms
# From ChatGPT
# vars: character vector of predictor names
# max_order: highest interaction order allowed (<= length(vars))
all_hierarchical_models <- function(vars, max_order = length(vars), response = "y") {
  stopifnot(max_order >= 1, max_order <= length(vars))
  
  # All terms up to max_order (no intercept in labels)
  full_form   <- as.formula(paste0("~ (", paste(vars, collapse = " + "), ")^", max_order))
  all_labels  <- attr(terms(full_form), "term.labels")  # e.g., "a", "a:b", "a:b:c", ...
  split_term  <- function(t) strsplit(t, ":", fixed = TRUE)[[1]]
  
  # For each term, precompute the set of lower-order terms required by hierarchy
  # (all non-empty proper subsets of its factors, including the main effects)
  subsets_required <- lapply(all_labels, function(t) {
    f <- split_term(t)
    if (length(f) == 1) return(character(0))  # main effects have no prerequisites
    # all non-empty proper subsets:
    req <- unlist(lapply(1:(length(f) - 1), function(k) {
      combn(f, k, FUN = function(x) paste(x, collapse = ":"), simplify = TRUE)
    }), use.names = FALSE)
    unique(req)
  })
  names(subsets_required) <- all_labels
  
  # Power set over all_labels, then filter by hierarchy
  n <- length(all_labels)
  out <- vector("list", 2^n)  # upper bound; we'll drop invalid later
  keep <- logical(2^n)
  idx <- 0L
  
  for (mask in 0:(2^n - 1L)) {
    # Select terms where the bit is on
    sel <- all_labels[as.logical(intToBits(mask))[seq_len(n)]]
    # hierarchy check: for every selected term, all its required subsets must be present
    ok <- TRUE
    for (t in sel) {
      req <- subsets_required[[t]]
      if (length(req) && !all(req %in% sel)) { ok <- FALSE; break }
    }
    if (!ok) next
    idx <- idx + 1L
    # Build a formula RHS; intercept is implicit, but we include 1 for clarity
    rhs <- if (length(sel) == 0) "1" else paste(c("1", sel), collapse = " + ")
    out[[idx]] <- paste(response, "~", rhs)
    keep[idx] <- TRUE
  }
  
  # Drop unused slots
  out <- out[keep]
  unlist(out, use.names = FALSE)
}
