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

# Helper function to prepare data for estimating effects of guard cell size and
# f_gmax on tau
prepare_tau_anatomy_data = function(joined_summary, logtau_threshold) {
  n_remove = sum(joined_summary$logtau_mean >= logtau_threshold, na.rm = TRUE)
  out = joined_summary |>
    filter(logtau_mean < logtau_threshold) |>
    mutate(loggcl = log(guard_cell_length_um),
           logitfgmax = qlogis(f_gmax)) |>
    rename(
      logtaumean = logtau_mean,
      logtausd = logtau_sd,
      loglambdamean = loglambda_mean,
      loglambdasd = loglambda_sd,
      accid = acc_id,
      lightintensity = light_intensity,
      lighttreatment = light_treatment
    ) |>
    set_attr("n_removed", n_remove)
  
  out
}

# Helpers for writing data dictionaries. From ChatGPT5.2

# ---- helper: infer "acceptable values" for categorical-like columns ----
acceptable_values = function(x, max_unique = 30) {
  if (is.list(x)) return(NA_character_)
  ux <- sort(unique(na.omit(x)))
  n <- length(ux)
  
  if (inherits(x, "factor") || inherits(x, "character")) {
    if (n == 0) return(NA_character_)
    if (n <= max_unique) return(paste(ux, collapse = "; "))
    return(paste0(n, " unique values (too many to list)"))
  }
  
  # For numeric/integer/logical: usually not enumerated
  NA_character_
}

# ---- helper: pretty type label ----
type_label = function(x) {
  if (inherits(x, "Date")) return("Date")
  if (inherits(x, "POSIXct")) return("POSIXct")
  if (inherits(x, "POSIXlt")) return("POSIXlt")
  if (inherits(x, "factor")) return("factor")
  if (inherits(x, "character")) return("character")
  if (inherits(x, "integer")) return("integer")
  if (inherits(x, "numeric")) return("numeric")
  if (inherits(x, "logical")) return("logical")
  paste(class(x), collapse = "/")
}

# Function to estimate phylogenetic h2 from fitted models
get_phy_h2 = function(fit) {
  vars = crossing(
    resp = fit$formula$responses,
    term = c("sd_phy__", "sd_accession__", "sigma_")
  ) |>
    mutate(var = paste0(term, resp, if_else(term == "sigma_", "", "_Intercept")), .keep = "unused") |>
    pull(var)
  
  fit |>
    as_draws_df() |>
    select(starts_with("."), all_of(vars)) |>
    rename_with(\(.x) {
      str_remove(.x, "_Intercept")
    }, .cols = ends_with("_Intercept")) |>
    rename_with(\(.x) {
      str_replace(.x, "sigma_", "sd_resid__")
    }, .cols = starts_with("sigma_")) |>
    pivot_longer(
      cols = -starts_with("."),
      names_to = c("component", "resp"),
      names_sep = "__",
      values_to = "sd"
    ) |>
    mutate(
      var = sd^2,
      component = str_remove(component, "^sd_"),
      .keep = "unused"
    ) |>
    pivot_wider(names_from = component, values_from = var) |>
    mutate(h2 = phy / (phy + accession + resid), .keep = "unused") |>
    summarize(
      estimate = median(h2),
      lowerCI = quantile(h2, 0.025),
      upperCI = quantile(h2, 0.975),
      .by = resp
    )
  
}

# Function from ChatGPT to calculate precision matrices from brms phylogenetic covariance parameters
make_precision_phy <- function(draws_df) {
  # Variable order (match your column names)
  vars <- c(
    "loglambdamean_Intercept",
    "logtaumean_Intercept",
    "loggcl_Intercept",
    "logitfgmax_Intercept"
  )
  
  # SD columns
  sd_cols <- paste0("sd_phy__", vars)
  
  # Helper to find the right cor column regardless of ordering in the name
  cor_col <- function(a, b) {
    nm1 <- paste0("cor_phy__", a, "__", b)
    nm2 <- paste0("cor_phy__", b, "__", a)
    if (nm1 %in% names(draws_df))
      return(nm1)
    if (nm2 %in% names(draws_df))
      return(nm2)
    stop("Missing correlation column for: ", a, " and ", b)
  }
  
  # Extract SDs as matrix: draws x 4
  S <- as.matrix(draws_df[, sd_cols, drop = FALSE])
  Dn <- nrow(S)
  P  <- length(vars)
  
  # Build correlation array R: draws x 4 x 4
  R <- array(0, dim = c(Dn, P, P))
  for (d in seq_len(Dn))
    diag(R[d, , ]) <- 1
  
  for (i in 1:(P - 1)) {
    for (j in (i + 1):P) {
      cc <- cor_col(vars[i], vars[j])
      R[, i, j] <- draws_df[[cc]]
      R[, j, i] <- draws_df[[cc]]
    }
  }
  
  # Build covariance array Sigma and precision array Omega
  Sigma <- array(NA_real_, dim = c(Dn, P, P))
  Omega <- array(NA_real_, dim = c(Dn, P, P))
  
  for (d in seq_len(Dn)) {
    Dmat <- diag(S[d, ], nrow = P)
    Sigma_d <- Dmat %*% R[d, , ] %*% Dmat
    Sigma[d, , ] <- Sigma_d
    Omega[d, , ] <- solve(Sigma_d)
  }
  
  dimnames(Omega) <- list(NULL, vars, vars)
  # dimnames(Sigma) <- list(NULL, vars, vars)
  Omega
}

# Function from ChatGPT to plot ellipses
ellipse_points = function(mu, Sigma, level = 0.95, n = 200) {
  stopifnot(length(mu) == 2, all(dim(Sigma) == c(2, 2)))
  r <- sqrt(qchisq(level, df = 2))           # radius for chosen level
  theta <- seq(0, 2*pi, length.out = n)
  
  # unit circle
  circle <- rbind(cos(theta), sin(theta))
  
  # transform circle -> ellipse: mu + A %*% circle, where A A^T = Sigma
  A <- chol(Sigma)                           # upper-triangular
  pts <- t(circle) %*% A                     # (n x 2)
  data.frame(x = mu[1] + pts[,1], y = mu[2] + pts[,2])
}

# Saturated vapor pressure
svp = function(T_leaf, Pa) (0.61365 * exp(17.502 * T_leaf / (240.97 + T_leaf) / Pa))

# Convert log-change to %-change
log_to_percent = function(.x) {
  100 * (exp(.x) - 1)
}
