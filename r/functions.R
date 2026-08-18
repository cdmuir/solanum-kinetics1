# Custom functions for the Solanum stomatal-kinetics analysis pipeline
# (r/00_ through r/38_) and ms/ms.qmd, grouped by purpose below. Note:
#   - The Ochoa et al. (2024) anatomical pore-conductance model, its
#     Kaiser (2009) series-resistance extension, and the phenomenological
#     nonlinear / Lehmann & Or (2015) crowding variants that were only
#     used by r/40-50 have been moved to
#     r/archive/functions-aperture-kinetics.R, alongside those archived
#     scripts.
#   - get_R(), get_parcor(), summarize_parcor(), get_phy_h2(),
#     close_polygon(), rescale_illustration(), and get_interia_prior(),
#     which were not called anywhere in the active pipeline, have been
#     moved to r/archive/functions-legacy-helpers.R.

# --- Recoding raw treatment/leaf-type codes -------------------------------
#
# Recode raw treatment/leaf-type codes into the factors used throughout
# downstream analyses, so the same labels, mapping direction, and factor
# levels are used everywhere these codes are recoded.
recode_lighttreatment = function(x) {
  factor(case_when(
    x == "high" ~ "sun",
    x == "low" ~ "shade"
  ), levels = c("shade", "sun"))
}

recode_lightintensity = function(x) {
  factor(case_when(
    x == "150" ~ "low",
    x == "2000" ~ "high"
  ), levels = c("low", "high"))
}

recode_leaftype = function(x) {
  factor(case_when(
    x == "1-sided RH" ~ "pseudohypo",
    x == "2-sided RH" ~ "amphi"
  ), levels = c("amphi", "pseudohypo"))
}

# --- Psychrometric / gas-exchange calculations ----------------------------

# Saturated vapor pressure (for RH calculation)
svp = function(T_leaf, Pa) (0.61365 * exp(17.502 * T_leaf / (240.97 + T_leaf) / Pa))

# Saturated vapor pressure (for VPD calculation)
li6800_svp = function(T_degreeC) {
  0.61365 * exp(17.502 * T_degreeC / (240.97 + T_degreeC))
}

# --- Weibull curve fitting (Bayesian pipeline: r/02-r/06, r/30) -----------
#
# Fit @eq-weibull to a single humidity-response curve with brms, check
# HMC convergence diagnostics against project-wide criteria, and re-fit
# with progressively more thinning/adapt_delta until a curve converges.

# Fit the Weibull closure-kinetics model to one curve via brms
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
    chains = 4,
    cores = 1,
    backend = "cmdstanr",
    control = list(adapt_delta = adapt_delta),
    seed = seed
  )
  
}

# Check a brms fit's Rhat, bulk ESS, and divergent-transition counts
# against the project's convergence_criteria (r/header.R)
check_convergence = function(fit, convergence_criteria) {
  diag = c(
    fit |>
      as_draws_df() |>
      summarize_draws() |>
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

# Re-fit a saved curve with more thinning/adapt_delta until it converges
# (or max_thin is reached), overwriting the file at `path`
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

# --- VPD trajectory fitting (r/27_summarize-vpd.R, r/29_plot-vpd-h2o-rate.R)
#
# Fit a saturating exponential to a curve-level trajectory (VPD or raw
# humidity channels), used to summarize how quickly the realized VPD
# stimulus approached its final value.

# Fit a saturating exponential to a single curve's trajectory of `yvar`,
# using data-driven starting values so nlsLM converges regardless of
# whether the trajectory rises (e.g. VPD, typically) or falls (e.g. H2O_s,
# typically) over the course of the curve.
fit_saturating = function(df, yvar) {
  y = df[[yvar]]
  t = df$t_sec
  
  R0_0 = y[which.min(t)]
  Rend = y[which.max(t)]
  rng = diff(range(y))
  if (rng == 0) rng = 1e-6 # avoid a degenerate zero-range start
  
  increasing = Rend >= R0_0
  Asym0 = if (increasing) max(y) + 0.01 * rng else min(y) - 0.01 * rng
  
  half_target = R0_0 + 0.63 * (Asym0 - R0_0)
  t_half = max(t[which.min(abs(y - half_target))], 1)
  lrc0 = log(1 / t_half)
  
  dat = tibble(y = y, t_sec = t)
  m = tryCatch(
    nlsLM(
      y ~ Asym + (R0 - Asym) * exp(-exp(lrc) * t_sec),
      data = dat,
      start = list(Asym = Asym0, R0 = R0_0, lrc = lrc0),
      control = nls.lm.control(maxiter = 200)
    ),
    error = function(e) NULL
  )
  
  if (is.null(m)) {
    return(tibble(
      converged = FALSE,
      r2_sat = NA_real_,
      sat_asym = NA_real_,
      sat_r0 = NA_real_,
      sat_rate = NA_real_,
      sat_rate_se = NA_real_
    ))
  }
  
  r2 = 1 - sum(residuals(m) ^ 2) / sum((y - mean(y)) ^ 2)
  cf = summary(m)$coefficients
  lrc_est = cf["lrc", "Estimate"]
  lrc_se = cf["lrc", "Std. Error"]
  rate_est = exp(lrc_est)
  # Delta method: rate = exp(lrc), so d(rate)/d(lrc) = rate.
  rate_se = rate_est * lrc_se
  
  tibble(
    converged = TRUE,
    r2_sat = r2,
    sat_asym = coef(m)[["Asym"]],
    sat_r0 = coef(m)[["R0"]],
    # rate constant (1/s): larger = faster approach to asymptote
    sat_rate = rate_est,
    # SE of sat_rate via the delta method, to propagate measurement error in
    # subsequent analyses (e.g., as `se(sat_rate) ~ 1` in a brms me() term)
    sat_rate_se = rate_se
  )
}

# --- Data preparation for the multiresponse model (r/10_fit-all.R) --------

# Rename/derive columns of the joined per-curve summary table (data/
# joined-summary.rds) into the names expected by the multiresponse brms
# formulas, and drop curves with implausibly large tau (see logtau_threshold
# in r/header.R), recording how many were removed as an attribute for
# reporting in the manuscript.
prepare_tau_anatomy_data = function(joined_summary, logtau_threshold) {
  n_remove = sum(joined_summary$logtau_mean >= logtau_threshold, na.rm = TRUE)
  out = joined_summary |>
    filter(logtau_mean < logtau_threshold) |>
    mutate(loggcl = log(guard_cell_length_um),
           loggi = log(ginit_mean), 
           loggmax = log(gmax)) |>
    rename(
      logtaumean = logtau_mean,
      logtausd = logtau_sd,
      loglambdamean = loglambda_mean,
      loglambdasd = loglambda_sd,
      accid = acc_id,
      leaftype = leaf_type,
      lightintensity = light_intensity,
      lighttreatment = light_treatment
    ) |>
    set_attr("n_removed", n_remove)
  
  out
}

# --- Posterior prediction & parameter extraction --------------------------

# Get posterior-predicted trait values (back-transformed via `inv`) for
# each accession/treatment combination in `newdata`, marginalizing over
# individual-level random effects but retaining phylogenetic and
# population random effects (r/16_make-tbl-estimates-accession.R)
get_posterior_epred = function(fit, newdata, resp, prefix = resp, inv) {
  posterior_epred(
    fit,
    newdata = newdata,
    resp = resp,
    re_formula = ~ (1 | accession) +
      (1 | gr(phy, cov = A))
  ) |>
    as_draws_df() |>
    summarize_draws(median, quantile2, .args = list(probs = c(0.025, 0.975))) |>
    mutate(
      !!paste0(prefix, "_estimate") := inv(median),
      !!paste0(prefix, "_lowerCI")  := inv(`q2.5`),
      !!paste0(prefix, "_upperCI")  := inv(`q97.5`)
    ) |>
    select(variable, starts_with(prefix))
}

# Extract a fixed set of parameters (`pars_of_interest`, defined by the
# calling script) from a brms fit's posterior and summarize each as a
# mean and 95% interval, tagged with `model_name` for cross-model
# comparison (r/32_plot-comparison.R)
extract_pars = function(fit, model_name) {
  as_draws_df(fit, variable = pars_of_interest) |>
    as_tibble() |>
    select(all_of(pars_of_interest)) |>
    pivot_longer(everything(), names_to = "parameter", values_to = "value") |>
    summarize(
      estimate = mean(value),
      lower = quantile(value, 0.025),
      upper = quantile(value, 0.975),
      .by = parameter
    ) |>
    mutate(model = model_name)
}

# --- Phylogenetic partial correlations from multivariate brms models -----
# (r/14_get-partial-cor.R)
#
# make_precision_phy() builds the phylogenetic covariance matrix for each
# posterior draw directly from a fixed set of trait SD/correlation
# columns, then inverts it to a precision matrix, from which
# r/14_get-partial-cor.R computes partial correlations. An earlier, more
# general implementation of the same idea (get_R()/get_parcor()/
# summarize_parcor(), plus the related get_phy_h2()) has been moved to
# r/archive/functions-legacy-helpers.R since it is not currently called
# anywhere in r/00-38 or ms.qmd.

# Function from ChatGPT to calculate precision matrices from brms
# phylogenetic covariance parameters, for a fixed set of five traits
# (lambda, tau, gcl, gi, gmax)
make_precision_phy <- function(draws_df) {
  # Variable order (match your column names)
  vars <- c(
    "loglambdamean_Intercept",
    "logtaumean_Intercept",
    "loggcl_Intercept",
    "loggi_Intercept",
    "loggmax_Intercept"
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

# --- Trait/parameter labeling for tables and figures ----------------------
# (r/13, r/16, r/18, r/20)

# Canonical LaTeX labels for model response variables (traits), used
# consistently across tables and figures.
trait_latex_label = function(x) {
  case_when(
    x == "logitfgmax" ~ "$\\mathrm{logit} \\left( f_\\mathrm{gmax} \\right)$",
    x == "loggcl" ~ "$\\log \\left( l_\\mathrm{gc} \\right)$",
    x == "loggi" ~ "$\\log \\left( g_\\mathrm{i} \\right)$",
    x == "loggmax" ~ "$\\log \\left( g_\\mathrm{max} \\right)$",
    x == "loglambdamean" ~ "$\\log \\left( \\lambda \\right)$",
    x == "logtaumean" ~ "$\\log \\left( \\tau \\right)$",
    TRUE ~ NA_character_
  )
}

# Same labels as trait_latex_label(), but matched against a raw brms
# parameter name that contains the trait code as a substring (e.g.
# "b_logtaumean_Intercept"), rather than requiring an exact match.
trait_latex_from_paramname = function(param) {
  case_when(
    str_detect(param, "_logitfgmax_") ~ trait_latex_label("logitfgmax"),
    str_detect(param, "_loggcl_") ~ trait_latex_label("loggcl"),
    str_detect(param, "_loggi_") ~ trait_latex_label("loggi"),
    str_detect(param, "_loggmax_") ~ trait_latex_label("loggmax"),
    str_detect(param, "_loglambdamean_") ~ trait_latex_label("loglambdamean"),
    str_detect(param, "_logtaumean_") ~ trait_latex_label("logtaumean"),
    TRUE ~ NA_character_
  )
}

# Canonical labels for variance-component / random-effect group codes,
# used consistently across the fixed/random-effects table
# (r/18_make-tbl-fit-summary.R) and the variance-decomposition figure
# (r/20_plot-variance.R). "accid" (the raw among-individual random
# effect SD reported in the table) and "ind" (the figure's combined
# among-individual + residual quantity, see Results) are both labeled
# "among-individual"; the manuscript text explains that the
# among-individual component in the variance-decomposition figure/table
# includes residual variance.
recode_variance_component = function(x) {
  case_when(
    x == "phy" ~ "phylogenetic",
    x == "accession" ~ "population (nonphylogenetic)",
    x %in% c("accid", "ind") ~ "among-individual",
    TRUE ~ NA_character_
  )
}

# Standardize brms posterior draw column names: strip the "_Intercept"
# suffix and rename the residual SD prefix ("sigma_") to the same
# "sd_<component>__<resp>" naming convention used for other variance
# components (e.g. "sd_phy__", "sd_accession__").
clean_posterior_names = function(.data) {
  .data |>
    rename_with(\(.x) str_remove(.x, "_Intercept"), .cols = contains("_Intercept")) |>
    rename_with(\(.x) str_replace(.x, "^sigma_", "sd_resid__"), .cols = starts_with("sigma_"))
}

# --- Data-dictionary helpers (r/15, r/16; from ChatGPT 5.2) ---------------
#
# Infer a human-readable "acceptable values" description and type label
# for each column of a table, used to auto-generate the data dictionaries
# accompanying tables/tbl-estimates-curve.csv and
# tables/tbl-estimates-accession.csv.

# Infer "acceptable values" for categorical-like columns
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

# Pretty type label
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

# --- Mediation-plot helpers (r/24_plot-mediation.R) -----------------------
#
# Prepare edge (path-coefficient) and node data frames, and join them by
# endpoint coordinates, for the hand-drawn path-analysis diagram in
# @fig-mediation.

# Prepare edge data: effect size, significance, sign, and a formatted
# label, and recode raw variable names to display labels
prepare_edges = function(.x) {
  .x |>
    mutate(
      abs_est = abs(estimate),
      effect_size = abs_est / se,
      sig = as.factor(ifelse(`q2.5` > 0 |
                               `q97.5` < 0, "sig.", "n.s.")),
      sign = ifelse(estimate >= 0, "positive", "negative"),
      label = sprintf("%.2f (%.2f,~%.2f)", estimate, `q2.5`, `q97.5`)
    ) |>
    mutate(across(
      c(from, to),
      \(.x) recode(
        .x,
        logtaumean = "$\\tau$",
        logitfgmax = "$f_\\mathrm{gmax}$",
        loggcl = "$l_\\mathrm{gc}$",
        loggi = "$g_\\mathrm{i}$",
        lighttreatmentsun = "sun treatment",
        lightintensityhigh = "high light",
        leaftypepseudohypo = "pseudohypo"
      )
    ))
}

# Join edges to node (x, y) coordinates by name, for both endpoints
join_nodes_edges = function(df_edges, df_nodes) {
  df_edges |>
    left_join(df_nodes, by = c("from" = "name")) |>
    rename(xstart = x, ystart = y) |>
    left_join(df_nodes, by = c("to" = "name")) |>
    rename(xend = x, yend = y) 
}

# --- Conceptual-figure kinetics helpers (r/26_plot-conceptual.R) ---------
#
# The Weibull closure equation and simple wrong-way/right-way response
# curves used to draw the schematic gsw(t) trace in @fig-conceptual.
# wwr_fun()/rwr_fun() read several parameters (g_pre, g_peak, t_wwr, gf0,
# tau_A, lambda_A) from the calling script's environment rather than as
# arguments. Also used by the archived r/archive/27_plot-conceptual-BD.R.

# Weibull closure-kinetics equation (@eq-weibull)
weibull_gsw <- function(t, gi, gf, tau, lambda) {
  gf + (gi - gf) * exp(-(t / tau)^lambda)
}

# Schematic wrong-way response: gsw rises from g_pre to g_peak
wwr_fun <- function(t) {
  g_pre + (g_peak - g_pre) * (1 - exp(-t / (t_wwr / 3)))
}

# Schematic right-way response: gsw decays from g_peak via weibull_gsw()
rwr_fun <- function(t) {
  weibull_gsw(t - t_wwr, gi = g_peak, gf = gf0, tau = tau_A, lambda = lambda_A)
}

# --- Plotting helpers (generic geometry) -----------------------------------
#
# Guard-cell illustration helpers (close_polygon(), rescale_illustration())
# have been moved to r/archive/functions-legacy-helpers.R, since they are
# only used by the archived r/archive/24_make-illustrations.R.

# Compute points along a 2D confidence ellipse for a bivariate normal
# with mean `mu` and covariance `Sigma`, for drawing phylogenetic
# covariance ellipses (r/21_plot-gcl-tau.R). Function from ChatGPT.
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

# --- Null-simulation helpers (r/33, r/34, r/35) ---------------------------
#
# Test whether the gi-tau association could be a curve-fitting artifact:
# simulate synthetic gsw(t) trajectories that combine one curve's real
# tau/lambda/sigma with another (randomly chosen) curve's real gi/gf, so
# there is no true association between gi and the kinetic parameters by
# construction, then re-fit and check whether a spurious correlation
# emerges.

# Simulate one synthetic curve: real time design/tau/lambda/sigma from
# curve `i`, but gi/gf swapped in from curve `j`
sim_one_dataset = function(real_pars, time_design, i, j) {
  pars_i = real_pars |>
    filter(id == i) |> 
    select(curve = id, tau, lambda, sigma) 
  pars_j = real_pars |>
    filter(id == j) |> 
    select(gi, gf)
  
  bind_cols(pars_i, pars_j) |>
    left_join(time_design, by = join_by(curve)) |>
    mutate(
      gsw_true = gf + (gi - gf) * exp(-(t_sec / tau) ^ lambda),
      gsw_sim = gsw_true + rnorm(n(), 0, sigma)
    )
  
}

# Fast re-fit of the Weibull functional form to one simulated curve via
# nls() (used in place of the full Bayesian pipeline for speed)
fit_nls_one = function(df) {
  gf0 = min(df$gsw_sim)
  gi0 = max(df$gsw_sim)
  m = tryCatch(
    nls(
      gsw_sim ~ gf + (gi - gf) * exp(-(t_sec / tau) ^ lambda),
      data = df,
      start = list(gf = gf0, gi = gi0, tau = 200, lambda = 1.2),
      control = nls.control(maxiter = 200, warnOnly = TRUE)
    ),
    error = function(e) NULL
  )
  if (is.null(m) || !m$convInfo$isConv) {
    return(tibble(converged = FALSE, gi_hat = NA_real_, gf_hat = NA_real_,
                  tau_hat = NA_real_, lambda_hat = NA_real_))
  }
  cf = coef(m)
  tibble(
    converged = TRUE,
    gi_hat = cf[["gi"]],
    gf_hat = cf[["gf"]],
    tau_hat = cf[["tau"]],
    lambda_hat = cf[["lambda"]]
  )
}

# One full null-simulation replicate: randomly pair curves, simulate,
# re-fit all of them, and compute the resulting gi-log(tau) correlation
run_one_replicate = function(replicate_id, real_pars, time_design, quiet = FALSE) {
  
  i = sample(real_pars$id, nrow(real_pars))
  j = sample(real_pars$id, nrow(real_pars))
  
  if (!quiet) {
    cli_alert_info(glue("Simulating {n} synthetic curves...", n = prettyNum(length(i), big.mark = ",")))
  }
  
  sim_data = map2(
    i,
    j,
    sim_one_dataset,
    real_pars = real_pars,
    time_design = time_design,
    .progress = !quiet
  )
  
  if (!quiet) {
    cli_alert_info(glue("Fitting {n} synthetic curves...", n = prettyNum(length(i), big.mark = ",")))
  }
  
  fits = sim_data |>
    map_dfr(fit_nls_one, .progress = !quiet)
  
  fits_valid = fits |>
    filter(converged, gi_hat > 0, gf_hat > 0, gf_hat < gi_hat, gi_hat < 1.4,
           tau_hat > 0, lambda_hat > 0) |>
    mutate(logtau_hat = log(tau_hat))
  
  ct = suppressWarnings(cor.test(fits_valid$gi_hat, fits_valid$logtau_hat))
  
  tibble(
    replicate = replicate_id,
    n_curves = nrow(fits_valid),
    n_converged = sum(fits$converged),
    n_total = nrow(fits),
    cor = unname(ct$estimate),
    cor_lower = ct$conf.int[1],
    cor_upper = ct$conf.int[2]
  )
  
}

# --- General utilities -----------------------------------------------------

# Convert log-change to %-change
log_to_percent = function(.x) {
  100 * (exp(.x) - 1)
}

# Spell out small integers (1-10) as words for use in running text; falls
# back to the numeral for anything else
num2word = function(x) {
  if (x == 1) {
    "one"
  } else if (x == 2) {
    "two"
  } else if (x == 3) {
    "three"
  } else if (x == 4) {
    "four"
  } else if (x == 5) {
    "five"
  } else if (x == 6) {
    "six"
  } else if (x == 7) {
    "seven"
  } else if (x == 8) {
    "eight"
  } else if (x == 9) {
    "nine"
  } else if (x == 10) {
    "ten"
  } else {
    as.character(x)
  }
}

