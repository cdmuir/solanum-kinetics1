# Functions moved out of r/functions.R because they are not called anywhere
# in the active r/00-39 pipeline or ms/ms.qmd. Kept here for reference/reuse
# rather than deleted outright, mirroring how the aperture-kinetics
# functions were archived to r/archive/functions-aperture-kinetics.R.
# Scripts that use these (all themselves archived) assume these
# definitions are available via source("r/header.R") as they were
# originally; that is no longer true after this move, so those archived
# scripts would need this file explicitly sourced to run again.

# --- Phylogenetic partial correlations (superseded by make_precision_phy()
# in r/functions.R, used by r/14_get-partial-cor.R) -------------------------
#
# get_R()/get_parcor()/summarize_parcor() are an earlier, more general
# implementation of computing partial correlations from a multivariate
# brms fit's posterior (correlation matrix -> partial correlations from
# the scaled inverse) than the fixed-five-trait make_precision_phy()
# approach currently used. get_phy_h2() is a related helper for computing
# phylogenetic heritability directly from a fitted model's variance
# components (used by the archived r/archive/24?_get-phy-h2.R), rather
# than from the pre-computed variance-decomposition table used in the
# active r/20_plot-variance.R.

# Build a symmetric correlation matrix from a vector of pairwise
# correlations and their row/column indices
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

# Get partial correlations from off-diagonal elements of scaled inverse
# correlation matrix
get_parcor = function(Psi) {
  n = nrow(Psi)
  which(upper.tri(matrix(NA, n, n)), arr.ind = TRUE) |>
    as_tibble() |>
    mutate(parcor = map2_dbl(row, col, ~ -Psi[.x, .y] / sqrt(Psi[.x, .x] * Psi[.y, .y]))) |>
    rename(p1 = row, p2 = col)
}

# Summarize the partial correlations from a brms fit
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

# Estimate phylogenetic heritability (h2) directly from a fitted model's
# variance components, for each response variable in `fit`. Depends on
# clean_posterior_names() in r/functions.R.
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
    clean_posterior_names() |>
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

# --- Guard-cell illustration helpers ---------------------------------------
# (used by the archived r/archive/24_make-illustrations.R)

# Close a polygon of (x, y) coordinates by setting the last row equal to
# the first.
close_polygon <- function(x) {
  x[nrow(x), ] <- x[1, ]
  x
}

# Rescale a set of illustration coordinates to a common size (largest x
# or y range = 1).
rescale_illustration <- function(coords) {
  scale <- max(c(diff(range(coords$x)), diff(range(coords$y))))
  coords |>
    select(x, y, side, name) |>
    mutate(x = x / scale, y = y / scale)
}

# --- Inertia-model priors --------------------------------------------------
# (an alternative model of stomatal kinetics not part of the current
# r/00-39 pipeline; see e.g. the archived r/archive/04_fit-inertia.R,
# r/archive/05_refit-inertia.R)

# Set priors for the inertia model
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
