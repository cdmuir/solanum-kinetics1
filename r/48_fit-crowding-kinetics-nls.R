# --- Fit the Weibull model to the density-dependent crowding simulation
# (r/47) and ask how tau_hat/lambda_hat covary with gi and gmax at each
# crowding strength (kappa0), including PARTIAL correlations (each
# predictor controlling for the other), since the hypothesis under test is
# specifically about tau decoupling from gmax while remaining correlated
# with gi -- something a simple pairwise correlation can miss if gi and
# gmax are themselves correlated (they are, via f_gmax).

source("r/header.R")

# Partial correlation of x with y, controlling for z, via residuals of two
# linear regressions (avoids adding a new package dependency for this one
# quick check).
partial_cor <- function(x, y, z) {
  rx <- residuals(lm(x ~ z))
  ry <- residuals(lm(y ~ z))
  cor(rx, ry)
}

sim_crowding_data <- read_rds("objects/sim-crowding-aperture-data.rds")
sim_pars <- read_rds("objects/sim-crowding-aperture-pars.rds")$sim_pars

plan(multisession, workers = 9)

fits <- sim_crowding_data |>
  group_split(id, kappa0) |>
  future_map_dfr(
    \(.df) bind_cols(id = .df$id[1], kappa0 = .df$kappa0[1], fit_nls_one(.df)),
    .progress = TRUE
  )

tbl_crowding_kinetics_fits <- fits |>
  left_join(sim_pars |> select(id, kappa0, kappa, gi, gf, gmax, d_m2, s_m2),
            by = c("id", "kappa0")) |>
  filter(converged) |>
  mutate(
    fgmax = gi / gmax,
    logtau_hat = log(tau_hat),
    loglambda_hat = log(lambda_hat),
    loggi = log(gi),
    loggmax = log(gmax),
    logfgmax = log(fgmax)
  )

write_rds(tbl_crowding_kinetics_fits, "objects/tbl-crowding-kinetics-fits.rds")

cor_one_kappa0 <- function(.df) {
  cor_tau_gi    <- cor.test(.df$loggi, .df$logtau_hat)
  cor_tau_gmax  <- cor.test(.df$loggmax, .df$logtau_hat)
  cor_lam_gi    <- cor.test(.df$loggi, .df$loglambda_hat)
  cor_lam_gmax  <- cor.test(.df$loggmax, .df$loglambda_hat)

  # Partial correlations: gi controlling for gmax, and gmax controlling for gi
  pc_tau_gi   <- partial_cor(.df$loggi, .df$logtau_hat, .df$loggmax)
  pc_tau_gmax <- partial_cor(.df$loggmax, .df$logtau_hat, .df$loggi)
  pc_lam_gi   <- partial_cor(.df$loggi, .df$loglambda_hat, .df$loggmax)
  pc_lam_gmax <- partial_cor(.df$loggmax, .df$loglambda_hat, .df$loggi)

  tibble(
    response = rep(c("log(tau_hat)", "log(lambda_hat)"), each = 2),
    predictor = rep(c("log(gi)", "log(gmax)"), 2),
    cor = c(cor_tau_gi$estimate, cor_tau_gmax$estimate,
            cor_lam_gi$estimate, cor_lam_gmax$estimate),
    cor_lower = c(cor_tau_gi$conf.int[1], cor_tau_gmax$conf.int[1],
                  cor_lam_gi$conf.int[1], cor_lam_gmax$conf.int[1]),
    cor_upper = c(cor_tau_gi$conf.int[2], cor_tau_gmax$conf.int[2],
                  cor_lam_gi$conf.int[2], cor_lam_gmax$conf.int[2]),
    partial_cor = c(pc_tau_gi, pc_tau_gmax, pc_lam_gi, pc_lam_gmax)
  )
}

tbl_crowding_kinetics_cor <- tbl_crowding_kinetics_fits |>
  group_by(kappa0) |>
  group_modify(~ cor_one_kappa0(.x)) |>
  ungroup()

write_rds(tbl_crowding_kinetics_cor, "objects/tbl-crowding-kinetics-cor.rds")

# Side check: does crowding actually produce the diminishing-returns
# relationship between density and gmax implied by the hypothesis (i.e.
# does high-kappa0 compress the achievable gmax range at high density,
# relative to the pore-geometry-only case)?
tbl_crowding_gmax_check <- sim_pars |>
  mutate(
    gmax_crowded = gleaf_op_nonlinear(1, d_m2, s_m2, gec, x_prime, kappa)
  )

write_rds(tbl_crowding_gmax_check, "objects/tbl-crowding-gmax-check.rds")
