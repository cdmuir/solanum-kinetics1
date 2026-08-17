# --- Fit the phenomenological Weibull model to the synthetic nonlinear-
# aperture-kinetics data (r/45_simulate-nonlinear-aperture-kinetics.R) and
# ask how the resulting tau_hat and lambda_hat covary with gi and gmax at
# each level of aperture-conductance nonlinearity (kappa), under the
# hypothesis that every curve shares the same intrinsic pore-closure rate
# constant k. kappa = 0 reproduces r/41-42 exactly (pore-geometry-only
# model); larger kappa imposes increasingly strong saturation of
# conductance relative to aperture at large aperture.

source("r/header.R")

sim_nonlinear_data <- read_rds("objects/sim-nonlinear-aperture-data.rds")
sim_pars <- read_rds("objects/sim-nonlinear-aperture-pars.rds")$sim_pars

plan(multisession, workers = 9)

fits <- sim_nonlinear_data |>
  group_split(id, kappa) |>
  future_map_dfr(
    \(.df) bind_cols(id = .df$id[1], kappa = .df$kappa[1], fit_nls_one(.df)),
    .progress = TRUE
  )

tbl_nonlinear_kinetics_fits <- fits |>
  left_join(sim_pars |> select(id, kappa, gi, gf, gmax, amax_um, d_m2, s_m2),
            by = c("id", "kappa")) |>
  filter(converged) |>
  mutate(
    fgmax = gi / gmax,
    logtau_hat = log(tau_hat),
    loglambda_hat = log(lambda_hat),
    loggi = log(gi),
    loggmax = log(gmax),
    logfgmax = log(fgmax)
  )

write_rds(tbl_nonlinear_kinetics_fits, "objects/tbl-nonlinear-kinetics-fits.rds")

# --- Summary: predicted covariation of tau_hat, lambda_hat with gi and
# gmax, separately at each kappa (degree of aperture-conductance
# nonlinearity)

cor_one_kappa <- function(.df) {
  cor_tau_gi     <- cor.test(.df$loggi, .df$logtau_hat)
  cor_tau_gmax   <- cor.test(.df$loggmax, .df$logtau_hat)
  cor_tau_fgmax  <- cor.test(.df$logfgmax, .df$logtau_hat)
  cor_lambda_gi    <- cor.test(.df$loggi, .df$loglambda_hat)
  cor_lambda_gmax  <- cor.test(.df$loggmax, .df$loglambda_hat)
  cor_lambda_fgmax <- cor.test(.df$logfgmax, .df$loglambda_hat)

  tibble(
    response = rep(c("log(tau_hat)", "log(lambda_hat)"), each = 3),
    predictor = rep(c("log(gi)", "log(gmax)", "log(f_gmax)"), 2),
    cor = c(cor_tau_gi$estimate, cor_tau_gmax$estimate, cor_tau_fgmax$estimate,
            cor_lambda_gi$estimate, cor_lambda_gmax$estimate, cor_lambda_fgmax$estimate),
    cor_lower = c(cor_tau_gi$conf.int[1], cor_tau_gmax$conf.int[1], cor_tau_fgmax$conf.int[1],
                  cor_lambda_gi$conf.int[1], cor_lambda_gmax$conf.int[1], cor_lambda_fgmax$conf.int[1]),
    cor_upper = c(cor_tau_gi$conf.int[2], cor_tau_gmax$conf.int[2], cor_tau_fgmax$conf.int[2],
                  cor_lambda_gi$conf.int[2], cor_lambda_gmax$conf.int[2], cor_lambda_fgmax$conf.int[2])
  )
}

tbl_nonlinear_kinetics_cor <- tbl_nonlinear_kinetics_fits |>
  group_by(kappa) |>
  group_modify(~ cor_one_kappa(.x)) |>
  ungroup()

write_rds(tbl_nonlinear_kinetics_cor, "objects/tbl-nonlinear-kinetics-cor.rds")
