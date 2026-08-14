# --- Fit the phenomenological Weibull model to the synthetic aperture-
# kinetics data (r/41_simulate-aperture-kinetics.R) and ask how the
# resulting tau_hat and lambda_hat covary with gi and gmax, under the
# hypothesis that every curve shares the same intrinsic pore-closure rate
# constant k (i.e. no true dependence of closure kinetics on anatomy or
# aperture state).

source("r/header.R")

sim_aperture_data <- read_rds("objects/sim-aperture-data.rds")
sim_pars <- read_rds("objects/sim-aperture-pars.rds")$sim_pars

plan(multisession, workers = 9)

fits <- sim_aperture_data |>
  group_split(id) |>
  future_map_dfr(
    \(.df) bind_cols(id = .df$id[1], fit_nls_one(.df)),
    .progress = TRUE
  )

tbl_aperture_kinetics_fits <- fits |>
  left_join(sim_pars, by = "id") |>
  filter(converged) |>
  mutate(
    fgmax = gi / gmax,
    logtau_hat = log(tau_hat),
    loglambda_hat = log(lambda_hat),
    loggi = log(gi),
    loggmax = log(gmax),
    logfgmax = log(fgmax)
  )

write_rds(tbl_aperture_kinetics_fits, "objects/tbl-aperture-kinetics-fits.rds")

# --- Summary: predicted covariation of tau_hat, lambda_hat with gi and gmax

cor_tau_gi     <- cor.test(tbl_aperture_kinetics_fits$loggi, tbl_aperture_kinetics_fits$logtau_hat)
cor_tau_gmax   <- cor.test(tbl_aperture_kinetics_fits$loggmax, tbl_aperture_kinetics_fits$logtau_hat)
cor_tau_fgmax  <- cor.test(tbl_aperture_kinetics_fits$logfgmax, tbl_aperture_kinetics_fits$logtau_hat)
cor_lambda_gi    <- cor.test(tbl_aperture_kinetics_fits$loggi, tbl_aperture_kinetics_fits$loglambda_hat)
cor_lambda_gmax  <- cor.test(tbl_aperture_kinetics_fits$loggmax, tbl_aperture_kinetics_fits$loglambda_hat)
cor_lambda_fgmax <- cor.test(tbl_aperture_kinetics_fits$logfgmax, tbl_aperture_kinetics_fits$loglambda_hat)

tbl_aperture_kinetics_cor <- tibble(
  response = rep(c("log(tau_hat)", "log(lambda_hat)"), each = 3),
  predictor = rep(c("log(gi)", "log(gmax)", "log(f_gmax)"), 2),
  cor = c(cor_tau_gi$estimate, cor_tau_gmax$estimate, cor_tau_fgmax$estimate,
          cor_lambda_gi$estimate, cor_lambda_gmax$estimate, cor_lambda_fgmax$estimate),
  cor_lower = c(cor_tau_gi$conf.int[1], cor_tau_gmax$conf.int[1], cor_tau_fgmax$conf.int[1],
                cor_lambda_gi$conf.int[1], cor_lambda_gmax$conf.int[1], cor_lambda_fgmax$conf.int[1]),
  cor_upper = c(cor_tau_gi$conf.int[2], cor_tau_gmax$conf.int[2], cor_tau_fgmax$conf.int[2],
                cor_lambda_gi$conf.int[2], cor_lambda_gmax$conf.int[2], cor_lambda_fgmax$conf.int[2])
)

write_rds(tbl_aperture_kinetics_cor, "objects/tbl-aperture-kinetics-cor.rds")
