#!/usr/bin/env Rscript
# Fit ONE of the 64 models from r/10_fit-all.R on a CHTC execute node.
#
# Usage (as an HTCondor job): fit_model.R <index>
#   <index> is 0-63 (HTCondor's $(Process)) OR 1-64; both are accepted -- see
#   below. The index is mapped to a (lambda formula, tau formula) pair using
#   the SAME explicit scheme as chtc/reassemble_results.R:
#
#     lambda_idx <- (i - 1) %/% 8   # 0..7 -> bf_lambda0..bf_lambda7
#     tau_idx    <- (i - 1) %% 8    # 0..7 -> bf_tau0..bf_tau7
#
# This mirrors the 8x8 crossing of lambda/tau formulas in r/10_fit-all.R
# (lines 100-129), but computed explicitly rather than via tidyr::crossing()
# so the mapping is guaranteed stable across machines/package versions.
#
# This script intentionally does NOT source r/header.R: header.R calls
# rm(list = ls()) and loads ~30 packages (many for plotting/GIS) that are
# irrelevant to model fitting and would bloat the CHTC container. Only the
# packages actually needed to fit these brms models are loaded below.

suppressPackageStartupMessages({
  library(ape)
  library(brms)
  library(checkmate)
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(glue)
  library(readr)
  library(stringr)
  library(cmdstanr)
})

# The container's %environment sets CMDSTAN to the baked-in CmdStan
# installation (see chtc/image.def). Point cmdstanr at it explicitly rather
# than relying on cmdstanr's default (~/.cmdstan), which may not exist for
# the job's runtime user/home.
cmdstan_env <- Sys.getenv("CMDSTAN")
if (nzchar(cmdstan_env)) {
  cmdstanr::set_cmdstan_path(cmdstan_env)
}

# ---- CLI argument ----------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)
assert_true(length(args) == 1)

raw_index <- as.integer(args[[1]])
assert_true(!is.na(raw_index))

# HTCondor's $(Process) is always 0-based for `queue 64` (0..63); convert to
# the 1-based model index (1..64) used throughout this script.
i <- raw_index + 1L
assert_true(i >= 1L && i <= 64L)

lambda_idx <- (i - 1L) %/% 8L  # 0..7
tau_idx    <- (i - 1L) %% 8L   # 0..7

file_label <- str_pad(i, 2, "left", "0")
out_file   <- glue("fit_{file_label}.rds")

cli::cli_inform("Model index {i} (0-based Process {raw_index}): lambda_idx = {lambda_idx}, tau_idx = {tau_idx} -> {out_file}")

# ---- Minimal reimplementation of prepare_tau_anatomy_data() -----------
# (originally defined in r/functions.R; reproduced here to avoid pulling in
# the full functions.R / header.R dependency chain)

prepare_tau_anatomy_data <- function(joined_summary, logtau_threshold) {
  n_remove <- sum(joined_summary$logtau_mean >= logtau_threshold, na.rm = TRUE)
  out <- joined_summary |>
    filter(logtau_mean < logtau_threshold) |>
    mutate(
      loggcl = log(guard_cell_length_um),
      loggi = log(ginit_mean),
      loggmax = log(gmax)
    ) |>
    rename(
      logtaumean = logtau_mean,
      logtausd = logtau_sd,
      loglambdamean = loglambda_mean,
      loglambdasd = loglambda_sd,
      accid = acc_id,
      leaftype = leaf_type,
      lightintensity = light_intensity,
      lighttreatment = light_treatment
    )
  attr(out, "n_removed") <- n_remove
  out
}

logtau_threshold <- 7

# ---- Load data ----------------------------------------------------------

joined_summary <- readRDS("joined-summary.rds") |>
  prepare_tau_anatomy_data(logtau_threshold)

assert_true(all(!is.na(joined_summary$loglambdamean)))
assert_true(all(!is.na(joined_summary$loglambdasd)))
assert_true(all(!is.na(joined_summary$logtausd)))
assert_true(all(!is.na(joined_summary$loggcl)))
assert_true(all(!is.na(joined_summary$loggi)))
assert_true(all(!is.na(joined_summary$loggmax)))
assert_true(all(!is.na(joined_summary$lighttreatment)))
assert_true(all(!is.na(joined_summary$lightintensity)))
assert_true(all(!is.na(joined_summary$leaftype)))

phy <- readRDS("phylogeny.rds")
A <- vcv(phy, corr = TRUE)
thin <- 6

# ---- Formulas (mirrors r/10_fit-all.R lines 36-98 exactly) --------------

bf_lambda0 <- bf(
  loglambdamean | se(loglambdasd, sigma = TRUE) ~
    lighttreatment +
    lightintensity +
    leaftype +
    loggcl +
    loggi +
    loggmax +
    (1 | accid) +
    (1 | a | accession) +
    (1 | b | gr(phy, cov = A))
)

bf_lambda1 <- update(bf_lambda0, . ~ . - loggcl)
bf_lambda2 <- update(bf_lambda0, . ~ . - loggi)
bf_lambda3 <- update(bf_lambda0, . ~ . - loggmax)
bf_lambda4 <- update(bf_lambda0, . ~ . - loggcl - loggi)
bf_lambda5 <- update(bf_lambda0, . ~ . - loggcl - loggmax)
bf_lambda6 <- update(bf_lambda0, . ~ . - loggi - loggmax)
bf_lambda7 <- update(bf_lambda0, . ~ . - loggcl - loggi - loggmax)

bf_tau0 <- bf(
  logtaumean | se(logtausd, sigma = TRUE) ~
    lighttreatment +
    lightintensity +
    leaftype +
    loggcl +
    loggi +
    loggmax +
    (1 | accid) +
    (1 | a | accession) +
    (1 | b | gr(phy, cov = A))
)

bf_tau1 <- update(bf_tau0, . ~ . - loggcl)
bf_tau2 <- update(bf_tau0, . ~ . - loggi)
bf_tau3 <- update(bf_tau0, . ~ . - loggmax)
bf_tau4 <- update(bf_tau0, . ~ . - loggcl - loggi)
bf_tau5 <- update(bf_tau0, . ~ . - loggcl - loggmax)
bf_tau6 <- update(bf_tau0, . ~ . - loggi - loggmax)
bf_tau7 <- update(bf_tau0, . ~ . - loggcl - loggi - loggmax)

bf_gcl <- bf(loggcl ~
               lighttreatment +
               leaftype +
               (1 | a | accession) +
               (1 | b | gr(phy, cov = A)))

bf_gi <- bf(
  loggi ~
    lighttreatment +
    lightintensity +
    leaftype +
    (1 | accid) +
    (1 | a | accession) +
    (1 | b | gr(phy, cov = A))
)

bf_gmax <- bf(loggmax ~
                lighttreatment +
                leaftype +
                (1 | a | accession) +
                (1 | b | gr(phy, cov = A)))

bf_lambda_list <- list(bf_lambda0, bf_lambda1, bf_lambda2, bf_lambda3,
                        bf_lambda4, bf_lambda5, bf_lambda6, bf_lambda7)
bf_tau_list <- list(bf_tau0, bf_tau1, bf_tau2, bf_tau3,
                     bf_tau4, bf_tau5, bf_tau6, bf_tau7)

form <- bf_lambda_list[[lambda_idx + 1L]] + bf_tau_list[[tau_idx + 1L]] +
  bf_gcl + bf_gi + bf_gmax + set_rescor(TRUE)

# ---- Fit -----------------------------------------------------------------
# Same settings as r/10_fit-all.R, except cores = 3 (matches request_cpus = 3
# in fit_all.sub) so the 3 chains sample in parallel instead of sequentially.

fit <- brm(
  formula = form,
  data = joined_summary |> mutate(phy = accession),
  data2 = list(A = A),
  cores = 3,
  chains = 3,
  iter = thin * 2e3,
  thin = thin,
  refresh = thin * 1e2,
  control = list(adapt_delta = 0.9),
  backend = "cmdstanr",
  family = student(),
  seed = 613135062 + i
) |> add_criterion("loo")

write_rds(fit, out_file)

cli::cli_inform("Wrote {out_file}")
