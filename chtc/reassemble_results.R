# chtc/reassemble_results.R
#
# Run this LOCALLY after downloading all 64 fit_<NN>.rds files from CHTC
# (e.g. via scp/rsync into a local directory) to:
#   1. Verify all 64 files are present and load as class "brmsfit".
#   2. Rebuild objects/df_forms.rds using the EXACT SAME (lambda_idx, tau_idx)
#      mapping used in chtc/fit_model.R, so it matches r/11_compare-models.R's
#      expectations.
#   3. Copy the verified fit_<NN>.rds files into objects/fits/.
#
# Usage:
#   Rscript chtc/reassemble_results.R /path/to/downloaded/fits
#
# If no path is given, defaults to chtc/downloaded_fits/.

source("r/header.R")

args <- commandArgs(trailingOnly = TRUE)
download_dir <- if (length(args) >= 1) args[[1]] else "chtc/downloaded_fits"

if (!dir.exists(download_dir)) {
  cli_abort("Download directory {.path {download_dir}} does not exist.")
}

expected_files <- glue("fit_{str_pad(1:64, 2, 'left', '0')}.rds")
found_files <- list.files(download_dir, pattern = "^fit_[0-9]{2}\\.rds$")

missing <- setdiff(expected_files, found_files)
if (length(missing) > 0) {
  cli_abort(c(
    "Missing {length(missing)} of 64 expected fit files in {.path {download_dir}}:",
    set_names(missing, "*")
  ))
}

# ---- Validate each file loads as a brmsfit -------------------------------

walk(expected_files, function(.f) {
  fit <- read_rds(file.path(download_dir, .f))
  assert_class(fit, "brmsfit")
})

cli_inform("All 64 fit files present and load as brmsfit objects.")

# ---- Rebuild df_forms.rds with the SAME formulas + index mapping used ----
# in chtc/fit_model.R (lambda_idx = (i-1) %/% 8, tau_idx = (i-1) %% 8), so
# that r/11_compare-models.R (and anything else reading objects/df_forms.rds)
# sees the same structure it would from a local r/10_fit-all.R run.

phy = read_rds("data/phylogeny.rds")
A = vcv(phy, corr = TRUE)

bf_lambda0 = bf(
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

bf_lambda1 = update(bf_lambda0, . ~ . - loggcl)
bf_lambda2 = update(bf_lambda0, . ~ . - loggi)
bf_lambda3 = update(bf_lambda0, . ~ . - loggmax)
bf_lambda4 = update(bf_lambda0, . ~ . - loggcl - loggi)
bf_lambda5 = update(bf_lambda0, . ~ . - loggcl - loggmax)
bf_lambda6 = update(bf_lambda0, . ~ . - loggi - loggmax)
bf_lambda7 = update(bf_lambda0, . ~ . - loggcl - loggi - loggmax)

bf_tau0 = bf(
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

bf_tau1 = update(bf_tau0, . ~ . - loggcl)
bf_tau2 = update(bf_tau0, . ~ . - loggi)
bf_tau3 = update(bf_tau0, . ~ . - loggmax)
bf_tau4 = update(bf_tau0, . ~ . - loggcl - loggi)
bf_tau5 = update(bf_tau0, . ~ . - loggcl - loggmax)
bf_tau6 = update(bf_tau0, . ~ . - loggi - loggmax)
bf_tau7 = update(bf_tau0, . ~ . - loggcl - loggi - loggmax)

bf_gcl = bf(loggcl ~
              lighttreatment +
              leaftype +
              (1 | a | accession) +
              (1 | b | gr(phy, cov = A)))

bf_gi = bf(
  loggi ~
    lighttreatment +
    lightintensity +
    leaftype +
    (1 | accid) +
    (1 | a | accession) +
    (1 | b | gr(phy, cov = A))
)

bf_gmax = bf(loggmax ~
               lighttreatment +
               leaftype +
               (1 | a | accession) +
               (1 | b | gr(phy, cov = A)))

bf_lambda_list = list(bf_lambda0, bf_lambda1, bf_lambda2, bf_lambda3,
                       bf_lambda4, bf_lambda5, bf_lambda6, bf_lambda7)
bf_tau_list = list(bf_tau0, bf_tau1, bf_tau2, bf_tau3,
                    bf_tau4, bf_tau5, bf_tau6, bf_tau7)

df_forms = tibble(i = 1:64) |>
  mutate(
    lambda_idx = (i - 1) %/% 8,
    tau_idx = (i - 1) %% 8,
    bf_lambda = bf_lambda_list[lambda_idx + 1],
    bf_tau = bf_tau_list[tau_idx + 1],
    form = map2(
      bf_lambda,
      bf_tau,
      ~ .x + .y + bf_gcl + bf_gi + bf_gmax + set_rescor(TRUE)
    ),
    file = glue("objects/fits/fit_{n}.rds", n = str_pad(i, 2, "left", "0"))
  ) |>
  select(-i, -lambda_idx, -tau_idx)

# ---- Copy files into objects/fits/ and write objects/df_forms.rds --------

dir.create("objects/fits", showWarnings = FALSE, recursive = TRUE)
walk(expected_files, function(.f) {
  file.copy(
    file.path(download_dir, .f),
    file.path("objects/fits", .f),
    overwrite = TRUE
  )
})

write_rds(df_forms, "objects/df_forms.rds")

cli_inform("Copied 64 fit files into objects/fits/ and wrote objects/df_forms.rds")
