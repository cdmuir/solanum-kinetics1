# Analyze effect of anatomy on tau in untreated (amphistomatous) leaves
# Note: I did some preliminary model comparisons suggesting no random effects
# of light_intensity, light_treatment, log_gcl, or log_fgmax
source("r/header.R")

fit_dir = "objects/fits"
if (!dir.exists(fit_dir)) {
  dir.create(fit_dir)
}

plan(multisession, workers = 19)

joined_summary = read_rds("data/joined-summary.rds") |>
  prepare_tau_anatomy_data(logtau_threshold)

write_rds(
  list(
    logtau_threshold = logtau_threshold,
    n_removed = attr(joined_summary, "n_removed")
  ),
  "objects/n_removed.rds"
)

assert_true(all(!is.na(joined_summary$loglambdamean)))
assert_true(all(!is.na(joined_summary$loglambdasd)))
assert_true(all(!is.na(joined_summary$logtausd)))
assert_true(all(!is.na(joined_summary$logtausd)))
assert_true(all(!is.na(joined_summary$loggcl)))
assert_true(all(!is.na(joined_summary$logitfgmax)))

assert_true(all(!is.na(joined_summary$lighttreatment)))
assert_true(all(!is.na(joined_summary$lightintensity)))
assert_true(all(!is.na(joined_summary$leaftype)))

phy = read_rds("data/phylogeny.rds")
A = vcv(phy, corr = TRUE)
thin = 1 # 6

# Define formulae
expl_var_combn = c(
  "",
  " + {.var}",
  " + {.var} + {.var}:lighttreatment",
  " + {.var} + {.var}:lightintensity",
  " + {.var} + {.var}:leaftype",
  " + {.var} + {.var}:lighttreatment + {.var}:lightintensity",
  " + {.var} + {.var}:lighttreatment + {.var}:leaftype",
  " + {.var} + {.var}:lightintensity + {.var}:leaftype",
  " + {.var} + {.var}:lighttreatment + {.var}:lightintensity + {.var}:leaftype"
)

template_form = "{.resp} ~ {.base}{.gcl}{.fgmax}{.rand}"

lhs = crossing(
  .base = "lighttreatment + lightintensity + leaftype",
  .gcl = map_chr(expl_var_combn, glue, .var = "loggcl"),
  .fgmax = map_chr(expl_var_combn, glue, .var = "logitfgmax"),
  .rand = " + (1|accid) + (1|a|accession) + (1|b|gr(phy, cov = A))"
)

n0 = floor(log10(nrow(lhs)^2)) + 1

bf_gcl = bf(loggcl ~
              lighttreatment +
              leaftype +
              (1 | a | accession) +
              (1 | b | gr(phy, cov = A)))

bf_fgmax = bf(
  logitfgmax ~
    lighttreatment +
    lightintensity +
    leaftype +
    (1 | accid) +
    (1 | a | accession) +
    (1 | b | gr(phy, cov = A))
)

crossing(
  lhs |>
    mutate(.resp = "loglambdamean | se(loglambdasd, sigma = TRUE)") |>
    mutate(form_lambda = glue(template_form), .keep = "unused"),
  lhs |>
    mutate(.resp = "logtaumean | se(logtausd, sigma = TRUE)") |>
    mutate(form_tau = glue(template_form), .keep = "unused")
) |>
  mutate(i = row_number(), fit_dir = fit_dir) |>
  # randomize order
  slice_sample(n = nrow(lhs) ^ 2) |>
  future_pwalk(\(form_lambda, form_tau, i, fit_dir) {
    file_name = glue("{fit_dir}/fit{.n}.rds",
                     .n = str_pad(
                       i,
                       width = n0,
                       pad = "0",
                       side = "left"
                     ))
    if (!file.exists(file_name)) {
      bf_lambda = bf(form_lambda)
      bf_tau = bf(form_tau)
      brm(
        formula = bf_lambda +
          bf_tau +
          bf_gcl +
          bf_fgmax +
          set_rescor(TRUE),
        data = joined_summary |>
          mutate(phy = accession),
        data2 = list(A = A),
        cores = 1,
        chains = 1,
        iter = thin * 2e3,
        thin = thin,
        refresh = 0, #thin * 1e2,
        # control = list(adapt_delta = 0.9),
        backend = "cmdstanr",
        family = student(),
        seed = 613135062 + i,
        file = file_name
      )
    } else {
      message("File ", file_name, " already exists. Skipping fit.")
    }
  }, .options = furrr_options(seed = TRUE), .progress = TRUE)

zipr(str_c(fit_dir, ".zip"),
     list.files(fit_dir, full.names = TRUE),
     recurse = TRUE)
