# Fit CDWeibulll mode to each curve separately
source("r/header.R")

sk_dir = "objects/sk-curves/"
if (!dir.exists(sk_dir)) {
  dir.create(sk_dir)
}

rh_curves = read_rds("https://github.com/cdmuir/solanum-aa/raw/refs/heads/main/data/trimmed_rh_curves.rds") |>
  mutate(ci = as.numeric(as.factor(curve))) |>
  mutate(t_sec = elapsed - min(elapsed), .by = ci)

form_vico = gsw ~ gf + dg * exp(-(t_sec / tau))  
form_cdweibull = gsw ~ gf + dg * exp(-(t_sec / tau) ^ lambda)

bform_cdweibull = bf(
  form_cdweibull,
  gf ~ 1,
  dg ~ 1,
  tau ~ 1,
  lambda ~ 1,
  nl = TRUE
)

# nls ----
safe_nls = safely(nls)

plan(multisession, workers = 10)
nls_summary = rh_curves |>
  split(~ curve) |>
  future_imap_dfr(\(df, curve_id) {

    init = list(gf = min(df$gsw),
                dg = diff(range(df$gsw)),
                tau = 50)
    
    fit_vico = safe_nls(form_vico, data = df, start = init)
    
    if (is.null(fit_vico$result)) {
      b_vico = tibble(
        parameter = names(init),
        value = NA,
        aic = NA,
        model = "vico"
      )
      init = c(init, lambda = 1)
    } else {
      b_vico = as_tibble(coef(fit_vico$result), rownames = "parameter") |>
        mutate(model = "vico", aic = AIC(fit_vico$result))
      init = c(as.list(coef(fit_vico$result)), lambda = 1)
    }
    
    fit_cdweibull = safe_nls(form_cdweibull, data = df, start = init)
    
    if (is.null(fit_cdweibull$result)) {
      b_cdweibull = tibble(
        parameter = names(init),
        value = NA,
        aic = NA,
        model = "cdweibull"
      )
    } else {
      b_cdweibull = as_tibble(coef(fit_cdweibull$result), rownames = "parameter") |>
        mutate(model = "cdweibull", aic = AIC(fit_cdweibull$result))
    }
    
    bind_rows(b_vico, b_cdweibull) |>
      mutate(curve = curve_id)
    
  }, .progress = TRUE, .options = furrr_options(seed = TRUE))

write_rds(nls_summary, "objects/nls_summary.rds")

# brms ----
rh_curves |>
  split( ~ curve) |>
  future_iwalk(\(df, curve_id) {
    file = paste0(sk_dir, curve_id, ".rds")
    if (!file.exists(file)) {
      
      # need to figure out somethign with priors/init
      # prior1 = prior(normal(1/3, 1/3), nlpar = "gf") +
      #   prior(normal(-1/4, 1/4), nlpar = "dg", ub = 0) +
      #   # prior(normal(-1, 1), nlpar = "lambda", ub = 0) +
      #   prior(normal(100, 100), nlpar = "tau", lb = 0)
      # df = filter(rh_curves, curve == "LA0107-C_pseudohypo_2000")
      # init = list(gf = min(df$gsw),
      #             dg = diff(range(df$gsw)),
      #             tau = 50)
      # fit_vico = nls(form_vico, data = df, start = init)
      # init = c(as.list(coef(fit_vico)), lambda = 1)
      # fit_cdweibull = nls(form_cdweibull, data = df, start = init)

      x = 2
      ad = 0.8
      n_divergent = Inf
      
      while (n_divergent > 10 & x < 10) {
        fit_curve = brm(
          formula = bform_cdweibull,
          iter = x * 2000,
          thin = x,
          data = df,
          chains = 1,
          cores = 1,
          # prior = prior1,
          backend = "cmdstanr",
          control = list(adapt_delta = ad),
          seed = 360036340 + df$ci[1]
        )
        n_divergent = nuts_params(fit_curve) |>
          subset(Parameter == "divergent__") |>
          pull(Value) |>
          sum()
        x = x + 1
        ad = min(ad * 1.1, 0.99)
      }
      
      write_rds(fit_curve, file)
    } else {
      message(paste0("Curve fit for ", curve_id, " already exists. Skipping."))
    }
  }, .progress = TRUE, .options = furrr_options(seed = TRUE))

# zip
zip::zipr(
  "objects/sk-curves.zip",
  list.files(sk_dir, full.names = TRUE),
  recurse = TRUE
)
