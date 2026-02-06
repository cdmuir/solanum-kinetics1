# Fit CDWeibulll mode to each curve separately using nls
source("r/header.R")

rh_curves = read_rds("data/rh_curves.rds")

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
