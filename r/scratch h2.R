fit_amphi_high = read_rds("objects/fit_amphi_high.rds")

lm(loggcl ~ lighttreatment + phy, data = fit_amphi_high$data) |> 
  aov() |>
  summary()

# phylo model
VarCorr(fit_amphi_high)$phy$sd["loggcl_Intercept",] # 0.15934622
VarCorr(fit_amphi_high)$residual$sd["loggcl",] # 0.112379640


# nonphylo model
VarCorr(fit_amphi_high2)$accession$sd["loggcl_Intercept",] # 0.09569342
VarCorr(fit_amphi_high2)$residual$sd["loggcl",] # 0.11238601

hyp_gcl <- paste(
  "sd_phy__loggcl_Intercept^2 /",
  "(sd_phy__loggcl_Intercept^2 + sd_accession__loggcl_Intercept^2 + sigma_loggcl^2) = 0"
)
hypothesis(fit_amphi_high0, hyp_gcl, class = NULL)

hyp_fgmax <- paste(
  "sd_phy__logfgmax_Intercept^2 /",
  "(sd_phy__logfgmax_Intercept^2 + sd_accession__logfgmax_Intercept^2 + sigma_logfgmax^2) = 0"
)
hypothesis(fit_amphi_high0, hyp_fgmax, class = NULL)

hyp <- paste(
  "sd_phy__logfgmax_Intercept^2 /",
  "(sd_phy__logfgmax_Intercept^2 + sd_accession__logfgmax_Intercept^2 + sigma_logfgmax^2) = sd_phy__loggcl_Intercept^2 /",
  "(sd_phy__loggcl_Intercept^2 + sd_accession__loggcl_Intercept^2 + sigma_loggcl^2)"
)

thin = 1
bf1 = bf(logtaumean | se(logtausd, sigma = TRUE) ~ lighttreatment + (-1 + lighttreatment|a|accession) + (-1 + lighttreatment|b|gr(phy, cov = A)))
bf2 = bf(loglambdamean | se(loglambdasd, sigma = TRUE) ~ lighttreatment + (-1 + lighttreatment|a|accession) + (-1 + lighttreatment|b|gr(phy, cov = A)))
bf3 = bf(loggcl ~ lighttreatment + (-1 + lighttreatment|a|accession) + (-1 + lighttreatment|b|gr(phy, cov = A)))
bf4 = bf(logfgmax ~ lighttreatment + (-1 + lighttreatment|a|accession) + (-1 + lighttreatment|b|gr(phy, cov = A)))


fit_amphi_highX = brm(
  bf1 + bf2 + bf3 + bf4 + set_rescor(TRUE),
  data = joined_summary |>
    filter(curve_type == "amphi", lightintensity == "high") |>
    mutate(phy = accession),
  data2 = list(A = A),
  cores = 1,
  chains = 1,
  iter = thin * 2e3,
  thin = thin,
  refresh = thin * 1e2,
  backend = "cmdstanr",
  seed = 613135062,
) |>
  add_criterion("loo")

hyp_gcl <- paste(
  "sd_phy__loggcl_lighttreatmentsun^2 /",
  "(sd_phy__loggcl_lighttreatmentsun^2 + sd_accession__loggcl_lighttreatmentsun^2 + sigma_loggcl^2) = 0"
)
hypothesis(fit_amphi_highX, hyp_gcl, class = NULL)

hyp_gcl <- paste(
  "sd_phy__loggcl_lighttreatmentshade^2 /",
  "(sd_phy__loggcl_lighttreatmentshade^2 + sd_accession__loggcl_lighttreatmentshade^2 + sigma_loggcl^2) = 0"
)
hypothesis(fit_amphi_highX, hyp_gcl, class = NULL)

hyp_fgmax <- paste(
  "sd_phy__logfgmax_lighttreatmentsun^2 /",
  "(sd_phy__logfgmax_lighttreatmentsun^2 + sd_accession__logfgmax_lighttreatmentsun^2 + sigma_logfgmax^2) = 0"
)
hypothesis(fit_amphi_highX, hyp_fgmax, class = NULL)

hyp_fgmax <- paste(
  "sd_phy__logfgmax_lighttreatmentshade^2 /",
  "(sd_phy__logfgmax_lighttreatmentshade^2 + sd_accession__logfgmax_lighttreatmentshade^2 + sigma_logfgmax^2) = 0"
)
hypothesis(fit_amphi_highX, hyp_fgmax, class = NULL)

hyp <- "sd_phy__logfgmax_lighttreatmentshade^2 / (sd_phy__logfgmax_lighttreatmentshade^2 + sd_accession__logfgmax_lighttreatmentshade^2 + sigma_logfgmax^2) = sd_phy__loggcl_lighttreatmentshade^2 / (sd_phy__loggcl_lighttreatmentshade^2 + sd_accession__loggcl_lighttreatmentshade^2 + sigma_loggcl^2)"

hypothesis(fit_amphi_highX, hyp, class = NULL)
