# testing out mediation calculations
library(mediation)
library(bayestestR)
library(brms)

data(jobs)

# Fit Bayesian mediation model in brms
f1 <- bf(job_seek ~ treat + econ_hard + sex + age)
f2 <- bf(depress2 ~ treat + job_seek + econ_hard + sex + age)
m2 <- brm(f1 + f2 + set_rescor(TRUE), data = jobs, refresh = 0)

mediation(m2)

# direct effect
m2 |> 
  as_draws_df() |>
  dplyr::select(b_depress2_treat) |>
  summarize_draws()

# indirect effect
m2 |> 
  as_draws_df() |>
  dplyr::select(b_jobseek_treat, b_depress2_job_seek) |>
  dplyr::mutate(indirect_effect = b_jobseek_treat * b_depress2_job_seek, .keep = "unused") |>
  summarize_draws()

# mediator effect
m2 |> 
  as_draws_df() |>
  dplyr::select(b_depress2_job_seek) |>
  summarize_draws()

mediation(m2)



source("r/header.R") 
library(bayestestR)

# suppose mfit is your brmsfit object
mfit = read_rds("objects/best_model.rds")
med_summary <- mediation(
  model     = mfit,
  treatment = "lighttreatmentsun",
  mediator  = "logitfgmax",
  response  = c(logitfgmax = "logitfgmax", logtaumean = "logtaumean"),
  centrality = "median",
  ci = 0.95
)

bayestestR:::.mediation
# insight::get_parameters(model)[c(coef_indirect, coef_mediator)]

print(med_summary)
med_summary


post <- as_draws_df(mfit)
colnames(post)[1:50]
post |>
  mutate(cov_lt = cor_phy__logtaumean_Intercept__loggcl_Intercept * sd_phy__logitfgmax_Intercept * sd_phy__logtaumean_Intercept) |>
  pull(cov_lt)
covariance_lt <- post[["cor_phy__logitfgmax_logtaumean"]] *
  post[["sd_phy__logitfgmax_Intercept"]] *
  post[["sd_phy__logtaumean_Intercept"]]
