# WORK IN PROGRESS
# Analyze effect of anatomy on tau
source("r/header.R")

joined_data = read_rds("data/joined-data.rds") |>
  # need to do more careful outlier removal
  filter(log_tau_mean < 7) |>
  mutate(
    log_gcl = log(guard_cell_length_um),
    log_invpi = log(1 / f_gmax)
  )

n = 1
fit = brm(
  log_tau_mean |
    se(log_tau_sd, sigma = TRUE) ~ light_intensity + light_treatment + log_gcl + log_invpi + light_intensity:log_gcl + light_intensity:log_invpi + light_treatment:log_gcl + light_treatment:log_invpi + (1 |
                                                                                                                                                                                                            acc_id) + (
                                                                                                                                                                                                              light_intensity + light_treatment + log_gcl + log_invpi |
                                                                                                                                                                                                                accession
                                                                                                                                                                                                            ),
  data = filter(joined_data, curve_type == "amphi"),
  cores = 1,
  chains = 1,
  iter = n * 2e3,
  thin = n,
  refresh = n * 1e2,
  backend = "cmdstanr",
  seed = 12345,
  control = list(max_treedepth = 10)
)

conditional_effects(fit)
bayes_R2(fit)

# checking to make sure I understand scaling correctly
n = 1e5
df1 = tibble(
  log_x1 = rnorm(n, 0, 0.1),
  log_x2 = rnorm(n, 0, 0.1),
  x1 = exp(log_x1),
  x2 = exp(log_x2),
  y = x1 / (4 * (3 - x2)),
  log_y = log(y) + rnorm(n, 0, 0.1)
)
ggplot(df1, aes(log_x2, log_y)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, color = "red") +
  theme_minimal()
lm(log_y ~ log_x1 + log_x2, data = df1) |> summary()

# fitting many models and comparing. not sure if I want to do this
set.seed(443448228)
df1 = tibble(form = all_hierarchical_models(
  c("light_intensity", "light_treatment", "log_gcl", "log_invpi"),
  response = "log_tau_mean | se(log_tau_sd, sigma = TRUE)"
)) |>
  mutate(seed = sample(1e9, n(), replace = FALSE))

plan(multisession, workers = 10)
tmp = future_map2(df1$form, df1$seed, \(.f, .s) {
  fit = brm(
    as.formula(paste0(.f, " + (1 | acc_id) + (1 | accession)")),
    data = filter(joined_data, curve_type == "amphi"),
    cores = 1,
    chains = 1,
    backend = "cmdstanr",
    seed = .s
  ) |>
    add_criterion(criterion = "loo")
  fit
}, .progress = TRUE, .options = furrr_options(seed = TRUE))

loos <- lapply(tmp, function(x) x$criteria$loo)
names(loos) <- str_c("Model_", seq_along(loos))
loo_compare(loos)

tmp[[5]]
tmp[[163]]


ggplot(joined_data, aes(lower_guard_cell_length_um, log_tau_mean, color = light_treatment)) +
  facet_wrap(~light_intensity) +
  geom_point() +
  scale_x_log10() 

ggplot(joined_data, aes(log_invpi, log_tau_mean, color = light_treatment)) +
  facet_wrap(~light_intensity) +
  geom_point() 

c(
  "1",
  "light_intensity",
  "light_treatment",
  "log_gcl",
  "log_invpi",
  "light_intensity + light_treatment",
  "light_intensity + log_gcl",
  "light_intensity + log_invpi",
  "light_treatment + log_gcl",
  "light_treatment + log_invpi",
  "log_gcl + log_invpi",
  "light_intensity + light_treatment + log_gcl",
  "light_intensity + light_treatment + log_invpi",
  "light_intensity + log_gcl + log_invpi",
  "light_treatment + log_gcl + log_invpi",
  "light_intensity + light_treatment + log_gcl + log_invpi"
) |> length()

n = 1
fit0 = brm(
  log_tau_mean | se(log_tau_sd, sigma = TRUE) ~ 1 + (1 | acc_id) + (1 | accession),
  data = filter(joined_data, curve_type == "amphi"),
  cores = 1,
  chains = 1,
  iter = n * 2e3,
  thin = n,
  refresh = n * 1e2,
  backend = "cmdstanr",
  seed = 12345,
  control = list(max_treedepth = 10)
)

# Test for main effect of measurement light intensity
fit1 = brm(
  log_tau_mean | se(log_tau_sd) ~ light_intensity + (1 | acc_id) + (1 | accession),
  data = filter(joined_data, curve_type == "amphi"),
  cores = 4,
  chains = 4,
  backend = "cmdstanr"
)

# Test for random effect of measurement light intensity
fit2 = brm(
  log_tau_mean | se(log_tau_sd) ~ light_intensity + (1 | acc_id) + (1 + light_intensity | accession),
  data = filter(joined_data, curve_type == "amphi"),
  cores = 4,
  chains = 4,
  backend = "cmdstanr"
)


# test code for differing within versus between group slopes
b0_between = 0
b1_between = 1
b1_within = -1
n_group = 1e1
n_rep = 1e1

df1 = crossing(nesting(
  group = factor(1:n_group),
  X = seq(-1, 1, length.out = n_group),
  Y = b0_between + b1_between * X + rnorm(n_group, 0, 1)
),
rep = 1:n_rep) |>
  mutate(
    x = X + rnorm(n(), 0, 1),
    y = Y + b1_within * (x - X) + rnorm(n(), 0, 1)
  )

ggplot(filter(df1, rep == 1), aes(X, Y)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, color = "red") +
  theme_minimal()

ggplot(df1, aes(x, y)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, color = "red") +
  theme_minimal()

ggplot(filter(df1, group == 2), aes(x, y)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, color = "red") +
  theme_minimal()

fit = brm(y ~ x + (1 + x | group), data = df1, cores = 1, chains = 1, backend = "cmdstanr")
summary(fit)
ranef(fit)
