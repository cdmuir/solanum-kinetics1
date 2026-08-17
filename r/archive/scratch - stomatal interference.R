# This basically shows that stomatal interference alone (using the Lehmann and Or 2015 equations) cannot explain variation in tau. It doesn't generate enough variance and only predicts very weak association between 
source("r/header.R")

k <- 1e-3 # D_wv / v (mol / s / m)
n <- 200 * 1e6 # stomatal density (m^-2)
r <- 10 / 1e6 # pore radius (m)
d <- 10 / 1e6  # pore depth (m)
# Eqn 10 from Lehman and Or (2015)
g_leaf = (n * k * 2 * pi * r ^ 2) / (2 * d + r * (pi - 2 * r * sqrt(n)))

imply_r = function(
    g_leaf,              # [mol / m^2 / s]
    stomatal_density_m2, # [1 / m^2]
    pore_depth,          # [m]
    D_wv = 2.49e-5,      # [m^2 / s]
    v = 2.24e-2          # [m^3 / mol]
) {
  
  # Rename variables to match Lehmann and Or (2015)
  k = D_wv / v
  n = stomatal_density_m2
  d = pore_depth

  # Solve for r
  A = 2 * pi * n * k + 2 * g_leaf * sqrt(n)
  B = -g_leaf * pi
  C = -2 * g_leaf * d

  r = (-B + sqrt(B^2 - 4 * A * C)) / (2 * A) # [m]
  
  r
  
}

r_to_gleaf = function(
    r,                   # [m]
    stomatal_density_m2, # [1 / m^2]
    pore_depth,          # [m]
    D_wv = 2.49e-5,      # [m^2 / s]
    v = 2.24e-2          # [m^3 / mol]
) {
  # Rename variables to match Lehmann and Or (2015)
  k = D_wv / v
  n = stomatal_density_m2
  d = pore_depth
  g_leaf = (n * k * 2 * pi * r ^ 2) / (2 * d + r * (pi - 2 * r * sqrt(n)))
  g_leaf
}

# assume pore depth = 0.25 * guard cell length
joined_data <- read_rds("data/joined-data.rds")
joined_summary <- read_rds("data/joined-summary.rds") 

sim_pars <- joined_summary |>
  mutate(
    pore_depth = 0.25 * guard_cell_length_um / 1e6,
    stomatal_density_m2 = 1e6 * stomatal_density_mm2,
    r_init = imply_r(ginit_mean, stomatal_density_m2, pore_depth),
    r_final = imply_r(gfinal_mean, stomatal_density_m2, pore_depth)
  ) |>
  full_join(joined_data |> select(id = curve, t_sec, gsw), by = join_by(id))

global_logtau = mean(joined_summary$logtau_mean)

sim_one_lehmannor_curve <- function(p, k) {
  p |>
    mutate(
      r = r_final + (r_init - r_final) * exp(-k * t_sec),
      gsw_sim = r_to_gleaf(r, stomatal_density_m2, pore_depth)
    ) 
}

id_val = unique(sim_pars$id)[1]
p = sim_pars |>
  filter(id == id_val)

fit = optim(
  first(p$logtau_mean),
  fn = function(.logtau, .p) {
    k = 1 / exp(.logtau)
    sum((.p$gsw - sim_one_lehmannor_curve(.p, k)$gsw_sim) ^ 2)
  }, .p = p)

p1 = sim_one_lehmannor_curve(p, 1 / exp(fit$par)) 
fit_nls_one(p1)

tmp = sim_one_lehmannor_curve(sim_pars, 1 / exp(global_logtau))
tmp1 = tmp |>
  split(~ id) |>
  # magrittr::extract(1:3) |>
  map(fit_nls_one) |>
  bind_rows(.id = "id")

joined_summary |>
  filter(logtau_mean < logtau_threshold) |>
  full_join(tmp1, by = join_by(id)) |>
  ggplot(aes(exp(logtau_mean), tau_hat)) +
  geom_point() +
  scale_x_log10() +
  scale_y_log10() +
  coord_equal()

joined_summary |>
  filter(logtau_mean < logtau_threshold) |>
  full_join(tmp1, by = join_by(id)) |>
  ggplot(aes(exp(loglambda_mean), lambda_hat)) +
  geom_point() +
  scale_x_log10() +
  scale_y_log10() +
  coord_equal()

joined_summary |>
  filter(logtau_mean < logtau_threshold) |>
  full_join(tmp1, by = join_by(id)) |>
  ggplot(aes(ginit_mean, tau_hat)) +
  geom_point() +
  scale_x_log10() +
  scale_y_log10()

# some stuff using units to make sure I understand the lehmann and or equations
library(units)
# Stomatal density (n per mm^2)
n = set_units(10, mm^-2)

# guard cell length (um)
gcl = set_units(30, um)

# r = pore radius at maximum using Sack and Buckley (2016) assumptions
# rmax = 0.5 * p
# p = c * gcl where c = 0.5
r = gcl / 4

# interstomatal distance (Muir 2020)
s = (2 / sqrt(3) / n ) ^ 0.5
s

# k
# biophysical_constant <- function(D_wv, v) set_units(D_wv / v, mol / m / s)
D_wv = set_units(2.49e-5, m^2 / s)
v = set_units(2.24e-2, m^3 / mol)
k = D_wv / v

# Eqn 2s - R_end
R_end = 1 / (4 * r * k * n)

set_units(1 / R_end, mol / m^2 / s)

# Eqn 2c - R_VS
R_VS = (1 / (4 * r) - 1 / (pi * s)) / (k * n)

R_VS / R_end

set_units(1 / R_VS, mol / m^2 / s)
