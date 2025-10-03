library(dplyr)
library(ggplot2)
library(purrr)
library(tibble)
library(tidyr)

rm(list = ls())

pars = list(
  # Main model
  chi = 0.105, # turgor-to-conductance scaling factor
  beta = 1.17, # Hydromechanical/biochemical response parameter
  D_s = 10, # Leaf-to-boundary layer H20 mole fraction gradient
  psi_s = 0, # Source water potential
  pi_e = 0.525, # Epidermal osmotic pressure
  pi_a = 0, # Apoplastic osmotic pressure
  r_bw = 1 / 3.3, # Boundary layer resistance to water vapor
  R = 0.0456, # Effective hydraulic resistance to the epidermis
  R_d = 2, # Not specified in Buckley, but this is what I used in photosynthesis package
  M = 0.98, # REsidual epidermal mechanical advantage
  rho = 0, # guard cell resistive advantage
  
  # ATP submodel
  K_c = 40.4, # Michaelis constant for RuBP carboxylation
  K_o = 2.48e4, # Michaelis constant for RuBP oxygenation
  J_m = 2.02 * 8.86e1, # Light-saturated potential electron transport rate (2.02 * V_m)
  theta_J = 0.908, # Curvature of the light response of J
  I = 1100, # Light intensity (umol m-2 s-1)
  f = 0.195, # Product of absorbance and effective quantum yield
  pO2 = 2.1e4, # Partial pressure of O2 (Pa)
  c_a = 365, # Ambient CO2 concentration (umol mol-1),
  p_t = 1e5, # Atmospheric pressure (Pa)
  tau_0 = 11.6, # Basal ATP level provided by other processes
  a_t = 12.6 * 8.86e1 / 1e3, # Total concentration of adenylates (tau + [ADP])
  p = 2.5 * 8.86e1 / 1e3, # Concentration of photophosphorylation sites
  V_m = 8.86e1, # Carboxylation rate Limited by Rubisco only
  V_r = 2.27 * 8.86e1 # Carboxylation rate limited by potential RuBP pool size only
)

init = list(
  p_i = pars$c_a / 2
)

get_Wc = function(p_i, pars) {
  # Eqn A19
  # 10x is to convert from Pa to ppm assuming total pressure is 1e5 Pa
  with(pars, (V_m * p_i) / (p_i + 10 * K_c * (1 + pO2 / K_o)))
}

get_J = function(pars) {
  # Eqn A21
  a = pars$theta_J
  b = -(pars$J_m + pars$f * pars$I)
  c = pars$f * pars$I * pars$J_m
  min((-b + c(-1, 1) * sqrt(b^2 - 4 * a * c)) / (2 * a))
}

get_gammastar = function(pars) {
  with(pars, 0.105 * K_c * pO2 / K_o)
}

get_Wj = function(p_i, pars) {
  # Eqn A20
  J = get_J(pars)
  # 10x is to convert from Pa to ppm assuming total pressure is 1e5 Pa
  gamma_star = 10 * get_gammastar(pars)
  (J * p_i) / (4 * (p_i + 2 * gamma_star))
}


get_A = function(p_i, pars) {
  # 10x is to convert from Pa to ppm assuming total pressure is 1e5 Pa
  gamma_star = 10 * get_gammastar(pars)
  W = pmin(get_Wc(p_i, pars), get_Wj(p_i, pars))
  (1 - gamma_star / p_i) * W - pars$R_d
}

get_g2 = function(p_i, pars) {
  # Eqn A25
  A = get_A(p_i, pars)
  A1 = with(pars, A / (c_a - p_i / p_t))
  # Derivation
  # g = A1 * (0.23 + 1.37 * omega)
  # g = A1 * (0.23 + 1.37 / (1 + g * pars$r_bw))
  # g = A1 * 0.23 + A1 * 1.37 / (1 + g * pars$r_bw)
  # 0 = g - A1 * 0.23 - A1 * 1.37 / (1 + g * pars$r_bw)
  # 0 = g * (1 + g * pars$r_bw) - A1 * 0.23 * (1 + g * pars$r_bw) - A1 * 1.37
  # 0 = g + g ^ 2 * pars$r_bw - 
  #   A1 * 0.23 - A1 * 0.23 * g * pars$r_bw - A1 * 1.37
  # 0 = g ^ 2 * pars$r_bw + g * (1 - A1 * 0.23 * pars$r_bw) - 
  #   A1 * 1.6
  a = pars$r_bw
  b = 1 - A1 * 0.23 * pars$r_bw
  c = -A1 * 1.6
  pmin((-b + sqrt(b^2 - 4 * a * c)) / (2 * a), (-b + sqrt(b^2 - 4 * a * c)) / (2 * a))
}

get_tau = function(p_i, pars) {
  # Eqn A24
  W_c = get_Wc(p_i, pars)
  W_j = get_Wj(p_i, pars)
  tau_c = with(pars, a_t - p * W_c / W_j)
  tau_j = with(pars, (a_t - p) * (V_r / V_m - 1) / (W_c / W_j * V_r / V_m - 1))
  pars$tau_0 + tau_c * (W_c < W_j) + tau_j * (W_c >= W_j)
}

get_g1 = function(p_i, pars) {
  tau = get_tau(p_i, pars)
  with(pars, chi * ((beta * tau - M) * (psi_s + pi_e) - pi_e + pi_a) / (1 + chi * R * D_s * (beta * tau - M + rho)))
}

pi_obj = function(p_i, pars) {
  # Objective function for optimization
  g1 = get_g1(p_i, pars)
  g2 = get_g2(p_i, pars)
  (g1 - g2)^2
}

optimize_p_i = function(pars) {
  tibble(
    p_i = seq_len(pars$c_a),
    g1 = get_g1(p_i, pars),
    g2 = get_g2(p_i, pars),
    d = (g1 - g2)^2
  ) |>
    filter(d == min(d)) |>
    pull(p_i) |>
    optim(
      pi_obj,
      pars = pars,
      method = "Brent",
      lower = 0,
      upper = pars$c_a,
      control = list(factr = 1e16)
    )
}

df1 = map_dfr(seq(40, 400, 1), \(.x, pars0) {
  
  pars0$I = 1000
  tibble(
    p_i = .x,
    W_c = get_Wc(.x, pars0),
    W_j = get_Wj(.x, pars0),
    tau = get_tau(.x, pars0)
  )
}, pars0 = pars)

# I am not sure why I made this figure. ignore?
df1 |>
  pivot_longer(c(W_c, W_j), names_to = "type", values_to = "W") |>
  ggplot(aes(p_i, W, color = type)) +
  geom_line() +
  labs(
    x = "Intercellular CO2 concentration (umol mol-1)",
    y = "Water flux (mol m-2 s-1)"
  ) +
  theme_minimal()

ggplot(df1, aes(p_i, tau)) +
  geom_line() +
  theme_minimal()


get_g = function(pars) {
  p_i = optimize_p_i(pars)$par
  get_g1(p_i, pars)
}

df1 = map_dfr(seq(10, 30, 1), \(.x, pars0) {
  pars0$D_s = .x
  tibble(
    D_s = .x,
    g = get_g(pars0)
  )
}, pars0 = pars) |>
  mutate(rel_g= g / max(g))

# Figure 3, lower-left panel
ggplot(df1, aes(D_s, rel_g)) +
  geom_line() +
  labs(
    x = "Leaf-to-boundary layer H20 mole fraction gradient",
    y = "Stomatal conductance (mol m-2 s-1)"
  ) +
  theme_minimal()

