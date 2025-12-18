library(cowplot)
library(dplyr)
library(ggplot2)
library(tibble)

theme_set(theme_cowplot())

# dP / dt = (P* - P) / tau_P
# P(t) = P* - (P* - P0) * exp(-t / tau_P)
# 
# da / dP = b0 * exp(-b1 * P)
# Assume a(P = 0) = a0
# a(P) = a0 + (b0 / b1) * (1 - exp(-b1 * P))

a1 = function(P, a0, b0, b1) {
  a0 + (b0 / b1) * (1 - exp(-b1 * P))
}

# Check function
tibble(P = seq(0, 6, 0.01), a1 = a1(P, 150, 800, 1)) |>
  ggplot(aes(x = P, y = a1)) +
  geom_line(size = 1, color = "blue") +
  theme(
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14)
  ) +
  ylim(105, 950)

# By chain rule:
# da / dP = (da / dt) / (dP / dt) = b0 * exp(-b1 * P)
# da / dt = (dP / dt) * b0 * exp(-b1 * P)
# 
# Subbing in dP / dt:
# da / dt = (P* - P) / tau_P * b0 * exp(-b1 * P)
# 
# Rearrnage equation above to get P in terms of a (assume a < a_max, b1 > 0):
# P(a) = log(b0/(-a*b1 + a0*b1 + b0))/b1

P1 = function(a, a0, b0, b1) {
  log(b0/(-a*b1 + a0*b1 + b0))/b1
}

# Check function
a = a1(P = 3, a0 = 150, b0 = 800, b1 = 1)
P1(a = a, a0 = 150, b0 = 800, b1 = 1)

# Use this equation to put da / dt in terms of a:
# P* - P = log(b0/(-a_star*b1 + a0*b1 + b0))/b1 - log(b0/(-a*b1 + a0*b1 + b0))/b1

# Simplify
P_star <- 1
a0 = 150
b0 = 800
b1 = runif(1)
a_star = a1(P_star, a0, b0, b1)
a = a_star * 1.1
log(b0/(-a_star*b1 + a0*b1 + b0))/b1 - log(b0/(-a*b1 + a0*b1 + b0))/b1
(log(b0/(-a_star*b1 + a0*b1 + b0)) - log(b0/(-a*b1 + a0*b1 + b0)))/b1
(log((b0/(-a_star*b1 + a0*b1 + b0)) / (b0/(-a*b1 + a0*b1 + b0))))/b1
(log((1/(-a_star*b1 + a0*b1 + b0)) / (1/(-a*b1 + a0*b1 + b0))))/b1
(log((-a*b1 + a0*b1 + b0) / (-a_star*b1 + a0*b1 + b0)))/b1
log((b0 + b1 * (a0 - a)) / (b0 + b1 * (a0 - a_star)))/b1
# Therefore:
# P* - P = log((b0 + b1 * (a0 - a)) / (b0 + b1 * (a0 - a_star)))/b1

# b0 * exp(-b1 * P) = b0 * exp(-b1 * (log(b0/(-a*b1 + a0*b1 + b0))/b1))

b0 * exp(-b1 * (log(b0/(-a*b1 + a0*b1 + b0))/b1))
b0 * exp(-log(b0/(-a*b1 + a0*b1 + b0)))
b0 * (-a*b1 + a0*b1 + b0) / b0
-a*b1 + a0*b1 + b0
b0 + b1 * (a0 - a)

# Therefore:
# b0 * exp(-b1 * P) = b0 + b1 * (a0 - a)

# Putting these together:
# da / dt = log((b0 + b1 * (a0 - a)) / (b0 + b1 * (a0 - a_star)))/ (b1 * tau_P) * (b0 + b1 * (a0 - a))

# Check
P <- 3
P_star <- 1
a0 = 150
b0 = 800
b1 = runif(1)
a_star = a1(P_star, a0, b0, b1)
a = a1(P = P, a0 = a0, b0 = b0, b1 = b1)
tau_P = 600

(P_star - P) / tau_P * b0 * exp(-b1 * P)
log((b0 + b1 * (a0 - a)) / (b0 + b1 * (a0 - a_star)))/ (b1 * tau_P) * (b0 + b1 * (a0 - a))

# From python

at = function(t, a_init, a0, b0, b1, a_star, tau_P) {
  # (-b0*exp(exp(-t/tau_P)*log((a0*b1 - a_init*b1 + b0)/(a0*b1 - a_star*b1 + b0))) + b0 + b1*(-a0*exp(exp(-t/tau_P)*log((a0*b1 - a_init*b1 + b0)/(a0*b1 - a_star*b1 + b0))) + a0 + a_star*exp(exp(-t/tau_P)*log((a0*b1 - a_init*b1 + b0)/(a0*b1 - a_star*b1 + b0)))))/b1
  A = log((a0*b1 - a_init*b1 + b0)/(a0*b1 - a_star*b1 + b0))
  B = log((a0*b1 - a_init*b1 + b0)/(a0*b1 - a_star*b1 + b0))
  (-b0*exp(exp(-t/tau_P)*A) + b0 + (-a0*exp(exp(-t/tau_P)*A) + a0 + a_star*exp(exp(-t/tau_P)*B)))
  
}

# simplified
at1 = function(t, a_init, a0, b0, b1, a_star, tau_P) {
  
  # option 1
  R = (a0*b1 - a_init*b1 + b0) / (a0*b1 - a_star*b1 + b0)
  # a_exact = a0 + b0/b1 + (a_star - a0 - b0/b1) * R**exp(-t/tau_P)

  # option 2
  # note that a_max = a0 + b0 / b1
  a_max = a0 + b0 / b1
  R = (a_max - a_init) / (a_max - a_star)
  a_exact = a_max + (a_star - a_max) * R ^ exp(-t/tau_P)
  
  # option 3
  # a_exact = (a0*b1 + b0 - ((a0*b1 - a_init*b1 + b0)/(a0*b1 - a_star*b1 + b0))**exp(-t/tau_P)*(b0 + b1*(a0 - a_star)))/b1

    a_exact
}
  

# Check function
tibble(
  t = seq(0, 3600, 1),
  P_init = 3,
  a0 = 150,
  b0 = 800,
  b1 = 2,
  a_init = a1(P_init, a0, b0, b1),
  P_star = 1,
  a_star = a1(P_star, a0, b0, b1),
  tau_P = 600,
  P = P_star - (P_star - P_init) * exp(-t / tau_P),
  a_t1 = a1(P, a0, b0, b1),
  a_t2 = at(t, a_init, a0, b0, b1, a_star, tau_P),
  a_t3 = at1(t, a_init, a0, b0, b1, a_star, tau_P)
) |>
  ggplot(aes(x = t)) +
  geom_line(aes(y = a_t1), color = "blue", size = 1) +
  geom_line(aes(y = a_t3), color = "red", linetype = "dashed", size = 1) +
  theme(
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14)
  ) +
  labs(y = "a(t)", title = "Comparison of a(t) from P(t) and direct integration")
  