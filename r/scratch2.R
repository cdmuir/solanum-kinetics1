# Deriving model with da / dP as power function
library(cowplot)
library(dplyr)
library(ggplot2)
library(tibble)

theme_set(theme_cowplot())

# dP / dt = (P* - P) / tau_P
# P(t) = P* - (P* - P0) * exp(-t / tau_P)
# 
# da / dP = k * (a - a_max) / (1 + P)

tibble(
  P = seq(0, 6, length.out = 100),
  a0 = 150, a_max = 950, k = -1, a = 250,
  da_dP = k * (a - a_max) / (1 + P)
) |>
  ggplot(aes(x = P, y = da_dP)) +
  geom_line(size = 1, color = "blue") +
  theme(
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14)
  ) 

# Assume a(P = 0) = a0
# a(P) = a_max + (P + 1)**k*(-a_max + a0)

a1 = function(P, a_max, a0, k) {
  a_max + (P + 1)**k*(-a_max + a0)
}

# Check function
tibble(P = seq(0, 6, 0.01), a1 = a1(P, 950, 150, -2)) |>
  ggplot(aes(x = P, y = a1)) +
  geom_line(size = 1, color = "blue") +
  theme(
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14)
  ) +
  ylim(150, 950)

# By chain rule:
# da / dP = (da / dt) / (dP / dt) = k * (a - a_max) / (1 + P)
# da / dt = (dP / dt) * k * (a - a_max) / (1 + P)
# 
# Subbing in dP / dt:
# da / dt = (P* - P) / tau_P * k * (a - a_max) / (1 + P)
# 
# Rearrange equation above to get P in terms of a (assume a0 < a < a_max, k < 0):
# P(a) = ((a - a_max)/(a0 - a_max))**(1/k) - 1

P1 = function(a, a_max, a0, k) {
  ((a - a_max)/(a0 - a_max))**(1/k) - 1
}

# Check function
a = a1(P = 3, a_max = 950, a0 = 150, k = -2)
P1(a = a, a_max = 950, a0 = 150, k = -2)

# Use this equation to put da / dt in terms of a:
# P* - P = (((a_star - a_max)/(a0 - a_max))**(1/k) - 1) - (((a - a_max)/(a0 - a_max))**(1/k) - 1)

# Simplify
P_star <- 1
a0 = 150
a_max = 950
k = -runif(1)
a_star = a1(P_star, a_max, a0, k)

(((a_star - a_max)/(a0 - a_max))**(1/k) - 1) - (((a - a_max)/(a0 - a_max))**(1/k) - 1)
-((a - a_max)/(a0 - a_max))**(1/k) + (-(a_max - a_star)/(a0 - a_max))**(1/k)
(-(a_max - a_star)/(a0 - a_max))**(1/k) - ((a - a_max)/(a0 - a_max))**(1/k)

# Therefore:
# P* - P = (-(a_max - a_star)/(a0 - a_max))**(1/k) - ((a - a_max)/(a0 - a_max))**(1/k)

# k * (a - a_max) / (1 + P)

P = P1(a = a, a_max = a_max, a0 = a0, k = k)

k * (a - a_max) / (1 + P)
k * (a - a_max) / (1 + (((a - a_max)/(a0 - a_max))**(1/k) - 1))
k*(a - a_max)/((a - a_max)/(a0 - a_max))**(1/k)

# Therefore:
# k * (a - a_max) / (1 + P) = k*(a - a_max)/((a - a_max)/(a0 - a_max))**(1/k)

# Putting these together:
# da / dt = ((-(a_max - a_star)/(a0 - a_max))**(1/k) - ((a - a_max)/(a0 - a_max))**(1/k)) / tau_P * k*(a - a_max)/((a - a_max)/(a0 - a_max))**(1/k)

# Check
P <- 3
P_star <- 1
a0 = 150
a_max = 950
k = -runif(1)
a_star = a1(P_star, a_max, a0, k)
a = a1(P, a_max, a0, k)
tau_P = 600

(P_star - P) / tau_P * k * (a - a_max) / (1 + P)
((-(a_max - a_star)/(a0 - a_max))**(1/k) - ((a - a_max)/(a0 - a_max))**(1/k)) / tau_P * k*(a - a_max)/((a - a_max)/(a0 - a_max))**(1/k)
-k*(a - a_max)*(((a - a_max)/(a0 - a_max))**(1/k) - (-(a_max - a_star)/(a0 - a_max))**(1/k))/(tau_P*((a - a_max)/(a0 - a_max))**(1/k))
(k / tau_P) * (a_max - a) * (1 - ((a_max - a) / (a_max - a_star)) ^ (-1/k))

# simplified - not correct for some reason
at1 = function(t, a_init, a0, a_max, k, a_star, tau_P) {
  
  (k / tau_P) * (a_max - a_init) * (1 - ((a_max - a_init) / (a_max - a_star)) ^ (-1/k))
  
}

# different approach (see python, just using a(P) to get P(t) in terms of a)
at1 = function(t, a_init, a0, a_max, k, a_star, tau_P) {
  a_max + ((((a_init - a_max)/(a0 - a_max))**(1/k) - ((-a_max + a_star)/(a0 - a_max))**(1/k) + (((-a_max + a_star)/(a0 - a_max))**(1/k) - 1)*exp(t/tau_P) + exp(t/tau_P))*exp(-t/tau_P))**k*(a0 - a_max)
  }  

# Check function
tibble(
  t = seq(0, 3600, 1),
  P_init = 3,
  a0 = 150,
  a_max = 950,
  k = -2,
  a_init = a1(P_init, a_max, a0, k),
  P_star = 1,
  a_star = a1(P_star, a_max, a0, k),
  tau_P = 600,
  P = P_star - (P_star - P_init) * exp(-t / tau_P),
  a_t1 = a1(P, a_max, a0, k),
  a_t2 = at1(t, a_init, a0, a_max, k, a_star, tau_P),
) |>
  ggplot(aes(x = t)) +
  geom_line(aes(y = a_t1), color = "blue", size = 1) +
  geom_line(aes(y = a_t2), color = "red", linetype = "dashed", size = 1) +
  theme(
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14)
  ) +
  labs(y = "a(t)", title = "Comparison of a(t) from P(t) and direct integration")
  