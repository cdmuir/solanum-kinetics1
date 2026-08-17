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

# eqn 10 for g_leaf (= apparent g_sw)
# rearrange it to get r as a function g_leaf
# calculate implied r at gi and gf
# assume exponential change in r from gi to gf
# figure out exponential rate that gets as close to observed lambda and tau as possible