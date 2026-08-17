# --- Ochoa et al. (2024) anatomical pore-conductance model ----------------
# (Plant Physiology, 10.1093/plphys/kiae292, Eq. 3), using Sack & Buckley
# (2016) non-grass shape factors (c = h = j = 0.5, kidney-shaped guard
# cells), the same convention already used for gmax in this study.

D_wv_ochoa <- 2.49e-5   # diffusivity of water vapor in air, m^2/s
v_ochoa    <- 2.24e-2   # molar volume of air, m^3/mol
c_SB <- 0.5; h_SB <- 0.5; j_SB <- 0.5   # Sack & Buckley non-grass constants

b_ochoa <- D_wv_ochoa / v_ochoa
m_ochoa <- pi * c_SB^2 / (sqrt(j_SB) * (4 * h_SB * j_SB + pi * c_SB))

mu_ochoa <- c_SB / sqrt(j_SB)      # Ochoa et al. 2024 mu, non-grass -> 1/sqrt(2)
j_ochoa  <- h_SB * j_SB / c_SB     # Ochoa et al. 2024 j, non-grass -> 0.5

# g_leaf,op as a function of fractional aperture (alpha = a / p)
gleaf_op <- function(alpha, d, s, gec, x_prime,
                      D_wa = D_wv_ochoa, V_ma = v_ochoa,
                      mu = mu_ochoa, j = j_ochoa) {
  gs  <- (D_wa * pi * alpha^2 * mu^2 * d * s) /
    (V_ma * (4 * j + pi * sqrt(alpha)) * (x_prime + alpha * mu * sqrt(s)))
  gns <- (2 - pi / 4) * alpha * mu^2 * d * s * gec
  gs + gns
}

# Numerically invert gleaf_op() for alpha given an observed conductance.
# Returns NA if the target conductance is not achievable within alpha in
# (0, 1] under the assumed gec/x_prime (i.e. exceeds the anatomical maximum).
solve_alpha <- function(g_target, d, s, gec, x_prime, mu = mu_ochoa, j = j_ochoa) {
  f <- function(alpha) gleaf_op(alpha, d, s, gec, x_prime, mu = mu, j = j) - g_target
  if (f(1) < 0) return(NA_real_)
  if (f(1e-8) > 0) return(NA_real_)
  uniroot(f, interval = c(1e-8, 1), tol = 1e-10)$root
}

# --- Series mesophyll/intercellular-airspace resistance ------------------
# (Kaiser 2009, Plant, Cell & Environment, 10.1111/j.1365-3040.2009.01990.x)
#
# Kaiser found that the pore-geometry-only relationship between aperture
# and gas exchange (as in gleaf_op() above) systematically underestimates
# how steeply conductance rises with aperture at small apertures, and that
# adding a diffusional resistance "in series" downstream of the stomatal
# pore (representing the substomatal cavity/intercellular airspace path,
# g_ias) substantially improved the fit. Physically, this additional
# resistance sits in series with the pore itself (both must be crossed by
# a diffusing molecule), but in parallel with (i.e. does not affect) the
# separate cuticular/epidermal pathway (gns), which bypasses the stomatal
# pore and substomatal cavity entirely. Hence:
#
#   1 / g_series(alpha) = 1 / gs(alpha) + 1 / g_ias
#   gleaf_op_series(alpha) = g_series(alpha) + gns(alpha)
#
# As g_ias -> Inf this reduces exactly to gleaf_op() above.
gleaf_op_series <- function(alpha, d, s, gec, x_prime, g_ias,
                             D_wa = D_wv_ochoa, V_ma = v_ochoa,
                             mu = mu_ochoa, j = j_ochoa) {
  gs  <- (D_wa * pi * alpha^2 * mu^2 * d * s) /
    (V_ma * (4 * j + pi * sqrt(alpha)) * (x_prime + alpha * mu * sqrt(s)))
  gns <- (2 - pi / 4) * alpha * mu^2 * d * s * gec
  g_series <- (gs * g_ias) / (gs + g_ias)
  g_series + gns
}

# Numerically invert gleaf_op_series() for alpha given an observed
# conductance and assumed g_ias. Returns NA if the target conductance is
# not achievable within alpha in (0, 1].
solve_alpha_series <- function(g_target, d, s, gec, x_prime, g_ias,
                                mu = mu_ochoa, j = j_ochoa) {
  f <- function(alpha) {
    gleaf_op_series(alpha, d, s, gec, x_prime, g_ias, mu = mu, j = j) - g_target
  }
  if (f(1) < 0) return(NA_real_)
  if (f(1e-8) > 0) return(NA_real_)
  uniroot(f, interval = c(1e-8, 1), tol = 1e-10)$root
}

# --- Phenomenological nonlinear aperture-conductance departure ------------
#
# Instead of a physically-motivated series resistor (gleaf_op_series() /
# Kaiser 2009 above), this represents a purely phenomenological departure
# from the Ochoa et al. (2024) pore-diffusion geometry: the pore-diffusion
# term gs is evaluated at an "effective" fractional aperture
#
#   alpha_eff = alpha / (1 + kappa * alpha)
#
# rather than at the true fractional aperture alpha. By construction,
# alpha_eff -> alpha to first order as alpha -> 0 (alpha_eff = alpha -
# kappa * alpha^2 + O(alpha^3)), so this exactly reproduces gleaf_op() at
# small aperture, but increasingly saturates below the Ochoa prediction as
# alpha grows (for kappa > 0), reaching alpha_eff = 1 / (1 + kappa) at full
# aperture (alpha = 1). kappa = 0 recovers gleaf_op() exactly. The
# cuticular/epidermal bypass term gns is left a function of the true alpha,
# since it does not pass through the stomatal pore and so has no reason to
# share the pore's nonlinearity.
#
# kappa is not measured; it is a free "how much departure" knob explored
# over a grid (see r/45_simulate-nonlinear-aperture-kinetics.R). See also
# gleaf_op_lehmann_or() below for a version where kappa is derived from
# stomatal geometry (density and pore size) rather than fitted freely.
gleaf_op_nonlinear <- function(alpha, d, s, gec, x_prime, kappa,
                                D_wa = D_wv_ochoa, V_ma = v_ochoa,
                                mu = mu_ochoa, j = j_ochoa) {
  alpha_eff <- alpha / (1 + kappa * alpha)
  gs  <- (D_wa * pi * alpha_eff^2 * mu^2 * d * s) /
    (V_ma * (4 * j + pi * sqrt(alpha_eff)) * (x_prime + alpha_eff * mu * sqrt(s)))
  gns <- (2 - pi / 4) * alpha * mu^2 * d * s * gec
  gs + gns
}

# Numerically invert gleaf_op_nonlinear() for alpha given an observed
# conductance and assumed kappa. Returns NA if the target conductance is
# not achievable within alpha in (0, 1].
solve_alpha_nonlinear <- function(g_target, d, s, gec, x_prime, kappa,
                                   mu = mu_ochoa, j = j_ochoa) {
  f <- function(alpha) {
    gleaf_op_nonlinear(alpha, d, s, gec, x_prime, kappa, mu = mu, j = j) - g_target
  }
  if (f(1) < 0) return(NA_real_)
  if (f(1e-8) > 0) return(NA_real_)
  uniroot(f, interval = c(1e-8, 1), tol = 1e-10)$root
}

# --- Geometry-derived crowding, informed by Lehmann & Or (2015) ----------
#
# Lehmann & Or (2015, New Phytologist 207:1015-1025) developed a resistor-
# network model in which interference between overlapping vapor
# concentration shells of neighboring stomata adds diffusive resistance,
# depending on stomatal spacing relative to pore size. Their exact model
# (periphery- vs interior-of-cluster resistances, an end-correction factor,
# a spacing-dependent vapor-shell resistance) needs cluster/spacing
# geometry we do not have for this dataset (only mean stomatal density);
# reproducing it exactly is not possible here.
#
# This function instead uses their qualitative mechanism -- crowding scales
# with the ratio of pore size to inter-stomatal spacing -- built from
# quantities already in this dataset, in place of the freely-chosen kappa
# in gleaf_op_nonlinear() / gleaf_op_crowding() above:
#   - Mean center-to-center stomatal spacing is approximated from density
#     alone, assuming a square lattice: lambda_s = 1 / sqrt(d).
#   - The instantaneous pore radius is approximated as alpha * amax / 2
#     (linear scaling with fractional aperture alpha, consistent with how
#     alpha already enters the pore-radius term of gleaf_op()'s gs
#     formula, i.e. "alpha * mu * sqrt(s)").
#   - chi(alpha) = pore diameter / spacing = 2 * alpha * (amax / 2) /
#     lambda_s is a dimensionless "how much of the inter-stomatal gap is
#     filled by the open pore" ratio -- 0 when closed, largest at kappa0 x
#     the true amax/spacing ratio at full aperture (median ~0.1, up to
#     ~0.2, across curves in this dataset; see r/49).
#   - As in gleaf_op_nonlinear(), the pore-diffusion term gs is evaluated
#     at alpha_eff = alpha / (1 + kappa0 * chi(alpha)), which reduces
#     algebraically to alpha / (1 + 2 * kappa0 * (amax / 2) * sqrt(d) *
#     alpha) -- i.e. the same rational saturating form as before, but with
#     the density/size dependence now derived from geometry rather than an
#     arbitrary function of d, and kappa0 left as the one unmeasured
#     "how many multiples of the geometric ratio actually matter"
#     calibration constant (kappa0 = 0 recovers gleaf_op() exactly).
gleaf_op_lehmann_or <- function(alpha, d, s, gec, x_prime, amax_m, kappa0,
                                 D_wa = D_wv_ochoa, V_ma = v_ochoa,
                                 mu = mu_ochoa, j = j_ochoa) {
  r_max <- amax_m / 2
  lambda_s <- 1 / sqrt(d)
  chi <- 2 * alpha * r_max / lambda_s
  alpha_eff <- alpha / (1 + kappa0 * chi)
  gs  <- (D_wa * pi * alpha_eff^2 * mu^2 * d * s) /
    (V_ma * (4 * j + pi * sqrt(alpha_eff)) * (x_prime + alpha_eff * mu * sqrt(s)))
  gns <- (2 - pi / 4) * alpha * mu^2 * d * s * gec
  gs + gns
}

# Numerically invert gleaf_op_lehmann_or() for alpha given an observed
# conductance and assumed kappa0. Returns NA if the target conductance is
# not achievable within alpha in (0, 1].
solve_alpha_lehmann_or <- function(g_target, d, s, gec, x_prime, amax_m, kappa0,
                                    mu = mu_ochoa, j = j_ochoa) {
  f <- function(alpha) {
    gleaf_op_lehmann_or(alpha, d, s, gec, x_prime, amax_m, kappa0, mu = mu, j = j) - g_target
  }
  if (f(1) < 0) return(NA_real_)
  if (f(1e-8) > 0) return(NA_real_)
  uniroot(f, interval = c(1e-8, 1), tol = 1e-10)$root
}

