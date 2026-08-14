# Calculate the theoretical maximum stomatal pore aperture (a_max) from
# guard cell length and stomatal density, following the same allometric
# assumptions used for gmax in this study (Sack & Buckley 2016, non-grass
# kidney-shaped guard cells: c = h = j = 0.5). Then invert Ochoa et al.
# (2024) Eq. 3 (g_leaf,op as a function of fractional aperture alpha) to
# find the aperture implied by the observed gi (ginit_mean) and gf
# (gfinal_mean) for each curve, across a grid of plausible cuticular
# conductance (gec) and mean free path (x_prime) values, since we do not
# have direct measurements of either.
#
# Key relationships (see Notes/response-to-reviewers for full derivation):
#   - p (pore length) = c * l_gc; a_max = p (circular pore at full aperture)
#   - s (stomatal area) = l_gc^2 / 2, consistent with this study's
#     gmax = b * m * d * l_gc / sqrt(2) (i.e. sqrt(s) = l_gc / sqrt(2))
#   - Total stomatal density d is recovered by inverting the study's own
#     gmax formula (gmax and l_gc are already stored in data/joined-summary.rds;
#     per-surface densities were not retained in r/01_join-data.R), rather
#     than re-deriving it from raw anatomy data. This is exact: at
#     gec = x_prime = 0, gleaf_op(alpha = 1) reproduces the stored gmax to
#     within floating-point error.
#   - Ochoa et al. (2024) mu and j (non-grass) are algebraically equivalent
#     to Sack & Buckley's c, h, j = 0.5 (mu = c/sqrt(j) = 1/sqrt(2),
#     j_Ochoa = h*j/c = 0.5), so the same shape-factor convention applies.

source("r/header.R")

# Constants, gleaf_op(), and solve_alpha() (Ochoa et al. 2024 Eq. 3, with
# Sack & Buckley 2016 non-grass shape factors) are defined in r/functions.R
# so they can be reused consistently across scripts (e.g. r/41_*, r/42_*).

joined_summary <- read_rds("data/joined-summary.rds")

# Per-curve anatomical quantities derived from stored gmax and guard cell length
anat <- joined_summary |>
  transmute(
    id, accession, replicate, leaf_type, light_intensity,
    ginit_mean, gfinal_mean, gmax,
    l_gc_m  = guard_cell_length_um * 1e-6,
    p_m     = c_SB * l_gc_m,                      # pore length = max aperture
    amax_um = p_m * 1e6,
    s_m2    = l_gc_m^2 / 2,
    d_m2    = gmax * sqrt(2) / (b_ochoa * m_ochoa * l_gc_m)   # recovered total stomatal density, m^-2
  )

# Grid of plausible cuticular conductance and mean free path values, since
# neither was measured directly.
#   gec: 0 (no cuticular conductance), then representative low/moderate/high
#     values spanning the compiled literature range (e.g. Kerstiens 1996),
#     mol m^-2 s^-1.
#   x_prime: 0 (no correction), then ~1x and ~2x the mean free path of air
#     molecules at ~20-25 degC, sea level (~6.8e-8 m), m.
gec_grid     <- c(0, 0.005, 0.02, 0.05)
xprime_grid  <- c(0, 6.8e-8, 1.36e-7)

grid <- expand_grid(gec = gec_grid, x_prime = xprime_grid)

tbl_implied_aperture <- anat |>
  cross_join(grid) |>
  mutate(
    alpha_gi = pmap_dbl(list(ginit_mean, d_m2, s_m2, gec, x_prime), solve_alpha),
    alpha_gf = pmap_dbl(list(gfinal_mean, d_m2, s_m2, gec, x_prime), solve_alpha),
    a_gi_um  = alpha_gi * amax_um,
    a_gf_um  = alpha_gf * amax_um
  )

write_rds(tbl_implied_aperture, "objects/tbl-implied-aperture.rds")
