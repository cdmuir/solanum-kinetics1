# --- Conceptual figure illustrating the two hypotheses linking anatomy to kinetics ---

source("r/header.R")

t_seq  <- seq(0, 400, by = 2)
gi0    <- 0.35
gf0    <- 0.04
tau_A  <- 150
lambda_A <- 2.8

base_theme <- theme_classic(base_size = 11) +
  theme(axis.title = element_text(size = 10))

# ---- Panel A: Annotated example time course with WWR + RWR ----
g_pre  <- 0.27
g_peak <- g_pre + 0.03  # subtle WWR overshoot, not the extreme bump used previously
t_step <- 0
t_wwr  <- 30            # short WWR duration

t_pre     <- seq(-100, -1, by = 2)
t_wwr_seq <- seq(0, t_wwr, by = 2)
t_rwr_seq <- seq(t_wwr + 2, 400, by = 2)

# Data logging did not begin at the WWR peak; it began once gsw declined
# back down through its pre-step value (g_pre), well into the RWR decline.
# Solve for that crossing time along the (still continuous, unlogged) true
# RWR trajectory so Panel A can mark the actual start of logging (t = 0 on
# the plotted axis, matching how curves are time-zeroed in the real data).
t_log_start <- t_wwr + tau_A * (-log((g_pre - gf0) / (g_peak - gf0))) ^ (1 / lambda_A)
t_logged_seq <- seq(t_log_start, 400, by = 2)

# Annotation heights, scaled relative to the (now shorter) WWR peak
y_bracket   <- g_peak + 0.04
y_bracket_label <- g_peak + 0.065
y_guide_top <- g_peak + 0.06
y_top_lim   <- g_peak + 0.09

df_A2 <- bind_rows(
  tibble(t = t_pre, gsw = g_pre, phase = "Pre-step"),
  tibble(
    t = t_wwr_seq,
    gsw = wwr_fun(t_wwr_seq),
    phase = "WWR"
  ),
  tibble(
    t = t_rwr_seq,
    gsw = rwr_fun(t_rwr_seq),
    phase = "RWR"
  )
)
df_A2_model <- tibble(t = t_logged_seq, gsw = rwr_fun(t_logged_seq))

tau_global   <- t_wwr + tau_A
gsw_at_tau2  <- rwr_fun(tau_global)

pA <- ggplot(df_A2, aes(t, gsw)) +
  annotate(
    "segment",
    x = t_step,
    xend = t_step,
    y = 0,
    yend = y_guide_top,
    linetype = "solid",
    color = "gray70",
    linewidth = 0.6
  ) +
  annotate(
    "text",
    x = t_step + 5,
    y = 0,
    label = "VPD step",
    hjust = 1.2,
    vjust = 0,
    size = 2.8,
    color = "gray50"
  ) +
  geom_line(linewidth = 1, color = "gray20") +
  geom_line(
    data = df_A2_model,
    aes(t, gsw),
    color = col_fast,
    linewidth = 1.1,
    inherit.aes = FALSE
  ) +
  geom_hline(yintercept = g_pre,
             linetype = "dotted",
             color = "gray55") +
  geom_hline(yintercept = gf0,
             linetype = "dotted",
             color = "gray55") +
  annotate(
    "segment",
    x = tau_global,
    xend = tau_global,
    y = 0,
    yend = gsw_at_tau2,
    linetype = "dashed",
    color = "gray40",
    linewidth = 0.7
  ) +
  annotate(
    "segment",
    x = t_log_start,
    xend = t_log_start,
    y = 0,
    yend = g_peak,
    color = col_fast,
    linewidth = 0.7
  ) +
  annotate(
    "text",
    x = t_log_start,
    y = g_peak,
    label = "Log data",
    hjust = 0.5,
    vjust = -1,
    size = 2.8,
    color = col_fast
  ) +
  geom_point(
    data = tibble(x = tau_global, y = gsw_at_tau2),
    aes(x, y),
    color = "gray30",
    size = 2.5,
    inherit.aes = FALSE
  ) +
  annotate(
    "text",
    x = -90,
    y = g_pre + 0.025,
    label = "$g_\\mathrm{i}$",
    size = 3.8,
    hjust = 1
  ) +
  annotate(
    "text",
    x = -90,
    y = gf0 + 0.025,
    label = "$g_\\mathrm{f}$",
    size = 3.8,
    hjust = 1
  ) +
  annotate(
    "text",
    x = tau_global + 25,
    y = 0.003,
    label = "$\\tau$",
    vjust = 0,
    hjust = 0,
    size = 4.5,
    color = "gray40"
  ) +
  annotate(
    "segment",
    x = t_step + 2,
    xend = t_wwr - 2,
    y = y_bracket,
    yend = y_bracket,
    color = "gray50",
    arrow = arrow(
      ends = "both",
      length = unit(0.05, "in"),
      type = "open"
    )
  ) +
  annotate(
    "text",
    x = (t_step + t_wwr) / 2,
    y = y_bracket_label,
    label = "WWR",
    size = 2.8,
    color = "gray50",
    hjust = 0.1,
    vjust = 1,
    angle = 33
  ) +
  # RWR begins where gsw starts to decline (the WWR peak), not where
  # logging/the Weibull model fit begins (marked separately below by the
  # "Log data" line at t_log_start).
  annotate(
    "segment",
    x = t_wwr + 2,
    xend = 380,
    y = y_bracket,
    yend = y_bracket,
    color = "gray50",
    arrow = arrow(
      ends = "both",
      length = unit(0.05, "in"),
      type = "open"
    )
  ) +
  annotate(
    "text",
    x = (t_wwr + 380) / 2,
    y = y_bracket_label,
    label = "RWR",
    size = 2.8,
    color = "gray50",
    hjust = 0.5
  ) +
  annotate(
    "text",
    x = 140,
    y = g_peak - 0.05,
    label = "$g_\\mathrm{sw,t} = g_\\mathrm{f} + (g_\\mathrm{i} - g_\\mathrm{f}) e^{-\\left(\\frac{t}{\\tau}\\right)^\\lambda}$",
    size = 2.9,
    hjust = 0,
    color = col_fast
  ) +
  scale_y_continuous(limits = c(0, y_top_lim), breaks = NULL) +
  scale_x_continuous(breaks = seq(0, 400, by = 100) + t_log_start,
                     labels = seq(0, 400, by = 100)) +
  labs(x = "time (s)", y = "stomatal conductance") +
  base_theme

# ---- Panel B: Effect of tau and lambda ----
df_B <- expand_grid(tau = c(80, 220),
                    lambda = c(1.2, 3.5),
                    t = t_seq) |>
  mutate(
    gsw        = weibull_gsw(
      t,
      gi = gi0,
      gf = gf0,
      tau = tau,
      lambda = lambda
    ),
    tau_lab    = ifelse(tau == 80, "Small $\\tau$ (fast)", "Large $\\tau$ (slow)"),
    lambda_lab = ifelse(lambda < 2, "Small $\\lambda$", "Large $\\lambda$ (lag)")
  )

df_B_legend <- distinct(df_B, tau_lab) |>
  mutate(t = 0, gsw = 0)

pB <- ggplot(df_B, aes(t, gsw, color = tau_lab, linetype = lambda_lab)) +
  # show.legend suppresses only the color legend on this layer -- the
  # default legend key glyph for a line geom's color aesthetic is a
  # colored line segment, which reviewers found distracting/ambiguous
  # next to the linetype legend for lambda. The color (tau) legend is
  # instead drawn from the invisible geom_point layer below, whose
  # default point glyph is overridden to a solid square swatch.
  geom_line(linewidth = 1, show.legend = c(color = FALSE, linetype = TRUE)) +
  geom_point(
    data = df_B_legend,
    aes(t, gsw, color = tau_lab),
    inherit.aes = FALSE,
    alpha = 0,
    show.legend = TRUE
  ) +
  scale_color_manual(
    values = c(
      "Small $\\tau$ (fast)" = col_fast,
      "Large $\\tau$ (slow)" = col_slow
    ),
    name = NULL,
    guide = guide_legend(override.aes = list(alpha = 1, shape = 15, size = 5))
  ) +
  scale_linetype_manual(
    values = c(
      "Small $\\lambda$" = "solid",
      "Large $\\lambda$ (lag)" = "dashed"
    ),
    name = NULL
  ) +
  geom_hline(yintercept = gi0,
             linetype = "dotted",
             color = "gray55") +
  geom_hline(yintercept = gf0,
             linetype = "dotted",
             color = "gray55") +
  scale_y_continuous(limits = c(-0.01, 0.42), breaks = NULL) +
  labs(x = "time (s)", y = "stomatal conductance") +
  base_theme +
  theme(
    legend.position    = c(0.75, 0.75),
    legend.background  = element_blank(),
    legend.key.width   = unit(1.1, "cm"),
    legend.text        = element_text(size = 9)
  )

# ---- Panel C: Guard cell size -> SA:V -> kinetics ----
df_C <- expand_grid(
  cell_size = c("Small guard cells (high SA:V)", "Large guard cells (low SA:V)"),
  t = t_seq
) |>
  mutate(
    tau_c     = ifelse(cell_size == "Small guard cells (high SA:V)", 80, 220),
    gsw       = weibull_gsw(
      t,
      gi = gi0,
      gf = gf0,
      tau = tau_c,
      lambda = 2.8
    ),
    cell_size = factor(
      cell_size,
      levels = c("Large guard cells (low SA:V)", "Small guard cells (high SA:V)")
    )
  )

pC <- ggplot(df_C, aes(t, gsw, color = cell_size)) +
  geom_line(linewidth = 1) +
  scale_color_manual(
    values = c(
      "Small guard cells (high SA:V)" = col_fast,
      "Large guard cells (low SA:V)" = col_slow
    ),
    name = NULL
  ) +
  scale_fill_manual(values = c("small" = col_fast, "large" = col_slow),
                    guide = "none") +
  annotate(
    "text",
    x = 300,
    y = 0.265,
    label = "SA:V $\\propto$ 1/radius",
    parse = FALSE,
    size = 3,
    hjust = 0.5,
    color = "gray30"
  ) +
  annotate(
    "text",
    x = 300,
    y = 0.225,
    label = "faster osmolyte flux",
    size = 2.9,
    hjust = 0.5,
    color = "gray30"
  ) +
  annotate(
    "text",
    x = 300,
    y = 0.185,
    label = "faster turgor change",
    size = 2.9,
    hjust = 0.5,
    color = "gray30"
  ) +
  geom_hline(yintercept = gi0,
             linetype = "dotted",
             color = "gray55") +
  geom_hline(yintercept = gf0,
             linetype = "dotted",
             color = "gray55") +
  scale_y_continuous(limits = c(-0.01, 0.42), breaks = NULL) +
  scale_x_continuous(limits = c(0, 400)) +
  labs(x = "time (s)", y = "stomatal conductance") +
  base_theme +
  theme(
    legend.position   = c(0.75, 0.9),
    legend.background = element_blank(),
    legend.text       = element_text(size = 9)
  )

# ---- Panel D: Franks turgor-aperture curve ----
K_D <- 0.5
df_D <- tibble(P = seq(0.01, 3.5, by = 0.01), g = P / (K_D + P))

fgmax_lo <- 0.25
fgmax_hi <- 0.70
P_lo <- K_D * fgmax_lo / (1 - fgmax_lo)
P_hi <- K_D * fgmax_hi / (1 - fgmax_hi)
sl_lo <- K_D / (K_D + P_lo)^2
sl_hi <- K_D / (K_D + P_hi)^2

dP <- 0.25 # length of tangent line
tang = tibble(
  grp = c("low", "high"),
  P_mid = c(P_lo, P_hi),
  g_mid = c(fgmax_lo, fgmax_hi),
  slope = c(sl_lo, sl_hi),
  intercept = g_mid - slope * P_mid,
  P = P_mid - sqrt(dP^2 / (1 + slope^2)),
  P_max = 2 * P_mid - P,
  g = intercept + slope * P,
  g_max = intercept + slope * P_max
)


pD <- ggplot() +
  geom_line(
    data = df_D,
    aes(P, g),
    linewidth = 1.2,
    color = "grey55"
  ) +
  geom_segment(
    data = filter(tang, grp == "low"),
    aes(P, g, xend = P_max, yend = g_max),
    color = col_fast,
    linewidth = 1.2
  ) +
  geom_segment(
    data = filter(tang, grp == "high"),
    aes(P, g, xend = P_max, yend = g_max),
    color = col_slow,
    linewidth = 1.2
  ) +
  scale_color_manual(values = c("low" = col_fast, "high" = col_slow)) +
  geom_point(aes(x = P_lo, y = fgmax_lo), color = col_fast, size = 3) +
  geom_point(aes(x = P_hi, y = fgmax_hi), color = col_slow, size = 3) +
  annotate(
    "text",
    x = 0.75,
    y = 0.13,
    label = "Low $f_\\mathrm{gmax}$ (steep slope)",
    color = col_fast,
    size = 2.9,
    hjust = 0
  ) +
  annotate(
    "text",
    x = 0.75,
    y = 0.07,
    label = "$~\\rightarrow$ fast closure",
    color = col_fast,
    size = 2.9,
    hjust = 0
  ) +
  annotate(
    "segment",
    x = 0.73,
    xend = P_lo + 0.09,
    y = 0.17,
    yend = fgmax_lo - 0.01,
    color = col_fast,
    linewidth = 0.5,
    arrow = arrow(length = unit(0.06, "in"), type = "open")
  ) +
  annotate(
    "text",
    x = 1.9,
    y = 0.50,
    label = "High $f_\\mathrm{gmax}$ (shallow slope)",
    color = col_slow,
    size = 2.9,
    hjust = 0
  ) +
  annotate(
    "text",
    x = 1.9,
    y = 0.44,
    label = "$~\\rightarrow$ slow closure",
    color = col_slow,
    size = 2.9,
    hjust = 0
  ) +
  annotate(
    "segment",
    x = 1.88,
    xend = P_hi + 0.06,
    y = 0.54,
    yend = fgmax_hi - 0.03,
    color = col_slow,
    linewidth = 0.5,
    arrow = arrow(length = unit(0.06, "in"), type = "open")
  ) +
  geom_hline(yintercept = 1,
             linetype = "dotted",
             color = "gray55") +
  annotate(
    "text",
    x = 3.0,
    y = 1.05,
    label = "maximum aperture",
    size = 3.5
  ) +
  labs(x = "Guard cell turgor pressure ($P_\\mathrm{g}$)", y = "stomatal conductance") +
  scale_x_continuous(breaks = NULL) +
  scale_y_continuous(breaks = NULL) +
  coord_cartesian(ylim = c(0, 1.12), xlim = c(0, 3.6)) +
  base_theme

# ---- Combine and save ----
gp1 <- (pA | pB) / (pC | pD) +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(face = "bold"))

tikz(
  "figures/conceptual.tex",
  standAlone = TRUE,
  width = 180 * 0.0393701,
  height = 160 * 0.0393701
)
print(gp1)
dev.off()

system("cd figures; pdflatex conceptual.tex; rm conceptual.aux conceptual.log")
