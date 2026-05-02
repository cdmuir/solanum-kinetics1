# Conceptual figure illustrating stomatal closure kinetics
source("r/header.R")

# ggpp (loaded via header.R) masks ggplot2::annotate; restore the ggplot2 version since this script relies on parse = TRUE support
annotate <- ggplot2::annotate

t_seq  <- seq(0, 400, by = 2)
gi0    <- 0.35
gf0    <- 0.04
tau_A  <- 150
lambda_A <- 2.8

base_theme <- theme_classic(base_size = 11) +
  theme(axis.title = element_text(size = 10))

# ---- Panel A: Annotated example time course with WWR + RWR ----
g_pre  <- 0.27
g_peak <- gi0
t_step <- 0
t_wwr  <- 65

t_pre     <- seq(-100, -1, by = 2)
t_wwr_seq <- seq(0, t_wwr, by = 2)
t_rwr_seq <- seq(t_wwr + 2, 400, by = 2)

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
df_A2_model <- tibble(t = t_rwr_seq, gsw = rwr_fun(t_rwr_seq))

tau_global   <- t_wwr + tau_A
gsw_at_tau2  <- rwr_fun(tau_global)

pA <- ggplot(df_A2, aes(t, gsw)) +
  annotate(
    "segment",
    x = t_step,
    xend = t_step,
    y = 0,
    yend = 0.41,
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
    linewidth = 1,
    inherit.aes = FALSE
  ) +
  geom_hline(yintercept = g_peak,
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
    x = t_wwr,
    xend = t_wwr,
    y = 0,
    yend = 0.41,
    color = col_fast,
    linewidth = 0.7
  ) +
  annotate(
    "text",
    x = t_wwr,
    y = 0,
    label = "Log data",
    hjust = -0.2,
    vjust = 0,
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
    y = g_peak + 0.025,
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
    x = tau_global + 8,
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
    y = 0.39,
    yend = 0.39,
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
    y = 0.415,
    label = "WWR",
    size = 2.8,
    color = "gray50",
    hjust = 0.5
  ) +
  annotate(
    "segment",
    x = t_wwr + 2,
    xend = 380,
    y = 0.39,
    yend = 0.39,
    color = col_fast,
    arrow = arrow(
      ends = "both",
      length = unit(0.05, "in"),
      type = "open"
    )
  ) +
  annotate(
    "text",
    x = (t_wwr + 380) / 2,
    y = 0.415,
    label = "RWR (Weibull model)",
    size = 2.8,
    color = col_fast,
    hjust = 0.5
  ) +
  annotate(
    "text",
    x = 160,
    y = 0.315,
    label = "$g_\\mathrm{sw,t} = g_\\mathrm{f} + (g_\\mathrm{i} - g_\\mathrm{f}) e^{-\\left(\\frac{t}{\\tau}\\right)^\\lambda}$",
    size = 2.9,
    hjust = 0,
    color = col_fast
  ) +
  scale_y_continuous(limits = c(0, 0.44), breaks = NULL) +
  scale_x_continuous(breaks = seq(0, 400, by = 100) + t_wwr,
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

pB <- ggplot(df_B, aes(t, gsw, color = tau_lab, linetype = lambda_lab)) +
  geom_line(linewidth = 1) +
  scale_color_manual(
    values = c(
      "Small $\\tau$ (fast)" = col_fast,
      "Large $\\tau$ (slow)" = col_slow
    ),
    name = NULL
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
  plot_annotation(tag_levels = "A")

# ggsave("figures/conceptual.pdf", width = 180, height = 160, units = "mm")

tikz(
  "figures/conceptual.tex",
  standAlone = TRUE,
  width = 180 * 0.0393701,
  height = 160 * 0.0393701
)
print(gp1)
dev.off()

system("cd figures; pdflatex conceptual.tex; rm conceptual.aux conceptual.log")
