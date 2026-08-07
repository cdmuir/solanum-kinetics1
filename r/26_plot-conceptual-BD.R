# Conceptual figure — panels B and D only (modified)
source("r/header.R")

annotate <- ggplot2::annotate

t_seq <- seq(0, 400, by = 2)
gi0   <- 0.35
gf0   <- 0.04

base_theme <- theme_classic(base_size = 11) +
  theme(axis.title = element_text(size = 10))

# ---- Panel B: Effect of tau only (lambda = 1) ----
df_B <- expand_grid(tau = c(80, 220), t = t_seq) |>
  mutate(
    gsw     = weibull_gsw(t, gi = gi0, gf = gf0, tau = tau, lambda = 1),
    tau_lab = ifelse(tau == 80, "Small $\\tau$ (fast)", "Large $\\tau$ (slow)")
  )

pB <- ggplot(df_B, aes(t, gsw, color = tau_lab)) +
  geom_line(linewidth = 1) +
  scale_color_manual(
    values = c(
      "Small $\\tau$ (fast)" = col_fast,
      "Large $\\tau$ (slow)" = col_slow
    ),
    name = NULL
  ) +
  geom_hline(yintercept = gi0, linetype = "dotted", color = "gray55") +
  geom_hline(yintercept = gf0, linetype = "dotted", color = "gray55") +
  scale_y_continuous(limits = c(-0.01, 0.42), breaks = NULL) +
  scale_x_continuous(breaks = NULL) +
  labs(x = "time (s)", y = "stomatal conductance") +
  base_theme +
  theme(
    axis.text.x        = element_blank(),
    axis.ticks.x       = element_blank(),
    legend.position    = c(0.75, 0.75),
    legend.background  = element_blank(),
    legend.key.width   = unit(1.1, "cm"),
    legend.text        = element_text(size = 9)
  )

# ---- Panel D: Franks turgor-aperture curve (reactive vs. buffered) ----
K_D      <- 0.5
df_D     <- tibble(P = seq(0.01, 3.5, by = 0.01), g = P / (K_D + P))

fgmax_lo <- 0.25
fgmax_hi <- 0.70
P_lo     <- K_D * fgmax_lo / (1 - fgmax_lo)
P_hi     <- K_D * fgmax_hi / (1 - fgmax_hi)
sl_lo    <- K_D / (K_D + P_lo)^2
sl_hi    <- K_D / (K_D + P_hi)^2

dP   <- 0.25
tang <- tibble(
  grp       = c("low", "high"),
  P_mid     = c(P_lo, P_hi),
  g_mid     = c(fgmax_lo, fgmax_hi),
  slope     = c(sl_lo, sl_hi),
  intercept = g_mid - slope * P_mid,
  P         = P_mid - sqrt(dP^2 / (1 + slope^2)),
  P_max     = 2 * P_mid - P,
  g         = intercept + slope * P,
  g_max     = intercept + slope * P_max
)

pD <- ggplot() +
  geom_line(
    data = df_D, aes(P, g),
    linewidth = 1.2, color = "grey55"
  ) +
  geom_segment(
    data = filter(tang, grp == "low"),
    aes(P, g, xend = P_max, yend = g_max),
    color = col_fast, linewidth = 1.2
  ) +
  geom_segment(
    data = filter(tang, grp == "high"),
    aes(P, g, xend = P_max, yend = g_max),
    color = col_slow, linewidth = 1.2
  ) +
  geom_point(aes(x = P_lo, y = fgmax_lo), color = col_fast, size = 3) +
  geom_point(aes(x = P_hi, y = fgmax_hi), color = col_slow, size = 3) +
  # Low fsmax: reactive
  annotate(
    "text", x = 0.75, y = 0.16,
    label = "Low $f_\\mathrm{gmax}$: \\textbf{reactive}",
    color = col_fast, size = 2.9, hjust = 0
  ) +
  annotate(
    "text", x = 0.75, y = 0.07,
    label = "$~\\rightarrow$ steep slope, fast closure",
    color = col_fast, size = 2.9, hjust = 0
  ) +
  annotate(
    "segment",
    x = 0.73, xend = P_lo + 0.09,
    y = 0.20, yend = fgmax_lo - 0.01,
    color = col_fast, linewidth = 0.5,
    arrow = arrow(length = unit(0.06, "in"), type = "open")
  ) +
  # High fsmax: buffered
  annotate(
    "text", x = 1.9, y = 0.52,
    label = "High $f_\\mathrm{gmax}$: \\textbf{buffered}",
    color = col_slow, size = 2.9, hjust = 0
  ) +
  annotate(
    "text", x = 1.9, y = 0.44,
    label = "$~\\rightarrow$ shallow slope, slow closure",
    color = col_slow, size = 2.9, hjust = 0
  ) +
  annotate(
    "segment",
    x = 1.88, xend = P_hi + 0.06,
    y = 0.56, yend = fgmax_hi - 0.03,
    color = col_slow, linewidth = 0.5,
    arrow = arrow(length = unit(0.06, "in"), type = "open")
  ) +
  geom_hline(yintercept = 1, linetype = "dotted", color = "gray55") +
  annotate(
    "text", x = 3.0, y = 1.05,
    label = "maximum aperture", size = 3.5
  ) +
  labs(
    x = "Guard cell turgor pressure ($P_\\mathrm{g}$)",
    y = "stomatal conductance"
  ) +
  scale_x_continuous(breaks = NULL) +
  scale_y_continuous(breaks = NULL) +
  coord_cartesian(ylim = c(0, 1.12), xlim = c(0, 3.6)) +
  base_theme

# ---- Combine and save ----
gp_BD <- pB | pD

tikz(
  "figures/conceptual_BD.tex",
  standAlone = TRUE,
  width = 180 * 0.0393701,
  height = 80 * 0.0393701
)
print(gp_BD)
dev.off()

system("cd figures; pdflatex conceptual_BD.tex; rm conceptual_BD.aux conceptual_BD.log")
