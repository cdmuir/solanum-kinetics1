# Conceptual figure illustrating stomatal closure kinetics
source("r/header.R")

# ggpp (loaded via header.R) masks ggplot2::annotate; restore the ggplot2 version
# since this script relies on parse = TRUE support
annotate <- ggplot2::annotate

# ---- Helper ----
weibull_gsw <- function(t, gi, gf, tau, lambda) {
  gf + (gi - gf) * exp(-(t / tau)^lambda)
}

col_fast <- "#0072B2"
col_slow <- "#D55E00"

t_seq  <- seq(0, 700, by = 2)
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
t_wwr  <- 55

wwr_fun <- function(t) {
  g_pre + (g_peak - g_pre) * (1 - exp(-t / (t_wwr / 3)))
}
rwr_fun <- function(t) {
  weibull_gsw(t - t_wwr, gi = g_peak, gf = gf0, tau = tau_A, lambda = lambda_A)
}

t_pre     <- seq(-100, -1,    by = 2)
t_wwr_seq <- seq(0,    t_wwr, by = 2)
t_rwr_seq <- seq(t_wwr + 2, 700, by = 2)

df_A2 <- bind_rows(
  tibble(t = t_pre,     gsw = g_pre,               phase = "Pre-step"),
  tibble(t = t_wwr_seq, gsw = wwr_fun(t_wwr_seq),  phase = "WWR"),
  tibble(t = t_rwr_seq, gsw = rwr_fun(t_rwr_seq),  phase = "RWR")
)
df_A2_model <- tibble(t = t_rwr_seq, gsw = rwr_fun(t_rwr_seq))

tau_global   <- t_wwr + tau_A
gsw_at_tau2  <- rwr_fun(tau_global)

pA <- ggplot(df_A2, aes(t, gsw)) +
  annotate("segment", x = t_step, xend = t_step, y = 0, yend = 0.41,
           linetype = "solid", color = "gray70", linewidth = 0.6) +
  annotate("text", x = t_step + 5, y = 0.40,
           label = "VPD\nstep", hjust = 0, vjust = 1, size = 2.8, color = "gray50") +
  geom_line(linewidth = 1, color = "gray20") +
  geom_line(data = df_A2_model, aes(t, gsw),
            color = col_fast, linewidth = 1, linetype = "dashed", inherit.aes = FALSE) +
  geom_hline(yintercept = g_peak, linetype = "dotted", color = "gray55") +
  geom_hline(yintercept = gf0,    linetype = "dotted", color = "gray55") +
  annotate("segment", x = tau_global, xend = tau_global, y = 0, yend = gsw_at_tau2,
           linetype = "dashed", color = "gray40", linewidth = 0.7) +
  geom_point(data = tibble(x = tau_global, y = gsw_at_tau2), aes(x, y),
             color = "gray30", size = 2.5, inherit.aes = FALSE) +
  annotate("text", x = 685, y = g_peak + 0.009,
           label = "g[i]", parse = TRUE, size = 3.8, hjust = 1) +
  annotate("text", x = 685, y = gf0 + 0.009,
           label = "g[f]", parse = TRUE, size = 3.8, hjust = 1) +
  annotate("text", x = tau_global + 8, y = 0.003, label = "tau",
           parse = TRUE, vjust = 0, hjust = 0, size = 4.5, color = "gray40") +
  annotate("segment", x = t_step + 2, xend = t_wwr - 2, y = 0.39, yend = 0.39,
           color = "gray50",
           arrow = arrow(ends = "both", length = unit(0.05, "in"), type = "open")) +
  annotate("text", x = (t_step + t_wwr) / 2, y = 0.415,
           label = "WWR", size = 2.8, color = "gray50", hjust = 0.5) +
  annotate("segment", x = t_wwr + 2, xend = 680, y = 0.39, yend = 0.39,
           color = col_fast,
           arrow = arrow(ends = "both", length = unit(0.05, "in"), type = "open")) +
  annotate("text", x = (t_wwr + 680) / 2, y = 0.415,
           label = "RWR (Weibull model)", size = 2.8, color = col_fast, hjust = 0.5) +
  annotate("text", x = 350, y = 0.245,
           label = "g[sw*','*t]==g[f]+(g[i]-g[f])~e^{-(t/tau)^lambda}",
           parse = TRUE, size = 2.9, hjust = 0, color = col_fast) +
  scale_y_continuous(limits = c(-0.01, 0.44)) +
  scale_x_continuous(breaks = seq(-100, 700, by = 100)) +
  labs(x = "Time (s)", y = expression(g[sw]~(mol~m^{-2}~s^{-1}))) +
  base_theme

# ---- Panel B: Effect of tau and lambda ----
df_B <- expand_grid(tau = c(80, 220), lambda = c(1.2, 3.5), t = t_seq) |>
  mutate(
    gsw        = weibull_gsw(t, gi = gi0, gf = gf0, tau = tau, lambda = lambda),
    tau_lab    = ifelse(tau == 80,  "Small \u03c4 (fast)", "Large \u03c4 (slow)"),
    lambda_lab = ifelse(lambda < 2, "Small \u03bb",        "Large \u03bb (lag)")
  )

pB <- ggplot(df_B, aes(t, gsw, color = tau_lab, linetype = lambda_lab)) +
  geom_line(linewidth = 1) +
  scale_color_manual(
    values = c("Small \u03c4 (fast)" = col_fast, "Large \u03c4 (slow)" = col_slow),
    name = NULL
  ) +
  scale_linetype_manual(
    values = c("Small \u03bb" = "solid", "Large \u03bb (lag)" = "dashed"),
    name = NULL
  ) +
  geom_hline(yintercept = gi0, linetype = "dotted", color = "gray55") +
  geom_hline(yintercept = gf0, linetype = "dotted", color = "gray55") +
  scale_y_continuous(limits = c(-0.01, 0.42)) +
  labs(x = "Time (s)", y = expression(g[sw]~(mol~m^{-2}~s^{-1}))) +
  base_theme +
  theme(
    legend.position    = c(0.68, 0.68),
    legend.background  = element_blank(),
    legend.key.width   = unit(1.1, "cm"),
    legend.text        = element_text(size = 9)
  )

# ---- Panel C: Guard cell size -> SA:V -> kinetics ----
df_C <- expand_grid(
  cell_size = c("Small guard cells", "Large guard cells"),
  t = t_seq
) |>
  mutate(
    tau_c     = ifelse(cell_size == "Small guard cells", 80, 220),
    gsw       = weibull_gsw(t, gi = gi0, gf = gf0, tau = tau_c, lambda = 2.8),
    cell_size = factor(cell_size, levels = c("Small guard cells", "Large guard cells"))
  )

theta <- seq(0, 2 * pi, length.out = 200)
circles <- bind_rows(
  tibble(x = 530 + 28 * cos(theta), y = 0.30 + 0.038 * sin(theta), grp = "small"),
  tibble(x = 530 + 62 * cos(theta), y = 0.14 + 0.083 * sin(theta), grp = "large")
)

pC <- ggplot(df_C, aes(t, gsw, color = cell_size)) +
  geom_line(linewidth = 1) +
  geom_polygon(data = circles, aes(x, y, group = grp, fill = grp),
               color = "gray30", linewidth = 0.5, alpha = 0.2, inherit.aes = FALSE) +
  scale_color_manual(
    values = c("Small guard cells" = col_fast, "Large guard cells" = col_slow),
    name = NULL
  ) +
  scale_fill_manual(values = c("small" = col_fast, "large" = col_slow), guide = "none") +
  annotate("text", x = 530, y = 0.30 + 0.038 + 0.025,
           label = "Small (high SA:V)", color = col_fast, size = 2.7, hjust = 0.5) +
  annotate("text", x = 530, y = 0.14 - 0.083 - 0.022,
           label = "Large (low SA:V)", color = col_slow, size = 2.7, hjust = 0.5) +
  annotate("text", x = 390, y = 0.225,
           label = "SA:V %prop% 1/radius",
           parse = TRUE, size = 3, hjust = 0.5, color = "gray30") +
  annotate("text", x = 390, y = 0.185,
           label = "faster solute flux", size = 2.9, hjust = 0.5, color = "gray30") +
  annotate("text", x = 390, y = 0.150,
           label = "faster turgor change", size = 2.9, hjust = 0.5, color = "gray30") +
  annotate("segment", x = 325, xend = 325, y = 0.205, yend = 0.148,
           arrow = arrow(length = unit(0.06, "in"), type = "open"),
           color = "gray40", linewidth = 0.6) +
  geom_hline(yintercept = gi0, linetype = "dotted", color = "gray55") +
  geom_hline(yintercept = gf0, linetype = "dotted", color = "gray55") +
  scale_y_continuous(limits = c(-0.01, 0.42)) +
  scale_x_continuous(limits = c(0, 700)) +
  labs(x = "Time (s)", y = expression(g[sw]~(mol~m^{-2}~s^{-1}))) +
  base_theme +
  theme(
    legend.position   = c(0.32, 0.84),
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

dP_t <- 0.5
tang <- bind_rows(
  tibble(P = P_lo + c(-1, 1) * dP_t,
         g = pmax(0.02, fgmax_lo + sl_lo * c(-1, 1) * dP_t), grp = "low"),
  tibble(P = P_hi + c(-1, 1) * dP_t,
         g = pmin(1.02, fgmax_hi + sl_hi * c(-1, 1) * dP_t), grp = "high")
)

dP_arrow <- 0.35
arr_df <- bind_rows(
  tibble(x = P_lo - dP_arrow / 2, xend = P_lo + dP_arrow / 2,
         y = fgmax_lo - 0.07, yend = fgmax_lo - 0.07, grp = "low"),
  tibble(x = P_hi - dP_arrow / 2, xend = P_hi + dP_arrow / 2,
         y = fgmax_hi - 0.07, yend = fgmax_hi - 0.07, grp = "high")
)

pD <- ggplot() +
  geom_line(data = df_D, aes(P, g), linewidth = 1.2) +
  geom_line(data = filter(tang, grp == "low"),  aes(P, g),
            color = col_fast, linewidth = 1.2, linetype = "dashed") +
  geom_line(data = filter(tang, grp == "high"), aes(P, g),
            color = col_slow, linewidth = 1.2, linetype = "dashed") +
  geom_segment(data = arr_df,
               aes(x = x, xend = xend, y = y, yend = yend, color = grp),
               arrow = arrow(ends = "both", length = unit(0.07, "in"), type = "closed"),
               linewidth = 0.8, show.legend = FALSE) +
  scale_color_manual(values = c("low" = col_fast, "high" = col_slow)) +
  geom_point(aes(x = P_lo, y = fgmax_lo), color = col_fast, size = 3) +
  geom_point(aes(x = P_hi, y = fgmax_hi), color = col_slow, size = 3) +
  annotate("text", x = P_lo + dP_arrow / 2 + 0.08, y = fgmax_lo - 0.07,
           label = "same~Delta*P[g]", parse = TRUE, color = col_fast, size = 2.7, hjust = 0) +
  annotate("text", x = P_hi + dP_arrow / 2 + 0.08, y = fgmax_hi - 0.07,
           label = "same~Delta*P[g]", parse = TRUE, color = col_slow, size = 2.7, hjust = 0) +
  annotate("text", x = 0.75, y = 0.13,
           label = "Low~f[gmax]~(steep~slope)", parse = TRUE,
           color = col_fast, size = 2.9, hjust = 0) +
  annotate("text", x = 0.75, y = 0.07,
           label = "phantom(0) %->% fast~closure", parse = TRUE,
           color = col_fast, size = 2.9, hjust = 0) +
  annotate("segment", x = 0.73, xend = P_lo + 0.03, y = 0.17, yend = fgmax_lo - 0.01,
           color = col_fast, linewidth = 0.5,
           arrow = arrow(length = unit(0.06, "in"), type = "open")) +
  annotate("text", x = 1.9, y = 0.50,
           label = "High~f[gmax]~(shallow~slope)", parse = TRUE,
           color = col_slow, size = 2.9, hjust = 0) +
  annotate("text", x = 1.9, y = 0.44,
           label = "phantom(0) %->% slow~closure", parse = TRUE,
           color = col_slow, size = 2.9, hjust = 0) +
  annotate("segment", x = 1.88, xend = P_hi + 0.06, y = 0.54, yend = fgmax_hi - 0.01,
           color = col_slow, linewidth = 0.5,
           arrow = arrow(length = unit(0.06, "in"), type = "open")) +
  geom_hline(yintercept = 1, linetype = "dotted", color = "gray55") +
  annotate("text", x = 3.4, y = 1.05, label = expression(g[max]), size = 3.5) +
  labs(
    x = expression(Guard~cell~turgor~pressure~(P[g])),
    y = expression(g[sw]/g[max])
  ) +
  coord_cartesian(ylim = c(0, 1.12), xlim = c(0, 3.6)) +
  base_theme

# ---- Combine and save ----
(pA | pB) / (pC | pD) +
  plot_annotation(tag_levels = "A")

ggsave("figures/conceptual.pdf", width = 180, height = 160, units = "mm")
