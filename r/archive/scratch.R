# modeling ideas (move elsewhere, material below is good)
# based on eq. A2 in Buckley et al. (2023)
# X is osmolality gradient (delta pi in buckley paper)
# X_f is target osmolality
# dX / dt = a (X - X_f)  # buckley paper, results in Vico eqn I think
#
# Now suppose that a is a function of X
# a = a0 + a1 * X
# dX / dt = (a0 + a1 * X) * (X - X_f)
# this only has implicit soln according to ChatGpt.
# alternate forms usually result in sigmoidal
# 
# what if we go back to dX / dt = a (X - X_f) but assume that d X_f / dt = b (X_f - X_inf)
# we get this according to chatGPT:
solve_X_dynamics <- function(t, X0, Xf0, X_inf, a, b) {
  if (a == b) {
    # Special case: a == b
    term1 <- X_inf
    term2 <- a * (Xf0 - X_inf) * t * exp(-a * t)
    term3 <- (X0 - X_inf) * exp(a * t)
    return(term1 + term2 + term3)
  } else {
    # General case: a != b
    term1 <- X_inf
    term2 <- (a * (Xf0 - X_inf)) / (b - a) * exp((b - 2 * a) * t)
    term3 <- (X0 - X_inf + (a * (Xf0 - X_inf)) / (b - a)) * exp(a * t)
    return(term1 - term2 + term3)
  }
}

# Not sure if this is useful
nls(gsw ~ gf + dg * exp(-(t_sec / tau)),
    data = df,
    start = list(gf = 0.3, dg = 0.2, tau = 60)
)

fit1 = read_rds(paste0(sk_dir, "LA0107-C_pseudohypo_2000.rds")) 
plot(fit1$data$t_sec, fit1$data$gsw)
ce = conditional_effects(fit1) 

ggplot(ce[[1]], aes(t_sec)) +
  geom_lineribbon(mapping = aes(y = estimate__, ymin = lower__, ymax = upper__),
                  color = "steelblue", fill = "grey") +
  geom_point(data = fit1$data, mapping = aes(y = gsw))


prior1 <- prior(normal(0.2251668, 0.01), nlpar = "gf") +
  prior(normal(-0.1514164, 0.02), nlpar = "dg", ub = 0) +
  prior(normal(600, 600), nlpar = "tau", lb = 0)

fit1 = brm(
  bf(
    gsw ~ gf + dg * exp(-(t_sec / tau)^lambda),
    gf ~ 1 + (1 | curve),
    dg ~ 1 + (1 | curve),
    tau ~ 1 + (1 | curve),
    lambda ~ 1 + (1 | curve),
    # autocor = ~ ar(time = t_sec, p = 1, cov = TRUE),
    nl = TRUE
  ),
  data = df1,
  backend = "cmdstanr",
  chains = 1,
  # prior = prior1,
  family = student
)

fit1
conditional_effects(fit1)

df2 = tibble(t_sec = seq(0, max(df1$t_sec, length.out = 1e2)))
p2 = predict(fit1, newdata = df2)
df2 = mutate(df2, gsw = p2[, "Estimate"], lower = p2[, "Q2.5"], upper = p2[, "Q97.5"])

ggplot(df1, aes(t_sec, gsw)) +
  geom_point() +
  geom_line(data = df2, aes(t_sec, gsw), color = "blue") +
  geom_ribbon(data = df2, aes(t_sec, ymin = lower, ymax = upper), alpha = 0.2) +
  labs(x = "Time (s)", y = "gsw")





