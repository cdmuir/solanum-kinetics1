p_to_a2 <- function(P, pars) {
  pars["a_max"] + (pars["y_0"] - pars["a_max"]) * (1 + P / pars["b"]) ^ pars["a"]
}
pt <- function(t, pars) {
  pars["P_f"] + (pars["P_0"] - pars["P_f"]) * exp(-pars["r"] * t)
}

fn <- function(pars, t, y) {
  P <- pt(t, pars[1:3])
  a2 <- p_to_a2(P, pars[4:7])
  sum((y - a2) ^ 2)
}

pars1 <- c(P_0 = 2, P_f = 1, r = 0.005)
pars2 <- c(y_0 = 150, a_max = 950, b = 1, a = -2)

t <- seq(0, 30 * 60, 1) # 30 minutes
P <- pt(t, pars1)
a2 <- p_to_a2(P, pars2)

plot(P, a2)
plot(t, P, type = "l", xlab = "Time (s)", ylab = "P", main = "Function pt")
plot(t, a2, type = "l", xlab = "P", ylab = "a2", main = "Function p_to_a2")

fit1 <- nls(a2 ~ g_f + (g_i - g_f) * exp(-(t/tau)^lambda), start = c(g_f = 750, g_i = 950, tau = 500, lambda = 1))
points(t, predict(fit1))
cor(a2, predict(fit1))
fit1

fit2 <- optim(c(pars1, pars2), fn, t = t, y = a2)
fit2
