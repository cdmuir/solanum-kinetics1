hist(rstable(1e5, alpha = 1, beta = 1, gamma = 10, delta = 0, pm = 2), breaks = 100, main = "Stable Distribution (alpha=2, beta=1)")
integrate(dstable, alpha = 2, beta = 1, gamma = 1, lower = -0, upper = Inf)

library(stabledist)

set.seed(1)

# --- choose parameters ---
lambda <- 0.6     # stability index in (0,1); also the stretched-exponential exponent
gamma  <- 1.0     # scale parameter for the stable distribution (controls overall timescale)
n_sys  <- 200000  # number of systems (large helps reduce Monte Carlo noise)

# --- draw random rates k (one-sided stable) ---
# For alpha<1 and beta=1, the stable is "totally skewed to the right";
# depending on parameterization, you may still get tiny negative numerical draws.
k <- rstable(n_sys, alpha = lambda, beta = 1, gamma = gamma, delta = 0, pm = 1)

# Guard against numerical negatives (rare): truncate or resample
k <- k[k > 0]
if (length(k) < n_sys) {
  k2 <- rstable(n_sys - length(k), alpha=lambda, beta=1, gamma=gamma, delta=0, pm=1)
  k2 <- k2[k2 > 0]
  k <- c(k, k2)
}

# --- time grid ---
t <- exp(seq(log(1e-3), log(10), length.out = 120))

# --- compute mean curve: E[exp(-k t)] ---
# Do this in chunks to avoid huge memory use
mean_y <- numeric(length(t))
chunk  <- 50000
idx <- split(seq_along(k), ceiling(seq_along(k)/chunk))

for (j in seq_along(t)) {
  s <- 0
  for (ii in idx) {
    s <- s + sum(exp(-k[ii] * t[j]))
  }
  mean_y[j] <- s / length(k)
}

# --- estimate lambda from linearization: log(-log(mean_y)) vs log(t) ---
# Avoid the ends where mean_y is too close to 1 or 0 (numerically unstable)
keep <- mean_y < 0.98 & mean_y > 1e-6
fit <- lm(log(-log(mean_y[keep])) ~ log(t[keep]))

lambda_hat <- unname(coef(fit)[2])
A_hat      <- unname(exp(coef(fit)[1]))  # so that -log(mean_y) ~ A_hat * t^lambda_hat

# Corresponding stretched exponential approximation:
# mean_y ≈ exp(-A_hat * t^lambda_hat) = exp(-(c_hat * t)^lambda_hat), where c_hat = A_hat^(1/lambda_hat)
c_hat <- A_hat^(1/lambda_hat)
y_stretch <- exp(-A_hat * t^lambda_hat)

cat("True lambda:", lambda, "\n")
cat("Estimated lambda:", signif(lambda_hat, 4), "\n")
cat("Estimated A:", signif(A_hat, 4), "\n")
cat("Estimated c:", signif(c_hat, 4), "\n")

# --- plot: mean curve vs fitted stretched exponential ---
plot(t, mean_y, log = "x", type = "l", lwd = 2,
     xlab = "t", ylab = "E[exp(-k t)]",
     main = "Mixture of exponentials with stable rates -> stretched exponential")
lines(t, y_stretch, lwd = 2, lty = 2)
legend("topright",
       legend = c("Monte Carlo mean", "Fitted stretched exp"),
       lty = c(1,2), lwd = 2, bty = "n")

# --- diagnostic plot: should be ~ straight line with slope lambda ---
plot(log(t[keep]), log(-log(mean_y[keep])),
     pch = 16, cex = 0.6,
     xlab = "log t", ylab = "log(-log(mean_y))",
     main = "Linearization: slope ≈ lambda")
abline(fit, lwd = 2)
mtext(sprintf("slope = %.3f (true %.3f)", lambda_hat, lambda), line = 0.5)
