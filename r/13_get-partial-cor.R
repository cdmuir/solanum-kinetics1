fit_amphi = read_rds("objects/best_amphi_model.rds")

# Function from ChatGPT
make_precision_phy <- function(draws_df) {
  # Variable order (match your column names)
  vars <- c("loglambdamean_Intercept",
            "logtaumean_Intercept",
            "loggcl_Intercept",
            "logfgmax_Intercept")
  
  # SD columns
  sd_cols <- paste0("sd_phy__", vars)
  
  # Helper to find the right cor column regardless of ordering in the name
  cor_col <- function(a, b) {
    nm1 <- paste0("cor_phy__", a, "__", b)
    nm2 <- paste0("cor_phy__", b, "__", a)
    if (nm1 %in% names(draws_df)) return(nm1)
    if (nm2 %in% names(draws_df)) return(nm2)
    stop("Missing correlation column for: ", a, " and ", b)
  }
  
  # Extract SDs as matrix: draws x 4
  S <- as.matrix(draws_df[, sd_cols, drop = FALSE])
  Dn <- nrow(S)
  P  <- length(vars)
  
  # Build correlation array R: draws x 4 x 4
  R <- array(0, dim = c(Dn, P, P))
  for (d in seq_len(Dn)) diag(R[d,,]) <- 1
  
  for (i in 1:(P - 1)) {
    for (j in (i + 1):P) {
      cc <- cor_col(vars[i], vars[j])
      R[, i, j] <- draws_df[[cc]]
      R[, j, i] <- draws_df[[cc]]
    }
  }
  
  # Build covariance array Sigma and precision array Omega
  Sigma <- array(NA_real_, dim = c(Dn, P, P))
  Omega <- array(NA_real_, dim = c(Dn, P, P))
  
  for (d in seq_len(Dn)) {
    Dmat <- diag(S[d, ], nrow = P)
    Sigma_d <- Dmat %*% R[d, , ] %*% Dmat
    Sigma[d, , ] <- Sigma_d
    Omega[d, , ] <- solve(Sigma_d)
  }
  
  dimnames(Omega) <- list(NULL, vars, vars)
  # dimnames(Sigma) <- list(NULL, vars, vars)
  Omega
}

# Usage:
# Omega_draws <- make_precision_phy(draws_df)
# Omega_draws[1, , ]   # precision matrix for first draw

draws_df = fit_amphi |> 
  as_draws_df() |>
  select(starts_with("."), starts_with("sd_phy__"), starts_with("cor_phy__")) 

Omega_draws <- make_precision_phy(draws_df)

-cov2cor(Omega_draws[1000,,])
