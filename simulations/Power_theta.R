###simulation to quantify the power to detect ascertainment

Sim_power <- function(R, N_base, M, Nr, K = 0.5, p_base= NULL) {
  
  #────────────────────────────────────────────
  # 1. Input validation
  #────────────────────────────────────────────
  stopifnot(
    length(R) == 1, is.numeric(R), is.finite(R), R > -1, R < 1,
    length(N_base) == 1, is.numeric(N_base), is.finite(N_base), N_base > 100,
    length(M) == 1, is.numeric(M), is.finite(M), M > 10,
    length(Nr) == 1, is.numeric(Nr), is.finite(Nr), Nr > 10, Nr < N_base,
    length(K) == 1, is.numeric(K), is.finite(K), K > 0, K < 1
  )
  
  #────────────────────────────────────────────
  # 2. Allele frequency generation / input
  #────────────────────────────────────────────
  if (is.null(p_base)) {
    p_base <- stats::runif(M, 0.05, 0.95)
  } else {
    stopifnot(
      is.numeric(p_base),
      length(p_base) == M,
      all(p_base > 0),
      all(p_base < 1)
    )
  }
  
  #────────────────────────────────────────────
  # 3. Simulate base population genotypes
  #────────────────────────────────────────────
  X_base <- do.call("cbind", lapply(1:M, function(j) {
    rbinom(N_base, 2, prob = p_base[j])
  }))
  
  #────────────────────────────────────────────
  # 4. Construct reference population (subset of base)
  #────────────────────────────────────────────
  idx <- sample.int(nrow(X_base), Nr)
  Xr <- X_base[idx, , drop = FALSE]
  
  #────────────────────────────────────────────
  # 5. Generate trait-increasing alllele effect sizes and polygenic score
  #────────────────────────────────────────────
  b <- abs(rnorm(M))
  g <- scale(c(X_base %*% b))
  
  #────────────────────────────────────────────
  # 6. Liability model with ascertainment
  #────────────────────────────────────────────
  l <- rnorm(n = N_base, mean = R * g, sd = sqrt(1 - R^2))
  
  # selection threshold based on K
  t <- qnorm(1 - K)
  kept <- which(l > t)
  Ns <- length(kept)
  
  #────────────────────────────────────────────
  # 7. Allele frequencies: selected vs reference
  #────────────────────────────────────────────
  ps <- colMeans(X_base[kept, , drop = FALSE]) / 2
  pr <- colMeans(Xr) / 2
  
  #────────────────────────────────────────────
  # 8. Theta1 estimation (score-based)
  #────────────────────────────────────────────
  u <- sum(2 * b * (ps - pr))
  v <- sum(b^2 * 2 * pr * (1 - pr))
  
  theta1 <- u / sqrt(v)
  
  #────────────────────────────────────────────
  # 9. Variance of theta1 (delta-method expansion)
  #────────────────────────────────────────────
  var_ps <- ps * (1 - ps) / (2 * Ns)
  var_pr <- pr * (1 - pr) / (2 * Nr)
  
  var_u <- sum(4 * b^2 * (var_ps + var_pr))
  var_v <- sum(4 * b^4 * (var_pr * (1 - 2 * pr)^2 + 2 * var_pr^2))
  cov_uv <- sum(-4 * b^3 * (1 - 2 * pr) * var_pr)
  
  var_theta1 <- (1 / v) * (var_u + (u^2 * var_v) / (4 * v^2) - (u * cov_uv / v))
  sd_theta1 <- sqrt(var_theta1)
  
  # Wald-type test
  ts1 <- (theta1 / sd_theta1)^2
  pval1 <- pchisq(ts1, df = 1, lower.tail = FALSE)
  
  #────────────────────────────────────────────
  # 10. Theta2 estimation (weighted regression)
  #────────────────────────────────────────────
  rhs <- b * pr * (1 - pr) / sqrt(v)
  lhs <- ps - pr
  iv_weight <- 1 / var_pr
  
  fit <- lm(lhs ~ rhs, weights = iv_weight)
  coef_tab <- summary(fit)$coefficients
  
  theta2 <- coef_tab[2, 1]
  pval2 <- coef_tab[2, 4]
  
  #────────────────────────────────────────────
  # 11. Output
  #────────────────────────────────────────────
  return(c(pval_theta1 = pval1, pval_theta2 = pval2))
}

#────────────────────────────────────────────
# 1. Simulation settings
#────────────────────────────────────────────
K      <- 0.5   # proportion selected into the study
R      <- 0.4   # ascertainment strength
N_base <- 2000  # base population sample size
M      <- 1000  # number of SNPs
Nr     <- 100   # reference sample size
alpha  <- 0.05  # significance threshold
n_rep  <- 1000    # number of simulation replicates

#────────────────────────────────────────────
# 2. Theoretical expectation
#────────────────────────────────────────────
t_cut  <- qnorm(1 - K)   # liability threshold
i_sel  <- dnorm(t_cut) / K
theta  <- i_sel * R

#────────────────────────────────────────────
# 3. Run simulations
#────────────────────────────────────────────
runs <- t(replicate(
  n = n_rep,
  expr = sim_power(R = R, N_base = N_base, M = M, Nr = Nr, K = K)
))

#────────────────────────────────────────────
# 4. Estimate empirical power
#────────────────────────────────────────────
power_theta1 <- mean(runs[, "pval_theta1"] < alpha)
power_theta2 <- mean(runs[, "pval_theta2"] < alpha)

#────────────────────────────────────────────
# 5. Summary output
#────────────────────────────────────────────
res <- c(
  theta        = theta,
  R            = R,
  power_theta1 = power_theta1,
  power_theta2 = power_theta2
)

print(res)