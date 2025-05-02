library(horseshoe)
simulate_once = function(n, tau_type, mu_type, sigma) {
  x1 = rnorm(n, 0, 1)
  x2 = rnorm(n, 0, 1)
  x3 = rnorm(n, 0, 1)
  x4 = rbinom(n, 1, 0.5)
  x5 = sample(1:3, n, replace = TRUE)
  g = 5 - 3 * x5
  u = runif(n, 0, 1)
  
  # tau
  if (tau_type == "homogeneous") {
    tau_true = 3
  } else if (tau_type == "heterogeneous") {
    tau_true = 1 + 2 * x2 * x4 # subject to change
  } else {
    stop("tau_type can only be 'homogeneous' or 'heterogeneous'.")
  }
  
  # mu
  if (mu_type == "linear") {
    mu_true = 1 + g + x1 * x3
  } else if (mu_type == "nonlinear") {
    mu_true = -6 + g + 6 * abs(x3 - 1)
  } else {
    stop("mu_type can only be 'linear' or 'nonlinear'.")
  }
  
  s = sd(mu_true)
  
  # pi
  pi = 0.8 * pnorm(3*mu_true/s-x1/2, 0, 1) + 0.05 + u / 10
  
  z = rbinom(n, 1, pi)
  err = rnorm(n, 0, sigma)
  y = mu_true + tau_true * z + err
  X = as.matrix(x1, x2, x3, x4, x5)
  X_full = cbind(X, X * z)
  
  hs_fit = horseshoe(y, X_full, method.tau = "truncatedCauchy", method.sigma = "Jeffreys", nmc = 1000)
  beta_samples = hs_fit$BetaSamples
  
  p = ncol(X)
  theta_samples = beta_samples[(p + 1):(2 * p), ]
  
  tau_draws = X %*% theta_samples
  
  # posterior summary
  tau_hat = rowMeans(tau_draws)
  tau_sd  = apply(tau_draws, 1, sd)
  lower   = apply(tau_draws, 1, quantile, probs = 0.025)
  upper   = apply(tau_draws, 1, quantile, probs = 0.975)
  
  # evaluation metrics
  rmse = sqrt(mean((tau_hat - tau_true)^2))
  coverage = mean(tau_true >= lower & tau_true <= upper)
  avg_len = mean(upper - lower)
  
  return(c(rmse = rmse, coverage = coverage, interval_length = avg_len))
}