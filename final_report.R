# data generating process

generate_data = function(n, tau_type, mu_type, sigma) {
  #' @param n Number of observations (rows in the dataset).
  #' @param tau_type Type of tau. Can be either "homogeneous" or "heterogeneous".
  #' @param mu_type Type of mu. Can be either "linear" or "nonlinear".
  #' @param sigma Standard deviation of the error term.
  
  x1 = rnorm(n, 0, 1)
  x2 = rnorm(n, 0, 1)
  x3 = rnorm(n, 0, 1)
  x4 = rbinom(n, 1, 0.5)
  x5 = sample(1:3, n, replace = TRUE)
  g = 5 - 3 * x5
  u = runif(n, 0, 1)
  
  # tau
  if (tau_type == "homogeneous") {
    tau = 3
  } else if (tau_type == "heterogeneous") {
    tau = 1 + 2 * x2 * x5
  } else {
    stop("tau_type can only be 'homogeneous' or 'heterogeneous'.")
  }
  
  # mu
  if (mu_type == "linear") {
    mu = 1 + g + x1 * x3
  } else if (mu_type == "nonlinear") {
    mu = -6 + g + 6 * abs(x3 - 1)
  } else {
    stop("mu_type can only be 'linear' or 'nonlinear'.")
  }
  
  s = sd(mu)
  
  # pi
  pi = 0.8 * pnorm(3*mu/s-x1/2, 0, 1) + 0.05 + u / 10
  
  z = rbinom(n, 1, pi)
  err = rnorm(n, 0, sigma)
  y = mu + tau * z + err
  
  return(data.frame(y, z, x1, x2, x3, x4, x5, tau))
}

sample_sizes <- c(250, 500)
effect_types <- c("homogeneous", "heterogeneous")
mu_types     <- c("nonlinear","linear")         # you can extend this if you want other μ‐settings
seeds        <- 2020:2039

rmse     <- function(true, est)  sqrt(mean((true - est)^2))
coverage <- function(true, L, U) mean(true >= L & true <= U)

# storage
nreps <- length(sample_sizes) * length(effect_types) * length(mu_types) * length(seeds)
out <- matrix(NA, nreps, 7,
              dimnames = list(NULL,
                              c("sample_size", "effect_type", "mu_type", "rmse_ate","cover_ate","len_ate",
                                "rmse_cate")))

n_sim    <- 5
i <- 1
for (ss in sample_sizes) {
  for (et in effect_types) {
    for (mt in mu_types) {
      for (seed in seeds) {
        set.seed(seeds)
dat       <- generate_data(ss, et, mt, 1)
X         <- model.matrix(~ x1 + x2 + x3 + x4 + x5, data = dat)[, -1]
Tt        <- dat$z
Y         <- dat$y
tau_true  <- dat$tau

# ---- fit bartc ----
fit_bc <- bartc(
  response    = Y,
  treatment   = Tt,
  confounders = data.matrix(X),
  n.samples   = 1000,
  n.burn      = 1000,
)

# ---- extract summaries ----

cate_draws <- extract(fit_bc,
                      type          = "icate",
                      sample        = "all",
                      combineChains = TRUE)
# ensure it's a matrix
cate_draws <- as.matrix(cate_draws)
tau_hat <- colMeans(cate_draws)

tau_lower <- apply(cate_draws, 2, quantile, probs = 0.025)
tau_upper <- apply(cate_draws, 2, quantile, probs = 0.975)

# ATE draws & interval
ate_draws <- extract(fit_bc,
                     type          = "sate",
                     sample        = "all",
                     combineChains = TRUE)
ate_hat   <- mean(ate_draws)
ate_lo    <- quantile(ate_draws, 0.025)
ate_hi    <- quantile(ate_draws, 0.975)

# 4) metrics
out[i,"sample_size"] <- ss
out[i,"effect_type"] <- et
out[i,"mu_type"]     <- mt
out[i, "rmse_cate"]  <-     rmse(tau_true, tau_hat)

out[i, "rmse_ate"]   <-     rmse(mean(tau_true), ate_draws)
out[i, "cover_ate"]  <- mean(tau_true >= tau_lower & tau_true <= tau_upper)
out[i, "len_ate"]    <-     ate_hi - ate_lo
i<-i+1
      }
    }
  }
}
str(out) 
out <- as.data.frame(out)# Check the structure of your data

# Convert relevant columns to numeric if they are not
out$rmse_ate <- as.numeric(as.character(out$rmse_ate))
out$cover_ate <- as.numeric(as.character(out$cover_ate))
out$len_ate <- as.numeric(as.character(out$len_ate))
out$rmse_cate <- as.numeric(as.character(out$rmse_cate))

summary_table <- aggregate(
  cbind(rmse_ate,cover_ate,len_ate,rmse_cate) ~ sample_size + effect_type + mu_type,
  data = out, FUN = mean
)
