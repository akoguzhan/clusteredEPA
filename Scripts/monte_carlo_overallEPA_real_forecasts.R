rm(list = ls())

library(dplyr)
library(here)
library(foreach)
library(doSNOW)

source(here("Tools", "clusteredEPA", "R", "overall_EPA_test.R"))
source(here("Tools", "clusteredEPA", "R", "forecasting.R"))
source(here("Tools", "clusteredEPA", "R", "tools.R"))

# Simulation parameters
T_eval_set <- c(20, 50, 100)
T_forecast_max <- max(T_eval_set)
n_sim <- 1000
N <- 20
T_est <- 20
lrv <- "EWC"
lrv_par <- NULL
burn <- 10

set.seed(1)

# Setup parallel backend with progress bar
n_cores <- parallel::detectCores(logical = FALSE)
cl <- makeCluster(n_cores)
registerDoSNOW(cl)

pb <- txtProgressBar(max = n_sim, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
op <- list(progress = progress)

# Monte Carlo loop
results_matrix <- foreach(sim = 1:n_sim, .combine = rbind, .packages = c("stats", "base", "matrixStats"), .options.snow = op) %dopar% {
  T_all <- burn + T_est + T_forecast_max
  T_total <- T_est + T_forecast_max
  
  # DGP
  F_t <- rnorm(T_all)
  Y <- matrix(0, nrow = N, ncol = T_all)
  U <- matrix(rnorm(N * T_all), N, T_all)
  for (t in 2:T_all) {
    Y[, t] <- 0.5 * Y[, t - 1] + F_t[t] + U[, t]
  }
  
  Y <- Y[, (burn + 1):T_all]
  X1 <- matrix(rnorm(N * T_total), N, T_total)
  X2 <- matrix(rnorm(N * T_total), N, T_total)
  
  # Forecasts
  forecasts <- array(NA, dim = c(N, T_forecast_max, 2))
  actuals <- Y[, (T_total - T_forecast_max + 1):T_total]
  
  for (t in 1:T_forecast_max) {
    train <- t:(T_est + t - 1)
    
    y_train_mat <- Y[, train]
    x1_train_mat <- X1[, train]
    x2_train_mat <- X2[, train]
    
    y_train <- as.vector(t(y_train_mat))
    x1_train <- matrix(as.vector(t(x1_train_mat)), ncol = 1)
    x2_train <- matrix(as.vector(t(x2_train_mat)), ncol = 1)
    id_train <- rep(1:N, each = T_est)
    time_train <- rep(1:T_est, times = N)
    
    forecasts[, t, 1] <- forecast_arx_heterogeneous(y_train, x1_train, id_train, time_train, p = 1, q = 1)
    forecasts[, t, 2] <- forecast_arx_heterogeneous(y_train, x2_train, id_train, time_train, p = 1, q = 1)
  }
  
  # Loss differentials
  L1 <- (actuals - forecasts[, , 1])^2
  L2 <- (actuals - forecasts[, , 2])^2
  DL_full <- L1 - L2
  
  row_result <- numeric(length(T_eval_set))
  for (j in seq_along(T_eval_set)) {
    T_eval <- T_eval_set[j]
    id_vec <- rep(1:N, each = T_eval)
    time_vec <- rep(1:T_eval, times = N)
    
    DL_subset <- matrix(DL_full[, 1:T_eval], ncol = 1)
    colnames(DL_subset) <- "DL"
    res_epa <- overall_EPA_test(Z = DL_subset, id = id_vec, time = time_vec, lrv = lrv, lrv_par = lrv_par)
    row_result[j] <- as.numeric(res_epa$p_oepa <= 0.05)
  }
  
  row_result
}

stopCluster(cl)

# Output rejection rates
rejection_rates <- colMeans(results_matrix)
print(round(rejection_rates, 3))
