rm(list = ls())

library(dplyr)
library(here)
library(foreach)
library(doSNOW)

source(here("Tools", "clusteredEPA", "R", "test_functions.R"))
source(here("Tools", "clusteredEPA", "R", "tools.R"))

# === Simulation parameters ===
N_set <- c(80, 120, 160, 200)
T_eval_set <- c(20, 50, 100, 150)
n_sim <- 2000
burn <- 50
lrv <- "EWC"
lrv_par <- NULL

set.seed(12)

# === Setup parallel backend with progress bar ===
n_cores <- parallel::detectCores(logical = FALSE)
cl <- makeCluster(n_cores)
registerDoSNOW(cl)
on.exit(stopCluster(cl)) # Ensure cluster is stopped on exit

pb <- txtProgressBar(max = length(N_set) * length(T_eval_set) * n_sim, style = 3)
counter <- 0
progress <- function(n) {
  counter <<- counter + 1
  setTxtProgressBar(pb, counter)
}
opts <- list(progress = progress)
# Note: Progress bar may not be accurate due to parallel updates.

# === Initialize result arrays ===
results_overall <- matrix(NA, nrow = length(N_set), ncol = length(T_eval_set),
                          dimnames = list(paste0("N=", N_set), paste0("T=", T_eval_set)))
results_clustered <- results_overall
results_overall_conditional <- results_overall
results_clustered_conditional <- results_overall

# === Loop over N and T, parallelize over sim ===
for (n_idx in seq_along(N_set)) {
  N <- N_set[n_idx]
  
  # Define gamma vector: N/4, N/4, 2N/4 cluster sizes
  gamma <- rep(NA, N)
  cut1 <- floor(N / 4)
  cut2 <- floor(N / 2)
  gamma[1:cut1] <- 1
  gamma[(cut1 + 1):cut2] <- 2
  gamma[(cut2 + 1):N] <- 3
  
  for (t_idx in seq_along(T_eval_set)) {
    T_eval <- T_eval_set[t_idx]
    T_adj <- T_eval - 1  # Adjust for lag
    
    cat(sprintf("      Running simulations for N = %d, T = %d\n", N, T_eval))
    
    # Precompute id_vec and time_vec outside the parallel loop
    id_vec <- rep(1:N, each = T_adj)
    time_vec <- rep(1:T_adj, times = N)
    
    res_mat <- foreach(sim = 1:n_sim, .combine = rbind,
                       .packages = c("stats", "base"),
                       .options.snow = opts) %dopar% {
                         
                         sim_data <- generate_forecast_simulation_data(
                           N, Tobs = T_eval + 2,
                           burn,
                           alpha = 1,
                           rho_vec = c(0.1, 0.2, 0.3),
                           gamma,
                           phi = 0.2,
                           lambda = 0.2,
                           delta_vec = 0.1 + c(-0.3, -0.1, 0.2)
                         )
                         
                         Y <- sim_data$Y
                         forecast_1 <- sim_data$forecast_1
                         forecast_2 <- sim_data$forecast_2
                         
                         L1 <- (Y[, 2:T_eval] - forecast_1[, 2:T_eval])^2
                         L2 <- (Y[, 2:T_eval] - forecast_2[, 2:T_eval])^2
                         DL <- matrix(as.vector(L1 - L2), ncol = 1)
                         
                         # Compute cond_mat once per simulation
                         cond_var <- as.vector(Y[, 1:(T_eval - 1)])
                         cond_mat <- matrix(cond_var, ncol = 1)
                         
                         # === Overall EPA Test ===
                         res_epa <- overall_EPA_test(Z = DL, id = id_vec, time = time_vec,
                                                     lrv = lrv, lrv_par = lrv_par)
                         rej_epa <- as.numeric(res_epa$p_oepa <= 0.05)
                         
                         # === Clustered EPA Test ===
                         res_clus <- epa_clustered_known(Z = DL, id = id_vec, time = time_vec, gamma = gamma,
                                                         lrv = lrv, lrv_par = lrv_par)
                         rej_clus <- as.numeric(res_clus$pval <= 0.05)
                         
                         # === Conditional Overall EPA Test ===
                         res_epa_cond <- overall_EPA_test(Z = cbind(DL, DL*cond_mat), id = id_vec, time = time_vec,
                                                          lrv = lrv, lrv_par = lrv_par)
                         rej_epa_cond <- as.numeric(res_epa_cond$p_oepa <= 0.05)
                         
                         # === Conditional Clustered EPA Test ===
                         res_clus_cond <- epa_clustered_known(Z = cbind(DL, DL*cond_mat), id = id_vec, time = time_vec,
                                                              gamma = gamma, lrv = lrv, lrv_par = lrv_par)
                         rej_clus_cond <- as.numeric(res_clus_cond$pval <= 0.05)
                         
                         return(c(rej_epa, rej_clus, rej_epa_cond, rej_clus_cond))
                       }
    
    results_overall[n_idx, t_idx]   <- mean(res_mat[, 1])
    results_clustered[n_idx, t_idx] <- mean(res_mat[, 2])
    results_overall_conditional[n_idx, t_idx]   <- mean(res_mat[, 3])
    results_clustered_conditional[n_idx, t_idx] <- mean(res_mat[, 4])
  }
}

# Cluster is stopped automatically by on.exit()

# === Output rejection rates ===
cat("\nOverall EPA Rejection Rates:\n")
print(round(results_overall, 3))

cat("\nClustered EPA Rejection Rates:\n")
print(round(results_clustered, 3))

cat("\nConditional Overall EPA Rejection Rates:\n")
print(round(results_overall_conditional, 3))

cat("\nConditional Clustered EPA Rejection Rates:\n")
print(round(results_clustered_conditional, 3))