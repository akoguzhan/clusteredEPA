rm(list = ls())

devtools::load_all("C:/Users/og7051ak/Nextcloud/Codes for Papers/Tests with Unknown Clusters/Tools/PanelKmeansInference")
devtools::load_all("C:/Users/og7051ak/Nextcloud/Codes for Papers/Tests with Unknown Clusters/Tools/clusteredEPA")

library(dplyr)
library(here)

# === Simulation parameters ===
N_set <- c(80, 120, 160, 200)
T_eval_set <- c(20, 50, 100, 150)
n_sim <- 5
burn <- 50
lrv <- "EWC"
lrv_par <- NULL

set.seed(1)

# === Initialize result arrays ===
results_overall <- matrix(NA, nrow = length(N_set), ncol = length(T_eval_set),
                          dimnames = list(paste0("N=", N_set), paste0("T=", T_eval_set)))
results_clustered <- results_overall
results_overall_conditional <- results_overall
results_clustered_conditional <- results_overall

results_naive <- results_overall
results_naive_conditional <- results_overall
results_split <- results_overall
results_split_conditional <- results_overall
results_selective <- results_overall
results_selective_conditional <- results_overall

# === Loop over N and T, serial version ===
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

    cat(sprintf("      Running simulations for N = %d, T = %d\n", N, T_eval))
    
    # Precompute id_vec and time_vec outside the sim loop
    id_vec <- rep(1:N, each = T_eval)
    time_vec <- rep(1:T_eval, times = N)
    
    res_mat <- matrix(NA, nrow = n_sim, ncol = 10)
    
    for (sim in 1:n_sim) {
      sim_data <- generate_forecast_simulation_data(
        N, Tobs = T_eval + 2,
        burn,
        alpha = 1,
        rho_vec = c(0.1, 0.2, 0.3),
        gamma,
        phi = 0.2,
        lambda = 0.2,
        delta_vec = c(0.1, 0.2, 0.3)
      )
      
      Y <- sim_data$Y
      forecast_1 <- sim_data$forecast_1
      forecast_2 <- sim_data$forecast_2
      
      L1 <- (Y[, 2:(T_eval+1)] - forecast_1[, 2:(T_eval+1)])^2
      L2 <- (Y[, 2:(T_eval+1)] - forecast_2[, 2:(T_eval+1)])^2
      DL <- matrix(as.vector(L1 - L2), ncol = 1)

      # Compute cond_mat once per simulation
      cond_var <- as.vector(Y[, 1:T_eval])
      cond_mat <- matrix(cond_var, ncol = 1)

      # === Overall EPA Test ===
      res_epa <- overall_EPA_test(Z = DL, id = id_vec, time = time_vec, lrv = lrv, lrv_par = lrv_par)
      rej_epa <- as.numeric(!is.na(res_epa$p_oepa) && res_epa$p_oepa <= 0.05)

      # === Clustered EPA Test ===
      res_clus <- epa_clustered_known(Z = DL, id = id_vec, time = time_vec, gamma = gamma, lrv = lrv, lrv_par = lrv_par)
      rej_clus <- as.numeric(!is.na(res_clus$pval) && res_clus$pval <= 0.05)

      # === Conditional Overall EPA Test ===
      res_epa_cond <- overall_EPA_test(Z = cbind(DL, DL*cond_mat), id = id_vec, time = time_vec, lrv = lrv, lrv_par = lrv_par)
      rej_epa_cond <- as.numeric(!is.na(res_epa_cond$p_oepa) && res_epa_cond$p_oepa <= 0.05)

      # === Conditional Clustered EPA Test ===
      res_clus_cond <- epa_clustered_known(Z = cbind(DL, DL*cond_mat), id = id_vec, time = time_vec, gamma = gamma, lrv = lrv, lrv_par = lrv_par)
      rej_clus_cond <- as.numeric(!is.na(res_clus_cond$pval) && res_clus_cond$pval <= 0.05)

      # === Naive Test: Estimate clusters, then EPA with known clusters ===
      est_gamma <- panel_kmeans_estimation(DL, id_vec, time_vec, K = 3, n_cores = parallel::detectCores(logical = FALSE) - 1)$final_cluster
      res_naive <- epa_clustered_known(Z = DL, id = id_vec, time = time_vec, gamma = est_gamma, lrv = lrv, lrv_par = lrv_par)
      rej_naive <- as.numeric(!is.na(res_naive$pval) && res_naive$pval <= 0.05)

      # === Conditional Naive Test ===
      est_gamma <- panel_kmeans_estimation(Z = cbind(DL, DL*cond_mat), id_vec, time_vec, K = 3, n_cores = parallel::detectCores(logical = FALSE) - 1)$final_cluster
      res_naive_cond <- epa_clustered_known(Z = cbind(DL, DL*cond_mat), id = id_vec, time = time_vec, gamma = est_gamma, lrv = lrv, lrv_par = lrv_par)
      rej_naive_cond <- as.numeric(!is.na(res_naive_cond$pval) && res_naive_cond$pval <= 0.05)

      # === Split Sample Test ===
      res_split <- epa_clustered_split(df = data.frame(DL = DL, id = id_vec, time = time_vec),
                            id = "id", time = "time", Z_names = "DL", K = 3, lrv = lrv, lrv_par = lrv_par)
      rej_split <- as.numeric(!is.na(res_split$pval) && res_split$pval <= 0.05)

      # === Conditional Split Sample Test ===
      res_split_cond <- epa_clustered_split(df = data.frame(DL = DL, cond = DL*cond_mat, id = id_vec, time = time_vec),
                            id = "id", time = "time", Z_names = c("DL", "cond"), K = 3, lrv = lrv, lrv_par = lrv_par)
      rej_split_cond <- as.numeric(!is.na(res_split_cond$pval) && res_split_cond$pval <= 0.05)

      # === Selective Inference Test ===
      r <- (1/4 + 4)/2
      res_sel <- epa_clustered_selective(Z = DL, id = id_vec, time = time_vec, K = 3, r = r, lrv = lrv, lrv_par = lrv_par)
      rej_sel <- as.numeric(!is.na(res_sel$pval) && res_sel$pval <= 0.05)
      
      # === Conditional Selective Inference Test ===
      r <- 
      res_sel_cond <- epa_clustered_selective(Z = cbind(DL, DL*cond_mat), id = id_vec, time = time_vec, K = 3, r = r, lrv = lrv, lrv_par = lrv_par)
      rej_sel_cond <- as.numeric(!is.na(res_sel_cond$pval) && res_sel_cond$pval <= 0.05)
      
      res_mat[sim, ] <- c(rej_epa, rej_clus, rej_epa_cond, rej_clus_cond,
                          rej_naive, rej_naive_cond,
                          rej_split, rej_split_cond,
                          rej_sel, rej_sel_cond)
    }
    
    results_overall[n_idx, t_idx]   <- mean(res_mat[, 1], na.rm = TRUE)
    results_clustered[n_idx, t_idx] <- mean(res_mat[, 2], na.rm = TRUE)
    results_overall_conditional[n_idx, t_idx]   <- mean(res_mat[, 3], na.rm = TRUE)
    results_clustered_conditional[n_idx, t_idx] <- mean(res_mat[, 4], na.rm = TRUE)
    results_naive[n_idx, t_idx] <- mean(res_mat[, 5], na.rm = TRUE)
    results_naive_conditional[n_idx, t_idx] <- mean(res_mat[, 6], na.rm = TRUE)
    results_split[n_idx, t_idx] <- mean(res_mat[, 7], na.rm = TRUE)
    results_split_conditional[n_idx, t_idx] <- mean(res_mat[, 8], na.rm = TRUE)
    results_selective[n_idx, t_idx] <- mean(res_mat[, 9], na.rm = TRUE)
    results_selective_conditional[n_idx, t_idx] <- mean(res_mat[, 10], na.rm = TRUE)
  }
}

# === Output rejection rates ===
cat("\nOverall EPA Rejection Rates:\n")
print(round(results_overall, 3))

cat("\nClustered EPA Rejection Rates:\n")
print(round(results_clustered, 3))

cat("\nConditional Overall EPA Rejection Rates:\n")
print(round(results_overall_conditional, 3))

cat("\nConditional Clustered EPA Rejection Rates:\n")
print(round(results_clustered_conditional, 3))

cat("\nNaive Clustered EPA Rejection Rates:\n")
print(round(results_naive, 3))

cat("\nConditional Naive Clustered EPA Rejection Rates:\n")
print(round(results_naive_conditional, 3))

cat("\nSplit Sample Clustered EPA Rejection Rates:\n")
print(round(results_split, 3))

cat("\nConditional Split Sample Clustered EPA Rejection Rates:\n")
print(round(results_split_conditional, 3))

cat("\nSelective Inference Clustered EPA Rejection Rates:\n")
print(round(results_selective, 3))

cat("\nConditional Selective Inference Clustered EPA Rejection Rates:\n")
print(round(results_selective_conditional, 3))