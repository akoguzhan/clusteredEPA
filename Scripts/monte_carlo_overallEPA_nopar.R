rm(list = ls())

devtools::load_all("C:/Users/og7051ak/Nextcloud/Codes for Papers/Tests with Unknown Clusters/Tools/PanelKmeansInference")
devtools::load_all("C:/Users/og7051ak/Nextcloud/Codes for Papers/Tests with Unknown Clusters/Tools/clusteredEPA")

library(dplyr)
library(here)
library(parallel)

# === Simulation parameters ===
N_set <- c(80, 120, 160, 200)
T_eval_set <- c(20, 50, 100, 150)
n_sim <- 2000
burn <- 10
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
    
    n_cores <- max(1, detectCores(logical = FALSE) - 1)
    cl <- makeCluster(n_cores)
    clusterEvalQ(cl, {
      devtools::load_all("C:/Users/og7051ak/Nextcloud/Codes for Papers/Tests with Unknown Clusters/Tools/PanelKmeansInference", quiet = TRUE)
      devtools::load_all("C:/Users/og7051ak/Nextcloud/Codes for Papers/Tests with Unknown Clusters/Tools/clusteredEPA", quiet = TRUE)
      library(mclust)
    })
    clusterExport(cl, varlist = c(
      "generate_forecast_simulation_data", "matrix_to_panel_data", "panel_kmeans_estimation",
      "clusteredEPA", "N", "T_eval", "burn", "gamma", "id_vec", "time_vec", "lrv", "lrv_par"
    ), envir = environment())
    
    sim_fun <- function(sim) {
      sim_data <- generate_forecast_simulation_data(
        N, Tobs = T_eval + 2,
        burn,
        alpha = 1,
        rho_vec = c(0.1, 0.2, 0.3),
        gamma,
        phi = 0.2,
        lambda = 0.2,
        delta_vec = 0 * c(0.1, 0.2, 0.3)
      )
      Y <- sim_data$Y
      forecast_1 <- sim_data$forecast_1
      forecast_2 <- sim_data$forecast_2
      
      L1 <- (Y[, 2:(T_eval+1)] - forecast_1[, 2:(T_eval+1)])^2
      L2 <- (Y[, 2:(T_eval+1)] - forecast_2[, 2:(T_eval+1)])^2
      DL <- matrix(as.vector(L1 - L2), ncol = 1)
      cond_var <- as.vector(Y[, 1:T_eval])
      cond_mat <- matrix(cond_var, ncol = 1)
      X <- DL * cond_mat
      
      df <- matrix_to_panel_data(cbind(DL, cond_mat), Z_names = c("DL", "X"), id_vec, time_vec)
      
      # === Overall EPA Test ===
      res_epa <- try(clusteredEPA(df, "DL", H = NULL, "id", "time", test = "overall_EPA_test", lrv = lrv, lrv_par = lrv_par), silent = TRUE)
      rej_epa <- if (!inherits(res_epa, "try-error") && !is.na(res_epa$p_oepa) && res_epa$p_oepa <= 0.05) 1 else 0
      
      # === Clustered EPA Test ===
      res_clus <- try(clusteredEPA(df, "DL", H = NULL, "id", "time", test = "epa_clustered_known", gamma = gamma, lrv = lrv, lrv_par = lrv_par), silent = TRUE)
      rej_clus <- if (!inherits(res_clus, "try-error") && !is.na(res_clus$pval) && res_clus$pval <= 0.05) 1 else 0
      
      # === Conditional Overall EPA Test ===
      res_epa_cond <- try(clusteredEPA(df, "DL", H = NULL, "id", "time", test = "overall_EPA_test", lrv = lrv, lrv_par = lrv_par), silent = TRUE)
      rej_epa_cond <- if (!inherits(res_epa_cond, "try-error") && !is.na(res_epa_cond$p_oepa) && res_epa_cond$p_oepa <= 0.05) 1 else 0
      
      # === Conditional Clustered EPA Test ===
      res_clus_cond <- try(clusteredEPA(df, "DL", H = "X", "id", "time", test = "epa_clustered_known", gamma = gamma, lrv = lrv, lrv_par = lrv_par), silent = TRUE)
      rej_clus_cond <- if (!inherits(res_clus_cond, "try-error") && !is.na(res_clus_cond$pval) && res_clus_cond$pval <= 0.05) 1 else 0
      
      # === Naive Test: Estimate clusters, then EPA with known clusters ===
      est_gamma <- try(panel_kmeans_estimation(DL, id_vec, time_vec, K = 3, n_cores = 1)$final_cluster, silent = TRUE)
      res_naive <- try(clusteredEPA(df, "DL", H = NULL, "id", "time", test = "epa_clustered_known", gamma = est_gamma, lrv = lrv, lrv_par = lrv_par), silent = TRUE)
      rej_naive <- if (!inherits(res_naive, "try-error") && !is.na(res_naive$pval) && res_naive$pval <= 0.05) 1 else 0
      
      # === Conditional Naive Test ===
      est_gamma_cond <- try(panel_kmeans_estimation(Z = cbind(DL, DL * cond_mat), id_vec, time_vec, K = 3, n_cores = 1)$final_cluster, silent = TRUE)
      res_naive_cond <- try(clusteredEPA(df, "DL", H = "X", "id", "time", test = "epa_clustered_known", gamma = est_gamma_cond, lrv = lrv, lrv_par = lrv_par), silent = TRUE)
      rej_naive_cond <- if (!inherits(res_naive_cond, "try-error") && !is.na(res_naive_cond$pval) && res_naive_cond$pval <= 0.05) 1 else 0
      
      # === Split Sample Test ===
      res_split <- try(clusteredEPA(df, "DL", H = NULL, "id", "time", K = 3, test = "epa_clustered_split", prop = 0.2, lrv = lrv, lrv_par = lrv_par), silent = TRUE)
      rej_split <- if (!inherits(res_split, "try-error") && !is.na(res_split$pval) && res_split$pval <= 0.05) 1 else 0
      
      # === Conditional Split Sample Test ===
      res_split_cond <- try(clusteredEPA(df, "DL", H = "X", "id", "time", K = 3, test = "epa_clustered_split", prop = 0.2, lrv = lrv, lrv_par = lrv_par), silent = TRUE)
      rej_split_cond <- if (!inherits(res_split_cond, "try-error") && !is.na(res_split_cond$pval) && res_split_cond$pval <= 0.05) 1 else 0
      
      # === Selective Inference Test ===
      res_sel <- try(clusteredEPA(df, "DL", H = NULL, "id", "time", K = 3, test = "epa_clustered_selective", pcombine_fun = "pCauchy", lrv = lrv, lrv_par = lrv_par), silent = TRUE)
      rej_sel <- if (!inherits(res_sel, "try-error") && !is.na(res_sel$selective_pval) && res_sel$selective_pval <= 0.05) 1 else 0
      
      # === Conditional Selective Inference Test ===
      res_sel_cond <- try(clusteredEPA(df, "DL", H = "X", "id", "time", K = 3, test = "epa_clustered_selective", pcombine_fun = "pCauchy", lrv = lrv, lrv_par = lrv_par), silent = TRUE)
      rej_sel_cond <- if (!inherits(res_sel_cond, "try-error") && !is.na(res_sel_cond$selective_pval) && res_sel_cond$selective_pval <= 0.05) 1 else 0
      
      c(rej_epa, rej_clus, rej_epa_cond, rej_clus_cond,
        rej_naive, rej_naive_cond,
        rej_split, rej_split_cond,
        rej_sel, rej_sel_cond)
    }
    
    res_list <- parLapply(cl, 1:n_sim, sim_fun)
    stopCluster(cl)
    
    res_mat <- do.call(rbind, res_list)
    
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

# Gather all results into a tidy data frame
results_df <- data.frame(
  N = rep(N_set, each = length(T_eval_set)),
  T = rep(T_eval_set, times = length(N_set)),
  overall_EPA = as.vector(results_overall),
  clustered_EPA = as.vector(results_clustered),
  overall_EPA_conditional = as.vector(results_overall_conditional),
  clustered_EPA_conditional = as.vector(results_clustered_conditional),
  naive = as.vector(results_naive),
  naive_conditional = as.vector(results_naive_conditional),
  split = as.vector(results_split),
  split_conditional = as.vector(results_split_conditional),
  selective = as.vector(results_selective),
  selective_conditional = as.vector(results_selective_conditional)
)

print(results_df)