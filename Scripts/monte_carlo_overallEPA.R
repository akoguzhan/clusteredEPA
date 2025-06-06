#rm(list = ls())

#devtools::install("C:/Users/og7051ak/Nextcloud/Codes for Papers/Tests with Unknown Clusters/Tools/PanelKmeansInference")
#devtools::install("C:/Users/og7051ak/Nextcloud/Codes for Papers/Tests with Unknown Clusters/Tools/clusteredEPA")

library(dplyr)
library(here)
library(foreach)
library(doSNOW)

# === Simulation parameters ===
N_set <- c(80, 120, 160, 200)
T_eval_set <- c(20, 50, 100, 150)
N_set <- c(80)
T_eval_set <- c(20)
n_sim <- 100
burn <- 50
lrv <- "EWC"
lrv_par <- NULL

set.seed(1)

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

# === Error log ===
error_log <- list()

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

    cat(sprintf("      Running simulations for N = %d, T = %d\n", N, T_eval))

    # Precompute id_vec and time_vec outside the parallel loop
    id_vec <- rep(1:N, each = T_eval)
    time_vec <- rep(1:T_eval, times = N)

    res <- foreach(sim = 1:n_sim, .combine = rbind,
                       .packages = c("PanelKmeansInference", "clusteredEPA", "stats", "base"),
                       .options.snow = opts) %dopar% {

      sim_error_log <- list()

      sim_data <- tryCatch({
        generate_forecast_simulation_data(
          N, Tobs = T_eval + 2,
          burn,
          alpha = 1,
          rho_vec = c(0.1, 0.2, 0.3),
          gamma,
          phi = 0.2,
          lambda = 0.2,
          delta_vec = c(0.1, 0.2, 0.3)
        )
      }, error = function(e) {
        sim_error_log[["generate_forecast_simulation_data"]] <<- e$message
        return(NULL)
      })

      if (is.null(sim_data)) {
        return(c(rep(NA, 10), list(sim_error_log)))
      }

      Y <- sim_data$Y
      forecast_1 <- sim_data$forecast_1
      forecast_2 <- sim_data$forecast_2

      L1 <- (Y[, 2:(T_eval+1)] - forecast_1[, 2:(T_eval+1)])^2
      L2 <- (Y[, 2:(T_eval+1)] - forecast_2[, 2:(T_eval+1)])^2
      DL <- matrix(as.vector(L1 - L2), ncol = 1)

      cond_var <- as.vector(Y[, 1:T_eval])
      cond_mat <- matrix(cond_var, ncol = 1)

      # === Overall EPA Test ===
      res_unc_oepa_known <- tryCatch({
        overall_EPA_test(Z = DL, id = id_vec, time = time_vec, lrv = lrv, lrv_par = lrv_par)
      }, error = function(e) {
        sim_error_log[["overall_EPA_test"]] <<- e$message
        return(list(p_oepa = NA))
      })
      rej_unc_oepa_known <- as.numeric(!is.na(res_unc_oepa_known$p_oepa) && res_unc_oepa_known$p_oepa <= 0.05)
      
      # === Conditional Overall EPA Test ===
      res_cond_oepa_known <- tryCatch({
        overall_EPA_test(Z = cbind(DL, DL*cond_mat), id = id_vec, time = time_vec, lrv = lrv, lrv_par = lrv_par)
      }, error = function(e) {
        sim_error_log[["overall_EPA_test_cond"]] <<- e$message
        return(list(p_oepa = NA))
      })
      rej_cond_oepa_known <- as.numeric(!is.na(res_cond_oepa_known$p_oepa) && res_cond_oepa_known$p_oepa <= 0.05)

      # === Clustered EPA Test ===
      res_unc_cepa_known <- tryCatch({
        epa_clustered_known(Z = DL, id = id_vec, time = time_vec, gamma = gamma, lrv = lrv, lrv_par = lrv_par)
      }, error = function(e) {
        sim_error_log[["epa_clustered_known"]] <<- e$message
        return(list(pval = NA))
      })
      rej_unc_cepa_known <- as.numeric(!is.na(res_unc_cepa_known$pval) && res_unc_cepa_known$pval <= 0.05)

      # === Conditional Clustered EPA Test ===
      res_cond_cepa_known <- tryCatch({
        epa_clustered_known(Z = cbind(DL, DL*cond_mat), id = id_vec, time = time_vec, gamma = gamma, lrv = lrv, lrv_par = lrv_par)
      }, error = function(e) {
        sim_error_log[["epa_clustered_known_cond"]] <<- e$message
        return(list(pval = NA))
      })
      rej_cond_cepa_known <- as.numeric(!is.na(res_cond_cepa_known$pval) && res_cond_cepa_known$pval <= 0.05)

      # === Naive Test: Estimate clusters, then EPA with known clusters ===
      est_gamma <- tryCatch({
        PanelKmeansInference::panel_kmeans_estimation(DL, id_vec, time_vec, K = 3)$final_cluster
      }, error = function(e) {
        sim_error_log[["panel_kmeans_estimation"]] <<- e$message
        return(rep(NA, N))
      })
      res_unc_cepa_naive <- tryCatch({
        epa_clustered_known(Z = DL, id = id_vec, time = time_vec, gamma = est_gamma, lrv = lrv, lrv_par = lrv_par)
      }, error = function(e) {
        sim_error_log[["epa_clustered_known_naive"]] <<- e$message
        return(list(pval = NA))
      })
      rej_unc_cepa_naive <- as.numeric(!is.na(res_unc_cepa_naive$pval) && res_unc_cepa_naive$pval <= 0.05)

      # === Conditional Naive Test ===
      res_cond_cepa_naive <- tryCatch({
        epa_clustered_known(Z = cbind(DL, DL*cond_mat), id = id_vec, time = time_vec, gamma = est_gamma, lrv = lrv, lrv_par = lrv_par)
      }, error = function(e) {
        sim_error_log[["epa_clustered_known_naive_cond"]] <<- e$message
        return(list(pval = NA))
      })
      rej_cond_cepa_naive <- as.numeric(!is.na(res_cond_cepa_naive$pval) && res_cond_cepa_naive$pval <= 0.05)

      # === Split Sample Test ===
      res_unc_cepa_split <- tryCatch({
        epa_clustered_split(df = data.frame(DL = DL, id = id_vec, time = time_vec),
                            id = "id", time = "time", Z_names = "DL", K = 3, lrv = lrv, lrv_par = lrv_par)
      }, error = function(e) {
        sim_error_log[["epa_clustered_split"]] <<- e$message
        return(list(pval = NA))
      })
      rej_unc_cepa_split <- as.numeric(!is.na(res_unc_cepa_split$pval) && res_unc_cepa_split$pval <= 0.05)

      # === Conditional Split Sample Test ===
      res_split_cond <- tryCatch({
        epa_clustered_split(df = data.frame(DL = DL, cond = DL*cond_mat, id = id_vec, time = time_vec),
                            id = "id", time = "time", Z_names = c("DL", "cond"), K = 3, lrv = lrv, lrv_par = lrv_par)
      }, error = function(e) {
        sim_error_log[["epa_clustered_split_cond"]] <<- e$message
        return(list(pval = NA))
      })
      rej_cond_cepa_split <- as.numeric(!is.na(res_split_cond$pval) && res_split_cond$pval <= 0.05)

      # === Selective Inference Test ===
      res_unc_cepa_sel <- tryCatch({
        epa_clustered_selective(Z = DL, id = id_vec, time = time_vec, K = 3, lrv = lrv, lrv_par = lrv_par)
      }, error = function(e) {
        sim_error_log[["epa_clustered_selective"]] <<- e$message
        return(list(pval = NA))
      })
      rej_unc_cepa_sel <- as.numeric(!is.na(res_unc_cepa_sel$pval) && res_unc_cepa_sel$pval <= 0.05)

      # === Conditional Selective Inference Test ===
      res_sel_cond <- tryCatch({
        epa_clustered_selective(Z = cbind(DL, DL*cond_mat), id = id_vec, time = time_vec, K = 3, lrv = lrv, lrv_par = lrv_par)
      }, error = function(e) {
        sim_error_log[["epa_clustered_selective_cond"]] <<- e$message
        return(list(pval = NA))
      })
      rej_cond_cepa_sel <- as.numeric(!is.na(res_sel_cond$pval) && res_sel_cond$pval <= 0.05)

      c(rej_unc_oepa_known, rej_cond_oepa_known,
        rej_unc_cepa_known, rej_cond_cepa_known,
        rej_unc_cepa_naive, rej_cond_cepa_naive,
        rej_unc_cepa_split, rej_cond_cepa_split,
        rej_unc_cepa_sel,   rej_cond_cepa_sel)
    }

    # Separate results and error logs
    res_mat <- do.call(rbind, lapply(res, function(x) x[1:10]))
    print(res_mat)

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

# Cluster is stopped automatically by on.exit()

# === Output rejection rates ===
cat("\nOverall EPA Rejection Rates:\n")
print(round(results_overall, 3))

cat("\nConditional Overall EPA Rejection Rates:\n")
print(round(results_overall_conditional, 3))

cat("\nClustered EPA Rejection Rates:\n")
print(round(results_clustered, 3))

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

# # Install openxlsx if not already installed
# if (!requireNamespace("openxlsx", quietly = TRUE)) install.packages("openxlsx")
# library(openxlsx)
# 
# # Arrange result matrices into a 5x2 list matrix
# result_names <- c(
#   "results_overall", "results_overall_conditional",
#   "results_clustered", "results_clustered_conditional",
#   "results_naive", "results_naive_conditional",
#   "results_split", "results_split_conditional",
#   "results_selective", "results_selective_conditional"
# )
# result_matrices <- mget(result_names)
# 
# result_grid <- matrix(result_matrices, nrow = 5, ncol = 2, byrow = TRUE,
#                       dimnames = list(
#                         c("overall", "known", "naive", "split", "selective"),
#                         c("unconditional", "conditional")
#                       ))
# 
# # Export each matrix to a separate sheet in an Excel file
# wb <- createWorkbook()
# for (i in seq_along(result_names)) {
#   addWorksheet(wb, result_names[i])
#   writeData(wb, result_names[i], result_matrices[[i]])
# }
# saveWorkbook(wb, file = here("Tools", "clusteredEPA", "Results", "simulation_results.xlsx"), overwrite = TRUE)