rm(list = ls())

# === Load Packages and Project Code ===
devtools::load_all("C:/Users/og7051ak/Nextcloud/Codes for Papers/Tests with Unknown Clusters/Tools/PanelKmeansInference")
devtools::load_all("C:/Users/og7051ak/Nextcloud/Codes for Papers/Tests with Unknown Clusters/Tools/clusteredEPA")

library(dplyr)
library(here)
library(openxlsx)
library(mclust)
library(ggplot2)
library(tidyr)
library(parallel)
library(pbapply)

# === Simulation Parameters ===
N_set <- c(80, 120, 160)
T_eval_set <- c(20, 50, 100, 200)
d_vals <- c(0, 0.125, 0.25, 0.375, 0.5)
n_sim <- 1000
burn <- 10
lrv <- "EWC"
lrv_par <- NULL

set.seed(1)

results_clustered <- list()
results_misc <- list()
results_pcomb <- list()

for (d in d_vals) {
  if (d == 0) {
    N_loop <- N_set
    T_loop <- T_eval_set
    delta_list <- list(`overall_holds` = c(0, 0, 0))
  } else {
    N_loop <- 80
    T_loop <- c(50, 200)
    delta_list <- list(
      `overall_fails` = d / 2 + d * c(-1.2, -0.8, 1),
      `overall_holds` = d * c(-1.2, -0.8, 1)
    )
  }
  
  for (variant in names(delta_list)) {
    for (N in N_loop) {
      for (T_eval in T_loop) {
        delta_vec <- delta_list[[variant]]
        gamma <- rep(NA, N)
        cut1 <- floor(N / 4)
        cut2 <- floor(N / 2)
        gamma[1:cut1] <- 1
        gamma[(cut1 + 1):cut2] <- 2
        gamma[(cut2 + 1):N] <- 3
        id_vec <- rep(1:N, each = T_eval)
        time_vec <- rep(1:T_eval, times = N)
        
        cat(sprintf("Running N = %d, T = %d, d = %.3f, variant = %s\n", N, T_eval, d, variant))
        n_cores <- detectCores(logical = FALSE)
        cl <- makeCluster(n_cores)
        clusterSetRNGStream(cl, iseed = 20250708)
        clusterExport(cl, varlist = ls(), envir = environment())
        clusterEvalQ(cl, {
          devtools::load_all("C:/Users/og7051ak/Nextcloud/Codes for Papers/Tests with Unknown Clusters/Tools/PanelKmeansInference")
          devtools::load_all("C:/Users/og7051ak/Nextcloud/Codes for Papers/Tests with Unknown Clusters/Tools/clusteredEPA")
          library(mclust)
        })
        
        sim_fun <- function(sim) {
          sim_data <- generate_forecast_simulation_data(
            N, Tobs = T_eval + 2, burn,
            alpha = 1, rho_vec = c(0.1, 0.2, 0.3), gamma = gamma,
            phi = 0.2, lambda = 0.2, delta_vec = delta_vec
          )
          
          Y <- sim_data$Y
          forecast_1 <- sim_data$forecast_1
          forecast_2 <- sim_data$forecast_2
          
          L1 <- t((Y[, 2:(T_eval + 1)] - forecast_1[, 2:(T_eval + 1)])^2)
          L2 <- t((Y[, 2:(T_eval + 1)] - forecast_2[, 2:(T_eval + 1)])^2)
          DL <- matrix(L1 - L2, ncol = 1)
          cond_var <- t(Y[, 1:T_eval])
          X <- matrix(cond_var, ncol = 1)
          
          df <- matrix_to_panel_data(cbind(DL, X), Z_names = c("DL", "X"), id_vec, time_vec)
          
          out_all <- list(
            known = clusteredEPA(df, "DL", NULL, "id", "time", test = "epa_clustered_known", gamma = gamma, lrv = lrv, lrv_par = lrv_par),
            known_cond = clusteredEPA(df, "DL", "X", "id", "time", test = "epa_clustered_known", gamma = gamma, lrv = lrv, lrv_par = lrv_par),
            naive = clusteredEPA(df, "DL", NULL, "id", "time", test = "epa_clustered_known", gamma = panel_kmeans_estimation(DL, id_vec, time_vec, Kmax = 5, n_cores = 1)$final_cluster, lrv = lrv, lrv_par = lrv_par),
            naive_cond = clusteredEPA(df, "DL", "X", "id", "time", test = "epa_clustered_known", gamma = panel_kmeans_estimation(cbind(DL, DL * X), id_vec, time_vec, Kmax = 5, n_cores = 1)$final_cluster, lrv = lrv, lrv_par = lrv_par),
            split = clusteredEPA(df, "DL", NULL, "id", "time", Kmax = 5, test = "epa_clustered_split", prop = 0.2, lrv = lrv, lrv_par = lrv_par),
            split_cond = clusteredEPA(df, "DL", "X", "id", "time", Kmax = 5, test = "epa_clustered_split", prop = 0.2, lrv = lrv, lrv_par = lrv_par),
            selective = clusteredEPA(df, "DL", NULL, "id", "time", Kmax = 5, test = "epa_clustered_selective", pcombine_fun = "Genmean_neq", lrv = lrv, lrv_par = lrv_par),
            selective_cond = clusteredEPA(df, "DL", "X", "id", "time", Kmax = 5, test = "epa_clustered_selective", pcombine_fun = "Genmean_neq", lrv = lrv, lrv_par = lrv_par),
            selective_mode2 = clusteredEPA(df, "DL", NULL, "id", "time", Kmax = 5, test = "epa_clustered_selective_mode2", pcombine_fun = "Genmean_neq", lrv = lrv, lrv_par = lrv_par),
            selective_mode2_cond = clusteredEPA(df, "DL", "X", "id", "time", Kmax = 5, test = "epa_clustered_selective_mode2", pcombine_fun = "Genmean_neq", lrv = lrv, lrv_par = lrv_par)
          )
          
          rejection_rates <- c(sapply(out_all, function(x) as.numeric(x$pval <= 0.05)),
                               as.numeric(out_all[["selective"]]$hom_test_pval <= 0.05),
                               as.numeric(out_all[["selective_cond"]]$hom_test_pval <= 0.05),
                               as.numeric(out_all[["selective"]]$overall_pval <= 0.05),
                               as.numeric(out_all[["selective_cond"]]$overall_pval <= 0.05))
          names(rejection_rates)[(length(rejection_rates)-3):length(rejection_rates)] <- c("hom", "hom_cond", "overall", "overall_cond")
          
          out_df <- data.frame(N = N, Tobs = T_eval, d = d, variant = variant, as.list(rejection_rates))
          return(out_df)
        }
        
        sim_results <- pblapply(1:n_sim, sim_fun, cl = cl)
        stopCluster(cl)
        
        combined_df <- bind_rows(sim_results)
        numeric_means <- combined_df %>% select(where(is.numeric)) %>% summarise(across(everything(), ~ mean(.x, na.rm = TRUE)))
        variant <- combined_df[1, "variant"]
        res_df <- cbind(variant, numeric_means)
        
        results_clustered[[length(results_clustered) + 1]] <- res_df %>% select(N, Tobs, d, variant, known, known_cond, naive, naive_cond, split, split_cond, selective, selective_cond, selective_mode2, selective_mode2_cond)
        
        write.xlsx(bind_rows(results_clustered), here("Tools", "clusteredEPA", "Results", "results_clustered.xlsx"))
      }
    }
  }
}