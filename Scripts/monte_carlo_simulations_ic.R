rm(list = ls())

# === Load Packages and Project Code ===
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
d_vals <- rev(c(0, 0.125, 0.25, 0.375, 0.5))
n_sim <- 200
burn <- 10
lrv <- "EWC"
lrv_par <- NULL
kappa <- 5
Ninit <- 10

set.seed(1)

results_ic_raw <- list()
results_ic_summary <- list()

for (d in d_vals) {
  if (d == 0) {
    N_loop <- N_set
    T_loop <- T_eval_set
    delta_list <- list(`overall_holds` = c(0, 0, 0))
  } else {
    N_loop <- c(80, 120)
    T_loop <- c(50, 200)
    delta_list <- list(
      `overall_fails` = d/2 + d * c(-1.2, -0.8, 1),
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
          library(here)
          devtools::load_all(here("Tools", "clusteredEPA"))
          devtools::load_all(here("Tools", "PanelKmeansInference"))
          library(mclust)
        })
        
        sim_fun <- function(sim) {
          sim_data <- generate_forecast_simulation_data(
            N, Tobs = T_eval + 2, burn,
            alpha = c(1, 1, 1), rho_vec = c(0.1, 0.2, 0.3), gamma = gamma,
            phi = 0.2, lambda = 0.2, delta_vec = delta_vec
          )
          
          Y <- sim_data$Y
          forecast_1 <- sim_data$forecast_1
          forecast_2 <- sim_data$forecast_2
          
          L1 <- t((Y[, 2:(T_eval + 1)] - forecast_1[, 2:(T_eval + 1)])^2)
          L2 <- t((Y[, 2:(T_eval + 1)] - forecast_2[, 2:(T_eval + 1)])^2)
          DL <- matrix(L1 - L2, ncol = 1)
          X <- matrix(t(Y[, 1:T_eval]), ncol = 1)
          
          est_unc <- panel_kmeans_cv_time(DL, id_vec, time_vec, Kmax = 5, nfolds = 10, Ninit = Ninit)
          cl_hat_unc <- est_unc$final_cluster
          K_hat_unc <- length(unique(cl_hat_unc))
          
          rand_unc <- adjustedRandIndex(cl_hat_unc, gamma)
          recovery_unc <- same_cl(cl_hat_unc, gamma, 3)
          correct_K_unc <- as.numeric(K_hat_unc == 3)
          
          est_cond <- panel_kmeans_cv_time(cbind(DL, DL * X), id_vec, time_vec, Kmax = 5, nfolds = 10, Ninit = Ninit)
          cl_hat_cond <- est_cond$final_cluster
          K_hat_cond <- length(unique(cl_hat_cond))
          
          rand_cond <- adjustedRandIndex(cl_hat_cond, gamma)
          recovery_cond <- same_cl(cl_hat_cond, gamma, 3)
          correct_K_cond <- as.numeric(K_hat_cond == 3)
          
          return(data.frame(N = N, T_eval = T_eval, d = d, variant = variant,
                            K_hat_unc = K_hat_unc, correct_K_unc = correct_K_unc,
                            rand_unc = rand_unc, recovery_unc = recovery_unc,
                            K_hat_cond = K_hat_cond, correct_K_cond = correct_K_cond,
                            rand_cond = rand_cond, recovery_cond = recovery_cond))
        }
        
        sim_results <- pblapply(1:n_sim, sim_fun, cl = cl)
        stopCluster(cl)
        
        df_all <- bind_rows(sim_results)
        results_ic_raw[[length(results_ic_raw) + 1]] <- df_all
        
        df_summary <- df_all %>%
          summarise(across(where(is.numeric), ~ mean(.x, na.rm = TRUE))) %>%
          mutate(N = N, T_eval = T_eval, d = d, variant = variant)
        
        results_ic_summary[[length(results_ic_summary) + 1]] <- df_summary
        
        write.xlsx(bind_rows(results_ic_summary), here("Tools", "clusteredEPA", "Results", "results_ic_summary.xlsx"))
      }
    }
  }
}
