rm(list = ls())

library(dplyr)
library(here)
library(foreach)
library(doSNOW)
library(openxlsx)

source(here("Tools", "clusteredEPA", "R", "test_functions.R"))
source(here("Tools", "clusteredEPA", "R", "tools.R"))

# === Simulation Design Grid ===
param_grid <- expand.grid(
  rho_case = c("null", "alt"),
  delta_case = c("null", "alt"),
  phi_case = c("null", "alt"),
  stringsAsFactors = FALSE
) %>%
  filter(!(rho_case == "alt" & delta_case == "null"))

# === Simulation parameters ===
N_set <- c(80, 120, 160, 200)
T_eval_set <- c(20, 50, 100, 150)
n_sim <- 200
burn <- 50
lrv <- "EWC"
lrv_par <- NULL
lambda <- 0.2

set.seed(12)

# === Setup parallel backend ===
n_cores <- parallel::detectCores(logical = FALSE)
cl <- makeCluster(n_cores)
registerDoSNOW(cl)

# === Initialize full result list ===
all_results <- list()

# === Loop over design grid ===
for (s in seq_len(nrow(param_grid))) {
  rho_case <- param_grid$rho_case[s]
  delta_case <- param_grid$delta_case[s]
  phi_case <- param_grid$phi_case[s]
  
  rho_vec <- if (rho_case == "null") rep(0, 3) else c(0.1, 0.2, 0.3)
  delta_vec <- if (delta_case == "null") rep(0, 3) else c(-0.3, -0.1, 0.2)
  phi <- if (phi_case == "null") 0 else 0.2
  
  test_case_label <- paste0("rho_", rho_case, "+delta_", delta_case, "+phi_", phi_case)
  cat("\n=== Running Test Case:", test_case_label, "===\n")
  
  results_matrix <- matrix(NA, nrow = length(N_set), ncol = length(T_eval_set) * 4)
  colnames(results_matrix) <- paste(rep(c("Overall", "Clustered", "Overall_Cond", "Clustered_Cond"), each = length(T_eval_set)),
                                    rep(paste0("T=", T_eval_set), times = 4), sep = ":")
  rownames(results_matrix) <- paste0("N=", N_set)
  
  pb <- txtProgressBar(max = length(N_set) * length(T_eval_set) * n_sim, style = 3)
  counter <- 0
  progress <- function(n) {
    counter <<- counter + 1
    setTxtProgressBar(pb, counter)
  }
  opts <- list(progress = progress)
  
  for (n_idx in seq_along(N_set)) {
    N <- N_set[n_idx]
    gamma <- rep(NA, N)
    cut1 <- floor(N / 4)
    cut2 <- floor(N / 2)
    gamma[1:cut1] <- 1
    gamma[(cut1 + 1):cut2] <- 2
    gamma[(cut2 + 1):N] <- 3
    
    for (t_idx in seq_along(T_eval_set)) {
      T_eval <- T_eval_set[t_idx]
      T_adj <- T_eval - 1  # Adjust for lag
      
      res_mat <- foreach(sim = 1:n_sim, .combine = rbind,
                         .packages = c("stats", "base"),
                         .options.snow = opts) %dopar% {
                           sim_data <- generate_forecast_simulation_data(
                             N, Tobs = T_eval + 1, burn,
                             alpha = 1,
                             rho_vec = rho_vec,
                             gamma = gamma,
                             phi = phi,
                             lambda = lambda,
                             delta_vec = delta_vec
                           )
                           
                           Y <- sim_data$Y
                           forecast_1 <- sim_data$forecast_1
                           forecast_2 <- sim_data$forecast_2
                           
                           L1 <- (Y[, 2:T_eval] - forecast_1[, 2:T_eval])^2
                           L2 <- (Y[, 2:T_eval] - forecast_2[, 2:T_eval])^2
                           DL <- matrix(as.vector(L1 - L2), ncol = 1)
                           colnames(DL) <- "DL"
                           
                           id_vec <- rep(1:N, each = T_adj)
                           time_vec <- rep(1:T_adj, times = N)
                           cond_var <- as.vector(Y[, 1:(T_eval - 1)])
                           cond_mat <- matrix(cond_var, ncol = 1)
                           colnames(cond_mat) <- "Z"
                           
                           res_epa <- overall_EPA_test(Z = DL, id = id_vec, time = time_vec,
                                                       lrv = lrv, lrv_par = lrv_par)
                           rej_epa <- as.numeric(res_epa$p_oepa <= 0.05)
                           
                           res_clus <- epa_clustered_known(Z = DL, id = id_vec, time = time_vec, gamma = gamma,
                                                           lrv = lrv, lrv_par = lrv_par)
                           rej_clus <- as.numeric(res_clus$pval <= 0.05)
                           
                           res_epa_cond <- overall_EPA_test(Z = cbind(DL, cond_mat), id = id_vec, time = time_vec,
                                                            lrv = lrv, lrv_par = lrv_par)
                           rej_epa_cond <- as.numeric(res_epa_cond$p_oepa <= 0.05)
                           
                           res_clus_cond <- epa_clustered_known(Z = cbind(DL, cond_mat), id = id_vec, time = time_vec,
                                                                gamma = gamma, lrv = lrv, lrv_par = lrv_par)
                           rej_clus_cond <- as.numeric(res_clus_cond$pval <= 0.05)
                           
                           return(c(rej_epa, rej_clus, rej_epa_cond, rej_clus_cond))
                         }
      
      results_matrix[n_idx, ((t_idx - 1) * 4 + 1):(t_idx * 4)] <- colMeans(res_mat)
    }
  }
  
  all_results[[test_case_label]] <- results_matrix
}

stopCluster(cl)

# === Write to Excel ===
out_path <- here("Tools", "clusteredEPA", "Results", "conditional_epa_results_all_in_one.xlsx")
wb <- createWorkbook()
addWorksheet(wb, "All Results")

# Format and combine into single sheet
long_result <- do.call(rbind, lapply(names(all_results), function(name) {
  df <- all_results[[name]]
  df_case <- as.data.frame(df)
  df_case$N <- rownames(df)
  df_case$Design <- name
  df_long <- tidyr::pivot_longer(df_case, cols = -c(N, Design), names_to = "Test_T", values_to = "Rejection")
  df_long
}))

writeDataTable(wb, sheet = "All Results", x = long_result)
saveWorkbook(wb, file = out_path, overwrite = TRUE)
