rm(list = ls())

# === Load Packages and Project Code ===
library(devtools)
library(dplyr)
library(openxlsx)
library(mclust)
library(ggplot2)
library(tidyr)
library(parallel)
library(pbapply)

# === Load Your Project Packages (adjust paths as needed) ===
devtools::load_all("~/clusteredEPA/Tools/PanelKmeansInference")
devtools::load_all("~/clusteredEPA/Tools/clusteredEPA")

# === Simulation Parameters ===
N_set <- c(80, 120, 160)
T_eval_set <- c(20, 50, 100, 200)
d_vals <- c(0, 0.125, 0.25, 0.375, 0.5)
n_sim <- 3
burn <- 10
lrv <- "EWC"
lrv_par <- NULL

# === Load existing results if available ===
results_clustered <- if (file.exists("clustered_partial.rds")) readRDS("clustered_partial.rds") else list()
results_misc <- if (file.exists("misc_partial.rds")) readRDS("misc_partial.rds") else list()
all_results <- if (file.exists("full_partial.rds")) readRDS("full_partial.rds") else list()

set.seed(1)

completed_keys <- lapply(all_results, function(x) paste(x$N, x$T, x$d, sep = ","))

progress_total <- length(N_set) * length(T_eval_set) * length(d_vals)
progress_counter <- 0

for (N in N_set) {
  gamma <- rep(NA, N)
  cut1 <- floor(N / 4)
  cut2 <- floor(N / 2)
  gamma[1:cut1] <- 1
  gamma[(cut1 + 1):cut2] <- 2
  gamma[(cut2 + 1):N] <- 3

  for (T_eval in T_eval_set) {
    id_vec <- rep(1:N, each = T_eval)
    time_vec <- rep(1:T_eval, times = N)

    for (d in d_vals) {
      key <- paste(N, T_eval, d, sep = ",")
      if (key %in% completed_keys) {
        next
      }

      set.seed(N + T_eval * 100 + round(d * 10000))

      progress_counter <- progress_counter + 1
      cat(sprintf("(%d/%d) Running N = %d, T = %d, d = %.3f\n", 
                  progress_counter, progress_total, N, T_eval, d))

      delta_vec <- d + d * c(-0.05, -0.1, 0.15)

      n_cores <- max(1, detectCores(logical = FALSE) - 1)
      cl <- makeCluster(n_cores)
      clusterExport(cl, varlist = ls(), envir = environment())
      clusterEvalQ(cl, {
        library(devtools)
        devtools::load_all("~/clusteredEPA/Tools/PanelKmeansInference")
        devtools::load_all("~/clusteredEPA/Tools/clusteredEPA")
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

        L1 <- (Y[, 2:(T_eval+1)] - forecast_1[, 2:(T_eval+1)])^2
        L2 <- (Y[, 2:(T_eval+1)] - forecast_2[, 2:(T_eval+1)])^2
        DL <- matrix(as.vector(L1 - L2), ncol = 1)
        cond_var <- as.vector(Y[, 1:T_eval])
        X <- matrix(cond_var, ncol = 1)

        df <- matrix_to_panel_data(cbind(DL, X), Z_names = c("DL", "X"), id_vec, time_vec)

        out_all <- list(
          known = clusteredEPA(df, "DL", NULL, "id", "time", test = "epa_clustered_known",gamma = gamma, lrv = lrv, lrv_par = lrv_par),
          known_cond = clusteredEPA(df, "DL", "X", "id", "time", test = "epa_clustered_known", gamma = gamma, lrv = lrv, lrv_par = lrv_par),
          naive = clusteredEPA(df, "DL", NULL, "id", "time", test = "epa_clustered_known", gamma = panel_kmeans_estimation(DL, id_vec, time_vec, Kmax = 5, n_cores = 1)$final_cluster, lrv = lrv, lrv_par = lrv_par),
          naive_cond = clusteredEPA(df, "DL", "X", "id", "time", test = "epa_clustered_known", gamma = panel_kmeans_estimation(cbind(DL, X), id_vec, time_vec, Kmax = 5, n_cores = 1)$final_cluster, lrv = lrv, lrv_par = lrv_par),
          split = clusteredEPA(df, "DL", NULL, "id", "time", Kmax = 5, test = "epa_clustered_split", prop = 0.2, lrv = lrv, lrv_par = lrv_par),
          split_cond = clusteredEPA(df, "DL", "X", "id", "time", Kmax = 5, test = "epa_clustered_split", prop = 0.2, lrv = lrv, lrv_par = lrv_par),
          selective = clusteredEPA(df, "DL", NULL, "id", "time", Kmax = 5, test = "epa_clustered_selective", pcombine_fun = "pGridIU", lrv = lrv, lrv_par = lrv_par),
          selective_cond = clusteredEPA(df, "DL", "X", "id", "time", Kmax = 5, test = "epa_clustered_selective", pcombine_fun = "pGridIU", lrv = lrv, lrv_par = lrv_par),
          overall = clusteredEPA(df, "DL", NULL, "id", "time", test = "overall_EPA_test", lrv = lrv, lrv_par = lrv_par),
          overall_cond = clusteredEPA(df, "DL", "X", "id", "time", test = "overall_EPA_test", lrv = lrv, lrv_par = lrv_par)
        )

        rejection_rates <- c(sapply(out_all, function(x) as.numeric(x$pval <= 0.05)), out_all[["selective"]]$hom_test_pval, out_all[["selective_cond"]]$hom_test_pval)
        names(rejection_rates)[(length(rejection_rates)-1):length(rejection_rates)] <- c("hom", "hom_cond")

        extract_clustering_stats <- function(obj, gamma) {
          if (is.null(obj) || is.null(obj$clustering)) {
            return(c(rand = NA, recovery = NA, avg_K = NA))
          }
          cl <- obj$clustering
          rand <- mclust::adjustedRandIndex(cl, gamma)
          recovery <- same_cl(cl, gamma, 3)
          Khat <- length(unique(cl))
          avg_K <- Khat
          return(c(rand = rand, recovery = recovery, avg_K = avg_K))
        }

        clustering_stats <- c(
          setNames(extract_clustering_stats(out_all$split, gamma), 
                   c("rand_split", "recovery_split", "avg_K_split")),
          setNames(extract_clustering_stats(out_all$split_cond, gamma), 
                   c("rand_split_cond", "recovery_split_cond", "avg_K_split_cond")),
          setNames(extract_clustering_stats(out_all$selective, gamma), 
                   c("rand_selective", "recovery_selective", "avg_K_selective")),
          setNames(extract_clustering_stats(out_all$selective_cond, gamma), 
                   c("rand_selective_cond", "recovery_selective_cond", "avg_K_selective_cond"))
        )

        out_df <- data.frame(N = N, T = T_eval, d = d, as.list(rejection_rates), as.list(clustering_stats))
        return(out_df)
      }

      sim_results <- pblapply(1:n_sim, sim_fun, cl = cl)
      stopCluster(cl)

      res_df <- as.data.frame(t(colMeans(do.call(rbind, sim_results), na.rm = TRUE)))
      results_clustered[[length(results_clustered) + 1]] <- res_df %>% 
        select(N, T, d, known, known_cond, naive, naive_cond, split, split_cond, selective, selective_cond)
      results_misc[[length(results_misc) + 1]] <- res_df %>% 
        select(N, T, d, overall, overall_cond, hom, hom_cond, starts_with("rand_"), starts_with("recovery_"), starts_with("avg_K_"))
      all_results[[length(all_results) + 1]] <- res_df

      saveRDS(results_clustered, file = "clustered_partial.rds")
      saveRDS(results_misc, file = "misc_partial.rds")
      saveRDS(all_results, file = "full_partial.rds")
    }
  }
}

final_clustered <- bind_rows(results_clustered)
final_misc <- bind_rows(results_misc)
final_all <- bind_rows(all_results)

write.xlsx(final_clustered, file = "results_clustered.xlsx")
write.xlsx(final_misc, file = "results_misc.xlsx")
write.xlsx(final_all, file = "results_full_raw.xlsx")