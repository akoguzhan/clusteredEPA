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
N_set <- 80
T_eval_set <- c(20, 50, 100, 200)
d_vals <- 0.25
n_sim <- 1000
burn <- 10
lrv <- "EWC"
lrv_par <- NULL

set.seed(1)

results_clustered <- list()
results_misc <- list()
results_pcomb <- list()

for (d in d_vals) {
  variant_list <- list(
    `scenario1` = -(d/2 + d * c(-1.2, -0.8, 1)),
    `scenario2` = c(0, 0, 0)
  )
  
  for (variant in names(variant_list)) {
    for (N in N_set) {
      for (T_eval in T_eval_set) {
        delta_vec2 <- variant_list[[variant]]
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
        clusterExport(cl, varlist = ls(), envir = environment())
        clusterEvalQ(cl, {
          devtools::load_all("C:/Users/og7051ak/Nextcloud/Codes for Papers/Tests with Unknown Clusters/Tools/PanelKmeansInference")
          devtools::load_all("C:/Users/og7051ak/Nextcloud/Codes for Papers/Tests with Unknown Clusters/Tools/clusteredEPA")
          library(mclust)
        })
        
        sim_fun <- function(sim) {
          sim_data_1 <- generate_forecast_simulation_data(
            N, Tobs = T_eval/2 + 1, burn,
            alpha = 1, rho_vec = c(0.1, 0.2, 0.3), gamma = gamma,
            phi = 0.2, lambda = 0.2, delta_vec = d/2 + d * c(-1.2, -0.8, 1)
          )
          
          sim_data_2 <- generate_forecast_simulation_data(
            N, Tobs = T_eval/2 + 2, burn,
            alpha = 1, rho_vec = c(0.1, 0.2, 0.3), gamma = gamma,
            phi = 0.2, lambda = 0.2, delta_vec = delta_vec2
          )
          
          Y <- cbind(sim_data_1$Y, sim_data_2$Y)
          forecast_1 <- cbind(sim_data_1$forecast_1, sim_data_2$forecast_1)
          forecast_2 <- cbind(sim_data_1$forecast_2, sim_data_2$forecast_2)
          
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
            selective_cond = clusteredEPA(df, "DL", "X", "id", "time", Kmax = 5, test = "epa_clustered_selective", pcombine_fun = "Genmean_neq", lrv = lrv, lrv_par = lrv_par)
          )
          
          rejection_rates <- c(sapply(out_all, function(x) as.numeric(x$pval <= 0.05)),
                               as.numeric(out_all[["selective"]]$hom_test_pval <= 0.05),
                               as.numeric(out_all[["selective_cond"]]$hom_test_pval <= 0.05),
                               as.numeric(out_all[["selective"]]$overall_pval <= 0.05),
                               as.numeric(out_all[["selective_cond"]]$overall_pval <= 0.05))
          names(rejection_rates)[(length(rejection_rates)-3):length(rejection_rates)] <- c("hom", "hom_cond", "overall", "overall_cond")
          
          allpvalues <- c(out_all[["selective"]]$pairwise_pvalues, out_all[["selective"]]$overall_pval)
          allpvalues_cond <- c(out_all[["selective_cond"]]$pairwise_pvalues, out_all[["selective_cond"]]$overall_pval)
          
          pval_combs <- as.numeric(
            c(cauchy_pcombine(allpvalues), cauchy_pcombine(allpvalues_cond),
              bonferroni_pcombine(allpvalues), bonferroni_pcombine(allpvalues_cond),
              iu_pcombine(allpvalues), iu_pcombine(allpvalues_cond),
              Genmean_rneg_pcombine(allpvalues), Genmean_rneg_pcombine(allpvalues_cond),
              Genmean_pcombine(allpvalues), Genmean_pcombine(allpvalues_cond),
              Geomean_pcombine(allpvalues), Geomean_pcombine(allpvalues_cond),
              bonferroni_compound_pcombine(allpvalues, bonferroni_pcombine, cauchy_pcombine), bonferroni_compound_pcombine(allpvalues_cond, bonferroni_pcombine, cauchy_pcombine),
              bonferroni_compound_pcombine(allpvalues, bonferroni_pcombine, iu_pcombine), bonferroni_compound_pcombine(allpvalues_cond, bonferroni_pcombine, iu_pcombine),
              bonferroni_compound_pcombine(allpvalues, bonferroni_pcombine, Genmean_rneg_pcombine), bonferroni_compound_pcombine(allpvalues_cond, bonferroni_pcombine, Genmean_rneg_pcombine),
              bonferroni_compound_pcombine(allpvalues, bonferroni_pcombine, Genmean_pcombine), bonferroni_compound_pcombine(allpvalues_cond, bonferroni_pcombine, Genmean_pcombine),
              bonferroni_compound_pcombine(allpvalues, bonferroni_pcombine, Geomean_pcombine), bonferroni_compound_pcombine(allpvalues_cond, bonferroni_pcombine, Geomean_pcombine))
            <= 0.05)
          
          names(pval_combs) <- c("cauchy_pcombine", "cauchy_cond_pcombine",
                                 "bonferroni_pcombine", "bonferroni_cond_pcombine",
                                 "iu_pcombine", "iu_cond_pcombine",
                                 "Genmean_rneg_pcombine", "Genmean_rneg_cond_pcombine",
                                 "Genmean_pcombine", "Genmean_cond_pcombine",
                                 "Geomean_pcombine", "Geomean_cond_pcombine",
                                 "cauchy_bonf_pcombine", "cauchy_bonf_cond_pcombine",
                                 "iu_bonf_pcombine", "iu_bonf_cond_pcombine",
                                 "Genmean_rneg_bonf_pcombine", "Genmean_rneg_bonf_cond_pcombine",
                                 "Genmean_bonf_pcombine", "Genmean_bonf_cond_pcombine",
                                 "Geomean_bonf_pcombine", "Geomean_bonf_cond_pcombine")
          
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
            setNames(extract_clustering_stats(out_all$split, gamma), c("rand_split", "recovery_split", "avg_K_split")),
            setNames(extract_clustering_stats(out_all$split_cond, gamma), c("rand_split_cond", "recovery_split_cond", "avg_K_split_cond")),
            setNames(extract_clustering_stats(out_all$selective, gamma), c("rand_selective", "recovery_selective", "avg_K_selective")),
            setNames(extract_clustering_stats(out_all$selective_cond, gamma), c("rand_selective_cond", "recovery_selective_cond", "avg_K_selective_cond"))
          )
          
          out_df <- data.frame(N = N, Tobs = T_eval, d = d, variant = variant,
                               as.list(rejection_rates), as.list(clustering_stats), as.list(pval_combs))
          return(out_df)
        }
        
        sim_results <- pblapply(1:n_sim, sim_fun, cl = cl)
        stopCluster(cl)
        
        combined_df <- bind_rows(sim_results)
        numeric_means <- combined_df %>% select(where(is.numeric)) %>% summarise(across(everything(), ~ mean(.x, na.rm = TRUE)))
        variant <- combined_df[1, "variant"]
        res_df <- cbind(variant, numeric_means)
        
        results_clustered[[length(results_clustered) + 1]] <- res_df %>% select(N, Tobs, d, variant, known, known_cond, naive, naive_cond, split, split_cond, selective, selective_cond)
        results_misc[[length(results_misc) + 1]] <- res_df %>% select(N, Tobs, d, variant, overall, overall_cond, hom, hom_cond, starts_with("rand_"), starts_with("recovery_"), starts_with("avg_K_"))
        results_pcomb[[length(results_pcomb) + 1]] <- res_df %>% select(N, Tobs, d, variant, ends_with("_pcombine"))
        
        write.xlsx(bind_rows(results_clustered), here("Tools", "clusteredEPA", "Results", "results_clustered_sbreaks.xlsx"))
        write.xlsx(bind_rows(results_misc), here("Tools", "clusteredEPA", "Results", "results_misc_sbreaks.xlsx"))
        write.xlsx(bind_rows(results_pcomb), here("Tools", "clusteredEPA", "Results", "results_pcomb_sbreaks.xlsx"))
      }
    }
  }
}
