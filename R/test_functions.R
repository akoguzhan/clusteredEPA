#' Overall Equal Predictive Ability (O-EPA) Test for Panel Data
#'
#' @description Computes the Overall EPA test from NT × P matrix input (panel structure),
#' averaging across units for each time period. The function supports both EWC and Newey-West LRV estimators.
#'
#' @param Z NT × P matrix containing loss differential and/or transformed moment columns.
#' @param id Vector of unit identifiers (length NT).
#' @param time Vector of time identifiers (length NT).
#' @param lrv Character: either "NeweyWest" or "EWC" for long-run variance estimation.
#' @param lrv_par Numeric tuning parameter: lag (NW) or degrees of freedom (EWC).
#'
#' @return A list with:
#'   \describe{
#'     \item{S_oepa}{Long-run variance estimator.}
#'     \item{W_oepa}{Wald test statistic.}
#'     \item{p_oepa}{P-value of the test.}
#'   }
#' @export
overall_EPA_test <- function(Z, id, time, lrv, lrv_par) {
  # Checks
  if (!is.matrix(Z)) stop("Z must be a matrix.")
  NT <- nrow(Z)
  if (length(id) != NT || length(time) != NT) {
    stop("Length of 'id' and 'time' must match number of rows in Z.")
  }
  
  # Unique panel dimension
  Tobs <- length(unique(time))
  N <- NT / Tobs
  P <- ncol(Z)
  
  # Edge case: not enough time periods for the number of moments
  if (Tobs < P) stop("Number of time periods (Tobs) must be at least as large as the number of moments (P).")
  
  # Average across cross-section for each time
  Zbar_Nt <- aggregate_matrix(Z, time)

  # Sample mean
  mu_hat_oepa <- colMeans(Zbar_Nt)

  # Test statistic and p-value
  if (identical(lrv, "NeweyWest")) {
    S_oepa <- NeweyWest(Zbar_Nt, lrv_par)
    # Use Cholesky-based inversion for numerical stability
    S_inv <- tryCatch(chol2inv(chol(S_oepa)), error = function(e) stop("S_oepa is not positive definite."))
    W_oepa <- as.numeric(Tobs * t(mu_hat_oepa) %*% S_inv %*% mu_hat_oepa)
    p_oepa <- pchisq(W_oepa, df = P, lower.tail = FALSE)
  } else if (identical(lrv, "EWC")) {
    EWC_res <- EWC(Zbar_Nt, lrv_par)
    S_oepa <- EWC_res$S
    lrv_par <- EWC_res$lrv_par
    adj_term <- ((lrv_par - P + 1) / lrv_par) / P
    # Use Cholesky-based inversion for numerical stability
    S_inv <- tryCatch(chol2inv(chol(S_oepa)), error = function(e) stop("S_oepa is not positive definite."))
    W_oepa <- as.numeric(adj_term * Tobs * t(mu_hat_oepa) %*% S_inv %*% mu_hat_oepa)
    p_oepa <- pf(W_oepa, df1 = P, df2 = lrv_par - P + 1, lower.tail = FALSE)
  } else {
    stop("Invalid 'lrv' value. Use 'NeweyWest' or 'EWC'.")
  }
  
  return(list(S_oepa = S_oepa, W_oepa = W_oepa, p_oepa = p_oepa))
}

#' Clustered EPA (C-EPA) Test with Known Group Memberships
#'
#' @param Z Matrix (NT x P) of loss differentials or moment conditions.
#' @param id Vector of length NT for unit identifiers.
#' @param time Vector of length NT for time identifiers.
#' @param gamma Vector of length N (number of units), giving cluster assignments (values in {1, ..., K}).
#' @param lrv Character; "EWC" (default) or "NeweyWest".
#' @param lrv_par Integer; lag (NeweyWest) or number of spectral windows (EWC).
#'
#' @return A list with elements:
#'   \describe{
#'     \item{W}{Wald statistic}
#'     \item{pval}{P-value}
#'   }
#' @export
epa_clustered_known <- function(Z, id, time, gamma, lrv = "EWC", lrv_par = NULL) {
  if (!is.matrix(Z)) stop("Z must be a matrix.")
  NT <- nrow(Z)
  P <- ncol(Z)
  N <- length(unique(id))
  Tobs <- length(unique(time))
  
  # Validate gamma
  if (length(gamma) != N) stop("Length of gamma must match number of unique units.")
  
  # Combine and sort
  df <- data.frame(Z, id = id, time = time)
  df <- df[order(df$time, df$id), ]
  Z <- as.matrix(df[, 1:P])
  id <- df$id
  time <- df$time
  
  # Expand gamma to NT
  id_unique <- sort(unique(id))
  gamma_long <- gamma[match(id, id_unique)]
  
  # Generate time-cluster means
  K <- length(unique(gamma))
  KP <- K * P
  if (Tobs < KP) stop("Number of time periods (Tobs) must be at least as large as K*P.")
  Z_by_group <- matrix(NA, nrow = Tobs, ncol = KP)
  
  for (k in 1:K) {
    idx_k <- gamma_long == k
    Z_k <- Z[idx_k, , drop = FALSE]
    time_k <- time[idx_k]
    Zbar_k <- aggregate_matrix(Z_k, time_k)  # Tobs x P
    Z_by_group[, ((k - 1) * P + 1):(k * P)] <- Zbar_k
  }
  
  # Compute test statistic
  mu_hat <- colMeans(Z_by_group)
  
  if (lrv == "NeweyWest") {
    Omega <- NeweyWest(Z_by_group, lrv_par = lrv_par)
    Omega_inv <- tryCatch(chol2inv(chol(Omega)), error = function(e) stop("Omega is not positive definite."))
    W <- as.numeric(Tobs * t(mu_hat) %*% Omega_inv %*% mu_hat)
    pval <- pchisq(W, df = KP, lower.tail = FALSE)
  } else if (lrv == "EWC") {
    EWC_res <- EWC(Z_by_group, lrv_par = lrv_par)
    Omega <- EWC_res$S
    lrv_par <- EWC_res$lrv_par
    adj_term <- ((lrv_par - KP + 1) / lrv_par) / KP
    Omega_inv <- tryCatch(chol2inv(chol(Omega)), error = function(e) stop("Omega is not positive definite."))
    W <- as.numeric(adj_term * Tobs * t(mu_hat) %*% Omega_inv %*% mu_hat)
    pval <- pf(W, df1 = KP, df2 = lrv_par - KP + 1, lower.tail = FALSE)
  } else {
    stop("Unknown lrv. Use 'EWC' or 'NeweyWest'.")
  }
  
  return(list(W = W, pval = pval))
}

#' Split-Sample Clustered EPA
#'
#' @param df Balanced panel data.frame
#' @param id Character; unit identifier column name
#' @param time Character; time identifier column name
#' @param Z_names Character vector; variable names for clustering
#' @param K Integer; number of clusters for KMeans
#' @param prop Numeric; proportion of time periods for R (training) sample
#' @param lrv Character; "EWC" or "NeweyWest" for LRV estimator
#' @param param Integer; lag (NW) or B (EWC)
#' @param Ninit Integer; number of KMeans initializations
#' @param iter.max Integer; max iterations in KMeans
#'
#' @return List with cluster assignments and test results from P
#' @export 
epa_clustered_split <- function(df, id, time, Z_names,
                                K, prop = 0.5,
                                lrv = "EWC", lrv_par = NULL,
                                Ninit = 10, iter.max = 10) {
  # Step 1: Convert to matrices
  matrices <- panel_data_to_matrix(df, id = id, time = time, Z_names = Z_names)
  Z <- matrices$Z_panel
  id_vec <- matrices$id
  time <- matrices$time

  # Step 2: Split matrices into R and P sets
  split <- split_panel_matrix(Z, id = id_vec, time = time, prop = prop)
  Z_R <- split$Z_R
  Z_P <- split$Z_P
  id_R <- split$id_R
  id_P <- split$id_P
  time_P <- split$time_P
  time_R <- split$time_R

  # Step 3: Run k-means clustering on R sample
  km_res <- panel_kmeans_estimation(Z = Z_R, id = id_R, time = time_R, K = K,
                                    Ninit = Ninit, iter.max = iter.max)
  gamma_R <- km_res$final_cluster

  # Step 4: Map cluster assignments from R to P by id
  id_map <- setNames(gamma_R, unique(id_R))
  gamma_P <- id_map[as.character(unique(id_P))]
  if (any(is.na(gamma_P))) stop("Some units in P were not found in R clustering.")

  # Step 5: Run clustered EPA test on P using clusters from R
  test_res <- epa_clustered_known(
    Z = Z_P,
    id = id_P,
    time = time_P,
    gamma = gamma_P,
    lrv = lrv,
    lrv_par = lrv_par
  )

  return(c(list(clustering = gamma_R), unlist(test_res) ))
}

#' Selective Inference for Clustered EPA
#'
#' Performs selective inference for clustered Equal Predictive Ability (EPA):
#' 1. Runs panel_kmeans_estimation to get cluster assignments (with K or Kmax).
#' 2. Calculates the long-run variance matrix internally.
#' 3. Runs panel_homogeneity_test to get all pairwise cluster p-values (with flexible combination).
#' 4. Runs overall_EPA_test for the global null.
#' 5. Combines all [K(K-1)/2 + 1] p-values using the specified method (default: pharmonic).
#'
#' @param Z NT x P matrix of loss differentials or moments.
#' @param id Vector of unit identifiers (length NT).
#' @param time Vector of time identifiers (length NT).
#' @param K Integer; number of clusters for k-means. Required unless Kmax is provided.
#' @param Kmax Optional; if provided, perform BIC-based selection from 2 to Kmax clusters.
#' @param pairs Optional matrix of cluster index pairs to test (each row: c(k1, k2)). If NULL, all unique pairs are tested.
#' @param pcombine_fun Function to combine p-values: "pmean", "porder", "pSimes", "pharmonic", "pCauchy".
#' @param method Variant used by pmerge (e.g., "H1", "I", "G", etc.).
#' @param order_k Used by porder.
#' @param r Used by pmean.
#' @param epi Used by pCauchy.
#' @param lrv Character; "EWC" or "NeweyWest" for overall_EPA_test.
#' @param lrv_par Integer; lag (NW) or B (EWC) for overall_EPA_test.
#' @param Ninit Integer; number of KMeans initializations.
#' @param iter.max Integer; max iterations in KMeans.
#' @return List with all pairwise p-values, overall p-value, and selective p-value.
#' @export 
epa_clustered_selective <- function(
    Z, id, time,
    K = NULL,
    Kmax = NULL,
    pairs = NULL,
    pcombine_fun = "pmean",
    method = "A",
    order_k = NULL,
    r = NULL,
    epi = NULL,
    lrv = "EWC",
    lrv_par = NULL,
    Ninit = 10,
    iter.max = 10
) {
  # 1. Estimate clusters
  estimated_k_means <- PanelKmeansInference::panel_kmeans_estimation(
    Z = Z, id = id, time = time,
    K = K, Kmax = Kmax,
    Ninit = Ninit, iter.max = iter.max
  )
  gamma <- estimated_k_means$final_cluster
  K_used <- if (!is.null(K)) K else estimated_k_means$BIC_selected_K
  N <- length(unique(id))
  Tobs <- length(unique(time))
  P <- ncol(Z)
  KP <- K_used * P

  # 2. Calculate long-run variance matrix (OmegaHat) for the clustered means
  # Expand gamma to NT
  id_unique <- sort(unique(id))
  gamma_long <- gamma[match(id, id_unique)]
  # Generate time-cluster means
  Z_by_group <- matrix(NA, nrow = Tobs, ncol = KP)
  for (k in 1:K_used) {
    idx_k <- gamma_long == k
    Z_k <- Z[idx_k, , drop = FALSE]
    time_k <- time[idx_k]
    Zbar_k <- aggregate_matrix(Z_k, time_k)  # Tobs x P
    Z_by_group[, ((k - 1) * P + 1):(k * P)] <- Zbar_k
  }
  if (lrv == "NeweyWest") {
    OmegaHat <- NeweyWest(Z_by_group, lrv_par = lrv_par)
  } else if (lrv == "EWC") {
    OmegaHat <- EWC(Z_by_group, lrv_par = lrv_par)$S
  } else {
    stop("Unknown lrv. Use 'EWC' or 'NeweyWest'.")
  }

  # 3. Pairwise cluster p-values (with flexible combination)
  homo_res <- panel_homogeneity_test(
    Z = Z,
    id = id,
    time = time,
    OmegaHat = OmegaHat,
    estimated_k_means = estimated_k_means,
    pairs = pairs,
    pcombine_fun = pcombine_fun,
    method = method,
    order_k = order_k,
    r = r,
    epi = epi
  )

  pairwise_pvals <- homo_res$pairwise_pvalues
  pairs_out <- homo_res$pairs
  pval_combined = homo_res$pvalue_combination

  # 4. Overall EPA test p-value
  oepa <- overall_EPA_test(Z = Z, id = id, time = time, lrv = lrv, lrv_par = lrv_par)
  overall_pval <- oepa$p_oepa

  # 5. Combine all p-values (pairwise + overall)
  all_pvals <- c(as.numeric(pairwise_pvals), overall_pval)
  names(all_pvals) <- c(
    paste0("pair_", sapply(pairs_out, function(x) paste(x, collapse = "_"))),
    "overall"
  )

  # Use the same combination method as for pairwise, but now on all p-values
  selective_pval <- switch(
    pcombine_fun,
    pmean     = pmerge::pmean(p = all_pvals, r = r, dependence = method),
    porder    = pmerge::porder(p = all_pvals, k = order_k),
    pSimes    = pmerge::pSimes(p = all_pvals, method = method),
    pharmonic = pmerge::pharmonic(p = all_pvals, method = method),
    stop("Invalid pcombine_fun specified.")
  )

  return(list(
    pairwise_pvalues = pairwise_pvals,
    pairs = pairs_out,
    overall_pval = overall_pval,
    selective_pval = selective_pval,
    pcombine_fun = pcombine_fun,
    method = method,
    pval_combined_pairwise = pval_combined,
    cluster_assignments = gamma,
    estimated_k_means = estimated_k_means,
    OmegaHat = OmegaHat
  ))
}

#' Clark and West (2007) MSPE-Adjusted Forecast Comparison Test
#'
#' @description Implements the Clark and West (2007, JOE) test to compare forecast accuracy between nested models.
#'
#' @param y1 Numeric vector of forecasts from model 1 (length T).
#' @param y2 Numeric vector of forecasts from model 2 (length T).
#' @param y  Numeric vector of actual observations (length T).
#'
#' @return A list with elements \code{T} (test statistic) and \code{pval} (p-value).
#' @export
CWtest <- function(y1, y2, y, lrv_par = 0) {
  if (!(length(y1) == length(y2) && length(y2) == length(y))) {
    stop("Input vectors y1, y2, and y must have the same length.")
  }
  P <- length(y1)
  
  # Compute squared prediction errors and f statistic
  e1 <- (y - y1)^2
  e2 <- (y - y2)^2
  f <- e1 - (e2 - (y1 - y2)^2)
  
  # OLS estimate of mean of f (mu)
  mu <- sum(f) / P
  
  # Compute long-run variance of f
  sigma2 <- NeweyWest(f, lrv_par = lrv_par)
  
  # Test statistic and p-value
  T_stat <- mu / sqrt(sigma2 / P)
  pval <- 1 - pnorm(T_stat)
  
  return(list(T = T_stat, pval = pval))
}

#' Giacomini and White (2006) Predictive Ability Test
#'
#' @description Implements the conditional predictive ability (CPA) test for two competing forecasts
#' using the time series of out-of-sample loss differentials and a Newey-West long-run variance estimate.
#'
#' @param d_t A numeric vector of out-of-sample loss differentials (L1 - L2) of length P.
#' @param x_t Optional matrix or vector of instruments (predictors) to condition on. Defaults to 1 (unconditional test).
#' @param max_lag Integer; lag length for Newey-West LRV estimation.
#'
#' @return A list with the test statistic, p-value, and long-run variance estimate.
#' @export
GWtest <- function(d_t, x_t = NULL, max_lag = 0) {
  if (!is.vector(d_t)) stop("d_t must be a vector")
  if (anyNA(d_t)) stop("d_t contains missing values")
  P <- length(d_t)
  
  # Default to unconditional mean if no instruments are given
  if (is.null(x_t)) {
    x_t <- matrix(1, nrow = P, ncol = 1)
  } else {
    x_t <- as.matrix(x_t)
    if (nrow(x_t) != P) stop("x_t must have same number of rows as d_t")
    if (anyNA(x_t)) stop("x_t contains missing values")
  }
  
  # Cross-product matrix
  X <- x_t
  Y <- d_t
  # Use qr.solve for stability
  beta_hat <- qr.solve(t(X) %*% X, t(X) %*% Y)
  e_t <- Y - X %*% beta_hat
  
  # Long-run variance estimate of moment vector (Newey-West)
  S_hat <- NeweyWest(X * e_t, lrv_par = max_lag)
  
  # Test statistic
  stat <- t(beta_hat) %*% solve(S_hat) %*% beta_hat
  df <- ncol(X)
  pval <- 1 - pchisq(stat, df)
  
  return(list(statistic = as.numeric(stat), pval = as.numeric(pval), S_hat = S_hat))
}