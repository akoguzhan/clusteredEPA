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
#'     \item{S}{Long-run variance estimator.}
#'     \item{W}{Wald test statistic.}
#'     \item{pval}{P-value of the test.}
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
    S <- NeweyWest(Zbar_Nt, lrv_par)
    # Use Cholesky-based inversion for numerical stability
    S_inv <- chol2inv(chol(S))
    W <- as.numeric(Tobs * t(mu_hat_oepa) %*% S_inv %*% mu_hat_oepa)
    pval <- pchisq(W, df = P, lower.tail = FALSE)
  } else if (identical(lrv, "EWC")) {
    EWC_res <- EWC(Zbar_Nt, lrv_par)
    S <- EWC_res$S
    lrv_par <- EWC_res$lrv_par
    adj_term <- ((lrv_par - P + 1) / lrv_par) / P
    # Use Cholesky-based inversion for numerical stability
    S_inv <- chol2inv(chol(S))
    W <- as.numeric(adj_term * Tobs * t(mu_hat_oepa) %*% S_inv %*% mu_hat_oepa)
    pval <- pf(W, df1 = P, df2 = lrv_par - P + 1, lower.tail = FALSE)
  } else {
    stop("Invalid 'lrv' value. Use 'NeweyWest' or 'EWC'.")
  }
  
  return(list(S = S, W = W, pval = pval))
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
    S <- NeweyWest(Z_by_group, lrv_par = lrv_par)
    S_inv <- tryCatch(chol2inv(chol(S)), error = function(e) stop("S is not positive definite."))
    W <- as.numeric(Tobs * t(mu_hat) %*% S_inv %*% mu_hat)
    pval <- pchisq(W, df = KP, lower.tail = FALSE)
  } else if (lrv == "EWC") {
    EWC_res <- EWC(Z_by_group, lrv_par = lrv_par)
    S <- EWC_res$S
    lrv_par <- EWC_res$lrv_par
    adj_term <- ((lrv_par - KP + 1) / lrv_par) / KP
    S_inv <- tryCatch(chol2inv(chol(S)), error = function(e) stop("S is not positive definite."))
    W <- as.numeric(adj_term * Tobs * t(mu_hat) %*% S_inv %*% mu_hat)
    pval <- pf(W, df1 = KP, df2 = lrv_par - KP + 1, lower.tail = FALSE)
  } else {
    stop("Unknown lrv. Use 'EWC' or 'NeweyWest'.")
  }
  
  return(list(W = W, S = S, pval = pval))
}

#' Split-Sample Clustered EPA (matrix interface)
#'
#' @param Z NT x P matrix of loss differentials or moment conditions.
#' @param id Vector of length NT for unit identifiers.
#' @param time Vector of length NT for time identifiers.
#' @param K Integer; number of clusters for KMeans.
#' @param prop Numeric; proportion of time periods for R (training) sample.
#' @param lrv Character; "EWC" or "NeweyWest" for LRV estimator.
#' @param lrv_par Integer; lag (NW) or B (EWC).
#' @param Ninit Integer; number of KMeans initializations.
#' @param iter.max Integer; max iterations in KMeans.
#'
#' @return List with cluster assignments and test results from P
#' @export 
epa_clustered_split <- function(
  Z, id, time,
  K = NULL,
  Kmax = NULL,
  prop = 0.5,
  lrv = "EWC",
  lrv_par = NULL,
  Ninit = 10,
  iter.max = 10) {
                                  
  # Split matrices into R and P sets
  split <- split_panel_matrix(Z, id = id, time = time, prop = prop)
  Z_Tr <- split$Z_Tr
  Z_Te <- split$Z_Te
  id_Tr <- split$id_Tr
  id_Te <- split$id_Te
  time_Tr <- split$time_Tr
  time_Te <- split$time_Te

  # Run k-means clustering on R sample
  km_res <- panel_kmeans_estimation(Z = Z_Tr, id = id_Tr, time = time_Tr, K = K, Kmax = Kmax,
                                    Ninit = Ninit, iter.max = iter.max)
  gamma_Tr <- km_res$final_cluster

  # Map cluster assignments from R to P by id
  id_map <- setNames(gamma_Tr, unique(id_Tr))
  gamma_Te <- id_map[as.character(unique(id_Te))]
  if (any(is.na(gamma_Te))) stop("Some units in P were not found in R clustering.")

  # Run clustered EPA test on P using clusters from R
  test_res <- epa_clustered_known(
    Z = Z_Te,
    id = id_Te,
    time = time_Te,
    gamma = gamma_Te,
    lrv = lrv,
    lrv_par = lrv_par
  )
  
  return(append(list(clustering = gamma_Tr),test_res))
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
  pcombine_fun = "Genmean_neq",
  r = -20,
  lrv = "EWC",
  lrv_par = NULL,
  Ninit = 10,
  iter.max = 10
) {
  # Pairwise cluster p-values (with flexible combination)
  homo_res <- panel_homogeneity_test(
    Z = Z,
    id = id,
    time = time,
    K = K,
    Kmax = Kmax,
    lrv = lrv, lrv_par = lrv_par,
    pairs = pairs,
    pcombine_fun = pcombine_fun,
    r = r
  )

  pairwise_pvals <- homo_res$pairwise_pvalues
  pairs_out <- homo_res$pairs
  pval_combined = homo_res$pval

  # Overall EPA test p-value
  oepa <- overall_EPA_test(Z = Z, id = id, time = time, lrv = lrv, lrv_par = lrv_par)
  overall_pval <- oepa$pval

  # Combine all p-values (pairwise + overall)
  all_pvals <- c(pairwise_pvals, overall_pval)
  names(all_pvals) <- c(
    paste0("pair_", sapply(pairs_out, function(x) paste(x, collapse = "_"))),
    "overall"
  )

  # Use the same combination method as for pairwise, but now on all p-values
  selective_pval <- switch(
    pcombine_fun,
    Geomean     = Geomean_pcombine(all_pvals),
    Genmean     = Genmean_pcombine(all_pvals, r = r),
    Genmean_neq = Genmean_rneg_pcombine(all_pvals, r = r),
    iu          = iu_pcombine(all_pvals, r = r),
    bonferroni  = bonferroni_pcombine(all_pvals),
    cauchy      = cauchy_pcombine(all_pvals),
    stop("Invalid pcombine_fun specified.")
  )

  return(list(
    pairwise_pvalues = pairwise_pvals,
    pairs = pairs_out,
    overall_pval = overall_pval,
    hom_test_pval = pval_combined,
    pval = selective_pval,
    clustering = homo_res$clustering
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