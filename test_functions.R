#' EPA Wald Test from Cross-Sectional Averages of Z
#'
#' @param Zbar T x K matrix of cross-sectional averages
#' @param method Character; "EWC" (default) or "NeweyWest"
#' @param param Integer; lag (NeweyWest) or B (EWC)
#'
#' @return List with elements: W (test stat), df (degrees of freedom), pval
#' @keywords internal
#' 
epa_overall_test <- function(Zbar, method = "EWC", param = NULL) {
  Tobs <- nrow(Zbar)
  K <- ncol(Zbar)
  
  mu_hat <- colMeans(Zbar)
  
  if (method == "NeweyWest") {
    if (is.null(param)) param <- floor(1.3 * sqrt(Tobs))
    Omega <- NeweyWest(Zbar, maxlag = param)
    W <- Tobs * t(mu_hat) %*% solve(Omega) %*% mu_hat
    pval <- pchisq(W, df = K, lower.tail = FALSE)
    df <- K
  } else if (method == "EWC") {
    if (is.null(param)) param <- Tobs
    Omega <- EWC(Zbar, B = param)
    adj <- ((param - K + 1) / param) / K
    W <- adj * Tobs * t(mu_hat) %*% solve(Omega) %*% mu_hat
    df <- c(K, param - K + 1)
    pval <- pf(W, df1 = df[1], df2 = df[2], lower.tail = FALSE)
  } else {
    stop("Unknown method: choose 'EWC' or 'NeweyWest'")
  }
  
  return(list(W = as.numeric(W), df = df, pval = pval))
}

#' EPA Test with Known Clusters
#'
#' @param Z T x (K * G) matrix of forecast error averages stacked by cluster
#' @param G Integer; number of clusters
#' @param K Integer; dimension of each Z_k
#' @param Tobs Integer; number of time periods
#' @param method Character; "EWC" (default) or "NeweyWest"
#' @param param Integer; lag or number of spectral base estimators
#'
#' @return A list with elements:
#'   \describe{
#'     \item{W}{Numeric; Wald statistic}
#'     \item{df}{Degrees of freedom (either scalar or vector)}
#'     \item{pval}{P-value}
#'   }
#'   
#' @examples
#' Z <- matrix(rnorm(200), nrow = 20)
#' epa_overall_test(Z, G, K, method = "EWC", param = 10)
#'   
#' @keywords internal
#' 
epa_clustered_known <- function(Z, G, K, method = "EWC", param = NULL) {
  Tobs <- nrow(Z)
  GK <- ncol(Z)
  
  mu_hat <- colMeans(Z)
  
  if (method == "NeweyWest") {
    if (is.null(param)) param <- floor(1.3 * sqrt(Tobs))
    Omega <- NeweyWest(Z, maxlag = param)
    W <- Tobs * t(mu_hat) %*% solve(Omega) %*% mu_hat
    df <- GK
    pval <- pchisq(W, df = df, lower.tail = FALSE)
  } else if (method == "EWC") {
    if (is.null(param)) param <- Tobs
    Omega <- EWC(Z, B = param)
    adj <- ((param - GK + 1) / param) / GK
    W <- adj * Tobs * t(mu_hat) %*% solve(Omega) %*% mu_hat
    df <- c(GK, param - GK + 1)
    pval <- pf(W, df1 = df[1], df2 = df[2], lower.tail = FALSE)
  } else {
    stop("Unknown method. Please use either 'EWC' or 'NeweyWest'.")
  }
  
  return(list(W = as.numeric(W), df = df, pval = pval))
}

#' EPA Split Sample Test Using Cluster Labels (Matrix Input)
#'
#' @param Z NTxK matrix of forecast errors
#' @param time NTx1 numeric vector of time identifiers
#' @param clusters NTx1 numeric vector of cluster assignments
#' @param method Character; "EWC" (default) or "NeweyWest"
#' @param param Integer; smoothing parameter for LR variance
#'
#' @return List with elements: W, df, pval
#' @keywords internal
epa_split_test_matrix <- function(Z, time, clusters, method = "EWC", param = NULL) {
  # Check input
  if (!is.matrix(Z)) stop("Z must be a matrix.")
  if (length(time) != nrow(Z)) stop("Length of 'time' must match rows of Z.")
  if (length(clusters) != nrow(Z)) stop("Length of 'clusters' must match rows of Z.")
  
  Tobs <- length(unique(time))
  K <- ncol(Z)
  G <- length(unique(clusters))
  
  # Combine time and cluster labels
  df <- data.frame(time = time, cluster = clusters, Z)
  colnames(df)[3:ncol(df)] <- paste0("V", 1:K)  # Rename Z columns for reshape2
  
  # Aggregate cross-sectional means by time and cluster
  agg <- stats::aggregate(df[, 3:(2 + K)], 
                          by = list(time = df$time, cluster = df$cluster),
                          FUN = mean)
  
  Z_wide <- reshape2::dcast(agg, time ~ cluster, value.var = paste0("V", 1:K))
  Zmat <- as.matrix(Z_wide[, -1])  # Drop time column
  
  return(epa_clustered_known(Zmat, G = G, K = K, Tobs = Tobs, method = method, param = param))
}

#' Split-Sample Clustered EPA Pipeline (User Function)
#'
#' @param df Balanced panel data.frame
#' @param id Unit identifier (character)
#' @param time Time identifier (character)
#' @param Z_names Outcome variable(s) used for clustering (character vector)
#' @param G Number of clusters for KMeans
#' @param prop Proportion of data used in R (training) sample
#' @param method "EWC" or "NeweyWest"
#' @param param Integer; lag (NW) or B (EWC)
#' @param Ninit Number of KMeans initializations
#' @param iter.max Max iterations in KMeans
#'
#' @return List with cluster assignments and test results from P
#' @export
epa_split_pipeline <- function(df, id, time, Z_names,
                               G, prop = 0.5,
                               method = "EWC", param = NULL,
                               Ninit = 10, iter.max = 10) {
  # Step 1: Convert to matrices
  matrices <- panel_data_to_matrix(df, id = id, time = time, Z_names = Z_names)
  Z <- matrices$Z
  id_vec <- matrices$id
  time_vec <- matrices$time
  
  # Step 2: Split matrices into R and P sets
  split <- split_panel_matrix(Z, id = id_vec, time = time_vec, prop = prop)
  Z_R <- split$Z_R
  Z_P <- split$Z_P
  id_R <- split$id_R
  id_P <- split$id_P
  time_P <- split$time_P
  
  # Step 3: Run k-means clustering on R sample
  km_res <- panel_kmeans_estimation_matrix(Z = Z_R, id = id_R, G = G,
                                           Ninit = Ninit, iter.max = iter.max)
  cluster_assignments <- km_res$final_cluster
  
  # Step 4: Map cluster assignments from R to P by id
  id_map <- setNames(cluster_assignments, unique(id_R))
  cluster_P <- id_map[as.character(id_P)]
  
  if (any(is.na(cluster_P))) stop("Some units in P were not found in R clustering.")
  
  # Step 5: Run EPA test on P using clusters from R
  test_res <- epa_split_test_matrix(Z = Z_P, time = time_P,
                                    clusters = cluster_P,
                                    method = method, param = param)
  
  return(list(
    clustering = cluster_assignments,
    epa_test = test_res
  ))
}