#' Convert Balanced Panel Data to Matrix Forms for EPA Testing and KMeans
#'
#' @param df A balanced panel data frame
#' @param id Character; name of unit identifier column
#' @param time Character; name of time identifier column
#' @param Z_names Character vector of variable names to extract
#'
#' @return A list with:
#' \describe{
#'   \item{Zbar}{T x K matrix of time-wise cross-sectional averages}
#'   \item{Z_panel}{NT x K matrix for clustering (unit-time stacked)}
#'   \item{id_numeric}{Numeric vector of unit identifiers (1 to N)}
#'   \item{time_numeric}{Numeric vector of time identifiers (1 to T)}
#'   \item{id}{Original id vector, ordered}
#'   \item{time}{Original time vector, ordered}
#' }
#' @keywords internal
panel_data_to_matrix <- function(df, id, time, Z_names) {
  # Input checks
  if (!is.data.frame(df)) stop("df must be a data frame.")
  if (!all(c(id, time) %in% names(df))) stop("id and/or time columns not found in df.")
  if (!all(Z_names %in% names(df))) stop("Some Z_names not found in df.")
  
  # Order by id and time
  df <- df[order(df[[id]], df[[time]]), ]
  
  # Convert id and time to numeric (sequential, unique)
  id_levels <- unique(df[[id]])
  time_levels <- unique(df[[time]])
  df$id_numeric <- match(df[[id]], id_levels)
  df$time_numeric <- match(df[[time]], time_levels)
  
  # T x K: time-averaged cross-sectional means
  Zbar_df <- aggregate(df[, Z_names, drop = FALSE], by = list(df[[time]]), FUN = mean)
  Zbar <- as.matrix(Zbar_df[, -1, drop = FALSE])
  colnames(Zbar) <- Z_names
  
  # NT x K: stacked (unit-major) matrix
  Z_panel <- as.matrix(df[, Z_names, drop = FALSE])
  colnames(Z_panel) <- Z_names
  
  return(list(
    Zbar = Zbar,
    Z_panel = Z_panel,
    id_numeric = df$id_numeric,
    time_numeric = df$time_numeric,
    id = df[[id]],
    time = df[[time]]
  ))
}

#' Compute Group Means by Cluster Membership
#'
#' @description Computes column means of a numeric matrix grouped by a vector of categorical indicators.
#'
#' @param X A numeric matrix of size N x P.
#' @param v A vector of group membership indicators of length N.
#'
#' @return A K x P matrix of group means, where K is the number of unique values in \code{v}.
#' @keywords internal
aggregate_matrix <- function(X, v) {
  if (!is.matrix(X)) stop("X must be a matrix.")
  v <- as.factor(v)
  counts <- as.numeric(table(v))
  group_sums <- rowsum(X, v)
  Xbar <- sweep(group_sums, 1, counts, FUN = "/")
  return(Xbar)
}

# Handles the creation of R and P subsamples from a balanced panel dataset
#' Split a balanced panel matrix into estimation (R) and evaluation (P) samples
#'
#' @param Z Matrix of dimension (N*T) x K
#' @param id Numeric vector of length N*T (unit id)
#' @param time Numeric vector of length N*T (time id)
#' @param prop Proportion of time periods for R sample
#' @param maxlag Lags to remove from R sample
#'
#' @return List with matrices and aligned id/time: Z_R, Z_P, id_R, id_P, time_R, time_P, P, Tobs
#'
#' @keywords internal
#' 
split_panel_matrix <- function(Z, id, time, prop = 0.5, maxlag = 0) {
  df <- data.frame(Z)
  df$id <- id
  df$time <- time
  
  df <- df[order(df$id, df$time), ]
  Tobs <- length(unique(time))
  
  # Estimate and prediction splits
  library(dplyr)
  data_R <- df %>% dplyr::group_by(id) %>% dplyr::slice_head(prop = prop)
  data_P <- dplyr::anti_join(df, data_R, by = c("id", "time"))
  
  # Trim maxlag from R
  if (maxlag > 0) {
    data_R <- data_R %>% dplyr::group_by(id) %>% dplyr::slice(1:(dplyr::n() - maxlag))
  }
  
  # Extract matrices and aligned ids
  Z_R <- as.matrix(data_R[, 1:(ncol(Z))])
  Z_P <- as.matrix(data_P[, 1:(ncol(Z))])
  id_R <- data_R$id
  id_P <- data_P$id
  time_R <- data_R$time
  time_P <- data_P$time
  P <- length(unique(time_P))
  
  return(list(
    Z_R = Z_R,
    id_R = id_R,
    time_R = time_R,
    Z_P = Z_P,
    id_P = id_P,
    time_P = time_P,
    P = P,
    Tobs = Tobs
  ))
}