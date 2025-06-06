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