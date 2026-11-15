# Handles the creation of R and P subsamples from a balanced panel dataset
#' Split a balanced panel matrix into estimation (R) and evaluation (P) samples
#'
#' @description Sorts the panel by unit and time, assigns the first \code{prop} share of each unit's
#' periods to the training sample and the rest to the testing sample, and drops the first
#' \code{maxlag} periods of each unit from the testing sample.
#'
#' @param Z Matrix of dimension (N*T) x K
#' @param id Numeric vector of length N*T (unit id)
#' @param time Numeric vector of length N*T (time id)
#' @param prop Proportion of time periods for the training sample
#' @param maxlag Number of initial periods removed from the testing sample (default
#'   \code{floor(sqrt(prop * T))})
#'
#' @return List with matrices and aligned id/time: Z_Tr, Z_Te, id_Tr, id_Te, time_Tr, time_Te,
#'   P (number of testing periods), Tobs
#'
#' @keywords internal
#' 
split_panel_matrix <- function(Z, id, time, prop = 0.5, maxlag = NULL) {
  df <- data.frame(Z)
  df$id <- id
  df$time <- time
  
  df <- df[order(df$id, df$time), ]
  Tobs <- length(unique(time))
  
  # Estimate and prediction splits
  library(dplyr)
  data_Tr <- df %>% dplyr::group_by(id) %>% dplyr::slice_head(prop = prop)
  data_Te <- dplyr::anti_join(df, data_Tr, by = c("id", "time"))
  
  # Trim maxlag from R
  if (is.null(maxlag)) {
	maxlag <- floor(sqrt(prop * Tobs))
  }
  if (maxlag > 0) {
    data_Te <- data_Te %>% dplyr::group_by(id) %>% dplyr::slice((maxlag + 1):(dplyr::n()))
  }
  
  # Extract matrices and aligned ids
  Z_Tr <- as.matrix(data_Tr[, 1:(ncol(Z))])
  Z_Te <- as.matrix(data_Te[, 1:(ncol(Z))])
  id_Tr <- data_Tr$id
  id_Te <- data_Te$id
  time_Tr <- data_Tr$time
  time_Te <- data_Te$time
  P <- length(unique(time_Te))
  
  return(list(
    Z_Tr = Z_Tr,
    id_Tr = id_Tr,
    time_Tr = time_Tr,
    Z_Te = Z_Te,
    id_Te = id_Te,
    time_Te = time_Te,
    P = P,
    Tobs = Tobs
  ))
}