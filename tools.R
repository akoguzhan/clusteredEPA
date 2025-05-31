# clusteredEPA_utils.R

#' Calculate adjustment term for EPA test
#'
#' @param B Number of blocks
#' @param K Number of variables
#'
#' @return Adjustment term (scalar)
adjustment_term <- function(B, K) {
  ((B - K + 1) / B) / K
}

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
#' }
#' @keywords internal
panel_data_to_matrix <- function(df, id, time, Z_names) {
  df <- df[order(df[[id]], df[[time]]), ]
  
  # Convert id and time to numeric
  df$id_numeric <- as.numeric(factor(df[[id]]))
  df$time_numeric <- as.numeric(factor(df[[time]]))
  
  # T x K: time-averaged cross-sectional means
  Zbar_df <- aggregate(df[, Z_names], by = list(df[[time]]), FUN = mean)
  Zbar <- as.matrix(Zbar_df[, -1])
  
  # NT x K: stacked (unit-major) matrix
  Z_panel <- as.matrix(df[, Z_names])
  
  return(list(
    Zbar = Zbar,
    Z_panel = Z_panel,
    id_numeric = df$id_numeric,
    time_numeric = df$time_numeric
  ))
}
