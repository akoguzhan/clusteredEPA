#' Long-run Variance Estimation via Eigenvalue Weighted Covariance (EWC)
#'
#' @description Estimates the long-run variance using the Eigenvalue Weighted Covariance estimator.
#'
#' @param X A TxK numeric matrix.
#' @param lrv_par An integer; number of base estimators (e.g., harmonic terms). If NULL, uses floor(Tobs^0.5).
#'
#' @return A list with the KxK estimated long-run variance matrix and the lrv_par used.
#' @keywords internal
#' 
EWC <- function(X, lrv_par = NULL) {
  X <- as.matrix(X)
  Tobs <- NROW(X)
  P <- NCOL(X)
  
  # Rule of thumb: sqrt(Tobs) harmonics if not specified
  if (is.null(lrv_par)) lrv_par <- floor(Tobs^0.5)
  if (lrv_par < 1 || lrv_par > Tobs) stop("lrv_par must be between 1 and Tobs.")
  
  if (lrv_par == Tobs) {
    OmegaHat <- NeweyWest(X, 0)
  } else {
    if (P == 1) {
      X <- X - mean(X)
    } else {
      X <- scale(X, scale = FALSE)
    }
    OmegaHat <- 0
    for (j in 1:lrv_par) {
      Cj <- cos(pi * j * (seq(1:Tobs) - 1/2) / Tobs)
      Lj <- sqrt(2 / Tobs) * t(X) %*% Cj
      OmegaHat <- OmegaHat + (Lj %*% t(Lj)) / lrv_par
    }
  }
  return(list(S = OmegaHat, lrv_par = lrv_par))
}

#' Long-run Variance Estimation via Newey-West
#'
#' @description Estimates the long-run variance of a time series using the Newey-West estimator.
#'
#' @param X A TxK numeric matrix.
#' @param lrv_par An integer specifying the maximum lag order.
#'
#' @return A KxK estimated long-run variance matrix.
#' @keywords internal
#' 
NeweyWest <- function(X, lrv_par = NULL) {
  X <- as.matrix(X)
  Tobs <- NROW(X)
  P <- NCOL(X)
  
  if (Tobs < 2) stop("Not enough observations.")
  if (is.null(lrv_par)) lrv_par <- floor(4 * (Tobs / 100)^(2/9))
  if (lrv_par >= Tobs) stop("lrv_par must be less than the number of observations.")
  
  # Demean
  if (P == 1) {
    X <- X - mean(X)
  } else {
    X <- scale(X, scale = FALSE)
  }
  
  samplevar <- t(X) %*% X / Tobs
  OmegaHat <- samplevar
  if (lrv_par > 0) {
    for (h in 1:lrv_par) {
      if (P == 1) {
        Xlag <- c(rep(0, h), X[1:(Tobs - h), ])
      } else {
        Xlag <- rbind(matrix(0, h, P), X[1:(Tobs - h), ])
      }
      gamma <- (t(X) %*% Xlag + t(Xlag) %*% X) / Tobs
      weight <- 1 - (h / (lrv_par + 1))
      OmegaHat <- OmegaHat + weight * gamma
    }
  }
  return(OmegaHat)
}