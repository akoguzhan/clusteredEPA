#' One-Step-Ahead Forecast from ARX(p,q) Model
#'
#' @description Computes a one-step-ahead forecast from an ARX(p,q) model estimated via OLS.
#'
#' @param Y A numeric vector of the dependent variable.
#' @param X A numeric matrix of the same length as `Y`, with K columns of predictors.
#' @param p Integer; number of lags for Y.
#' @param q Integer; number of lags for X.
#'
#' @return A scalar numeric value giving the one-step-ahead forecast.
#' @export
forecast_arx <- function(Y, X, p, q) {
  Tobs <- length(Y)
  if (!is.matrix(X)) X <- matrix(X, ncol = 1)
  K <- ncol(X)
  max_lag <- max(p, q)
  start <- max_lag + 1
  
  n_reg <- Tobs - start
  Xmat <- matrix(NA, nrow = n_reg, ncol = 1 + p + q * K)
  Yvec <- rep(NA, n_reg)
  
  for (t in 1:n_reg) {
    tt <- start + t - 1
    y_lags <- rev(Y[(tt - p):(tt - 1)])
    x_lags <- as.vector(sapply(1:q, function(l) X[tt - l, ]))
    Xmat[t, ] <- c(1, y_lags, x_lags)
    Yvec[t] <- Y[tt]
  }
  
  beta_hat <- solve(t(Xmat) %*% Xmat, t(Xmat) %*% Yvec)
  
  y_lags_new <- rev(Y[(Tobs - p + 1):Tobs])
  x_lags_new <- as.vector(sapply(1:q, function(l) X[Tobs - l + 1, ]))
  Xnew <- c(1, y_lags_new, x_lags_new)
  
  forecast <- sum(Xnew * beta_hat)
  return(forecast)
}

#' Pooled Panel ARX(p,q) Forecasts from Stacked Matrices
#'
#' @description Computes one-step-ahead forecasts using a pooled ARX(p,q) model on panel data
#' where inputs are NT x 1 and NT x P stacked matrices. Homogeneous coefficients are assumed.
#'
#' @param Y NT x 1 matrix (or vector) of the dependent variable stacked by unit over time.
#' @param X NT x P matrix of exogenous regressors stacked the same way.
#' @param id Vector of unit identifiers (length NT).
#' @param time Vector of time indices (length NT).
#' @param p Integer; number of lags of Y.
#' @param q Integer; number of lags of X.
#'
#' @return A named numeric vector of length N (number of units) containing the one-step-ahead forecasts.
#' @export
forecast_panel_arx <- function(Y, X, id, time, p, q) {
  if (!is.matrix(Y)) Y <- as.matrix(Y)
  if (!is.matrix(X)) stop("X must be a matrix.")
  
  NT <- length(Y)
  N <- length(unique(id))
  Tobs <- length(unique(time))
  max_lag <- max(p, q)
  P <- ncol(X)
  
  # Combine and sort
  df <- data.frame(Y = Y, X, id = id, time = time)
  df <- df[order(df$id, df$time), ]
  Y_sorted <- as.matrix(df$Y)
  X_sorted <- as.matrix(df[, 2:(1 + P)])
  id_sorted <- df$id
  time_sorted <- df$time
  
  # Construct regression design matrix
  reg_data <- list()
  y_target <- c()
  
  for (i in unique(id_sorted)) {
    idx_i <- which(id_sorted == i)
    y_i <- Y_sorted[idx_i]
    x_i <- X_sorted[idx_i, , drop = FALSE]
    
    for (t in (max_lag + 1):(Tobs - 1)) {
      y_lags <- rev(y_i[(t - p):(t - 1)])
      x_lags <- as.vector(sapply(1:q, function(l) x_i[t - l, ]))
      reg_data[[length(reg_data) + 1]] <- c(1, y_lags, x_lags)
      y_target <- c(y_target, y_i[t])
    }
  }

  Xmat <- do.call(rbind, reg_data)
  beta_hat <- solve(t(Xmat) %*% Xmat, t(Xmat) %*% y_target)

  # Generate forecasts for next period (Tobs) for each unit
  forecasts <- rep(NA, N)
  names(forecasts) <- as.character(unique(id_sorted))
  
  for (i in unique(id_sorted)) {
    idx_i <- which(id_sorted == i)
    y_i <- Y_sorted[idx_i]
    x_i <- X_sorted[idx_i, , drop = FALSE]
    
    y_lags <- rev(y_i[(Tobs - p + 1):Tobs])
    x_lags <- as.vector(sapply(1:q, function(l) x_i[Tobs - l + 1, ]))
    Xnew <- c(1, y_lags, x_lags)
    forecasts[as.character(i)] <- sum(Xnew * beta_hat)
  }
  
  return(forecasts)
}

#' Fixed Effects Panel ARX(p,q) Forecasts with Estimated Intercepts
#'
#' @description Computes one-step-ahead forecasts using a fixed effects ARX(p,q) model on panel data.
#' Forecasts include estimated unit-specific intercepts from the regression.
#'
#' @param Y NT x 1 matrix (or vector) of the dependent variable stacked by unit over time.
#' @param X NT x P matrix of exogenous regressors stacked the same way.
#' @param id Vector of unit identifiers (length NT).
#' @param time Vector of time indices (length NT).
#' @param p Integer; number of lags of Y.
#' @param q Integer; number of lags of X.
#'
#' @return A named numeric vector of length N (number of units) containing the one-step-ahead forecasts.
#' @export
forecast_panel_arx_fe <- function(Y, X, id, time, p, q) {
  if (!is.matrix(Y)) Y <- as.matrix(Y)
  if (!is.matrix(X)) stop("X must be a matrix.")
  
  NT <- length(Y)
  N <- length(unique(id))
  Tobs <- length(unique(time))
  max_lag <- max(p, q)
  P <- ncol(X)
  
  # Combine and sort
  df <- data.frame(Y = Y, X, id = id, time = time)
  df <- df[order(df$id, df$time), ]
  Y_sorted <- as.matrix(df$Y)
  X_sorted <- as.matrix(df[, 2:(1 + P)])
  id_sorted <- df$id
  time_sorted <- df$time
  
  # Build regression design and target vectors
  reg_data <- list()
  y_target <- c()
  id_index <- c()
  
  for (i in unique(id_sorted)) {
    idx_i <- which(id_sorted == i)
    y_i <- Y_sorted[idx_i]
    x_i <- X_sorted[idx_i, , drop = FALSE]
    
    for (t in (max_lag + 1):(Tobs - 1)) {
      y_lags <- rev(y_i[(t - p):(t - 1)])
      x_lags <- as.vector(sapply(1:q, function(l) x_i[t - l, ]))
      reg_data[[length(reg_data) + 1]] <- c(y_lags, x_lags)
      y_target <- c(y_target, y_i[t])
      id_index <- c(id_index, i)
    }
  }
  
  Xmat <- do.call(rbind, reg_data)
  id_index <- factor(id_index)
  Dmat <- model.matrix(~ id_index - 1)  # Fixed effects dummies
  X_reg <- cbind(Dmat, Xmat)
  
  # Estimate full model with fixed effects
  beta_full <- solve(t(X_reg) %*% X_reg, t(X_reg) %*% y_target)
  alpha_hat <- beta_full[1:N]                       # Fixed effects
  beta_hat <- beta_full[(N + 1):length(beta_full)]  # Slope coefficients
  
  # Forecasts
  forecasts <- rep(NA, N)
  names(forecasts) <- as.character(levels(id_index))
  
  for (i in unique(id_sorted)) {
    idx_i <- which(id_sorted == i)
    y_i <- Y_sorted[idx_i]
    x_i <- X_sorted[idx_i, , drop = FALSE]
    
    y_lags <- rev(y_i[(Tobs - p + 1):Tobs])
    x_lags <- as.vector(sapply(1:q, function(l) x_i[Tobs - l + 1, ]))
    Xnew <- c(y_lags, x_lags)
    
    forecasts[as.character(i)] <- alpha_hat[i] + sum(Xnew * beta_hat)
  }
  
  return(forecasts)
}

#' Unit-Level ARX(p,q) Forecasts for Panel Data (Fully Heterogeneous)
#'
#' @description Computes one-step-ahead forecasts for panel data by estimating separate ARX(p,q)
#' models for each unit using OLS. Each unit has its own model and coefficients.
#'
#' @param Y A numeric vector of length NT (stacked panel data of the dependent variable).
#' @param X A numeric matrix of size NT x P (stacked panel data of regressors).
#' @param id A vector of unit identifiers of length NT.
#' @param time A vector of time identifiers of length NT.
#' @param p Number of lags for Y.
#' @param q Number of lags for X.
#'
#' @return A named numeric vector of length N containing one-step-ahead forecasts per unit.
#' @export
forecast_arx_heterogeneous <- function(Y, X, id, time, p, q) {
  if (!is.matrix(X)) stop("X must be a matrix.")
  if (!is.vector(Y)) Y <- as.vector(Y)
  
  NT <- length(Y)
  P <- ncol(X)
  max_lag <- max(p, q)
  
  # Sort by (id, time)
  df <- data.frame(Y = Y, X, id = id, time = time)
  df <- df[order(df$id, df$time), ]
  
  Y_sorted <- df$Y
  X_sorted <- as.matrix(df[, 2:(1 + P)])
  id_sorted <- df$id
  time_sorted <- df$time
  units <- unique(id_sorted)
  Tobs <- length(unique(time_sorted))
  N <- length(units)
  
  forecasts <- rep(NA, N)
  names(forecasts) <- as.character(units)
  
  for (i in units) {
    idx_i <- which(id_sorted == i)
    y_i <- Y_sorted[idx_i]
    x_i <- X_sorted[idx_i, , drop = FALSE]
    
    T_i <- length(y_i)
    n_reg <- T_i - max_lag - 1
    if (n_reg <= 0) next
    
    Xmat <- matrix(NA, nrow = n_reg, ncol = 1 + p + q * P)
    Yvec <- rep(NA, n_reg)
    
    for (t in 1:n_reg) {
      tt <- max_lag + t
      y_lags <- rev(y_i[(tt - p):(tt - 1)])
      x_lags <- as.vector(sapply(1:q, function(l) x_i[tt - l, ]))
      Xmat[t, ] <- c(1, y_lags, x_lags)
      Yvec[t] <- y_i[tt]
    }
    
    beta_hat <- tryCatch(
      solve(t(Xmat) %*% Xmat, t(Xmat) %*% Yvec),
      error = function(e) rep(NA, ncol(Xmat))
    )
    
    # Forecast for time Tobs + 1
    y_lags_new <- rev(y_i[(T_i - p + 1):T_i])
    x_lags_new <- as.vector(sapply(1:q, function(l) x_i[T_i - l + 1, ]))
    Xnew <- c(1, y_lags_new, x_lags_new)
    
    forecasts[as.character(i)] <- sum(Xnew * beta_hat)
  }
  
  return(forecasts)
}
