#' Generate Forecast Simulation Data with Cluster-Specific Forecast Variance Deviations
#'
#' @param N Integer; number of units.
#' @param Tobs Integer; number of time periods (after burn-in).
#' @param burn Integer; burn-in length.
#' @param alpha Intercept for AR(1) process (default = 1).
#' @param rho_vec Numeric vector of length 3; AR coefficients for 3 clusters.
#' @param gamma Integer vector of length N; cluster membership in {1,2,3}.
#' @param phi AR(1) coefficient for V_{it} innovations.
#' @param lambda Cross-sectional loading of the common shock.
#' @param delta_vec Numeric vector of length 3; deviations from optimal forecast variance.
#'
#' @return A list with Y, forecast_1, and forecast_2 matrices (N x (Tobs - 1)).
#' @export
generate_forecast_simulation_data <- function(N, Tobs, burn,
                                              alpha = 1,
                                              rho_vec = c(0.2, 0.3, 0.4),
                                              gamma,
                                              phi = 0.2,
                                              lambda = 0.2,
                                              delta_vec = c(0, -0.1, -0.2)) {
  if (length(rho_vec) != 3) stop("rho_vec must be a vector of length 3.")
  if (length(gamma) != N) stop("gamma must be a vector of length N with values in {1,2,3}.")
  if (length(delta_vec) != 3) stop("delta_vec must be a vector of length 3.")
  
  T_total <- burn + Tobs
  
  rho_i <- rho_vec[gamma]
  delta_i <- delta_vec[gamma]
  print(delta_i)
  sigma2_V_i <- alpha^2 * (1 - rho_i)^2 + delta_i
  if (any(sigma2_V_i <= 0)) stop("Resulting forecast variance is non-positive. Adjust delta_vec.")
  
  sigma2_xi_i <- sigma2_V_i * (1 - phi^2) - lambda^2
  if (any(sigma2_xi_i < 0)) stop("Choose smaller lambda or phi; some cluster variances invalid.")
  
  # === DGP ===
  Y <- matrix(0, nrow = N, ncol = T_total)
  U <- matrix(rnorm(N * T_total), N, T_total)
  
  for (t in 2:T_total) {
    Y[, t] <- alpha * (1 - rho_i) + rho_i * Y[, t - 1] + U[, t]
  }
  
  Y <- Y[, (burn + 1):T_total]  # Remove burn-in
  Y_lag <- Y[, 1:(Tobs - 1)]
  Y_outcome <- Y[, 2:Tobs]
  
  # === V_it construction ===
  F_common <- rnorm(burn + Tobs - 1)
  V <- matrix(0, nrow = N, ncol = burn + Tobs - 1)
  V[, 1] <- rnorm(N, mean = 0, sd = sqrt(sigma2_V_i))
  
  for (t in 2:(burn + Tobs - 1)) {
    xi_t <- rnorm(N, 0, sqrt(sigma2_xi_i))
    V[, t] <- phi * V[, t - 1] + lambda * F_common[t] + xi_t
  }
  
  V <- V[, (burn + 1):(burn + Tobs - 1)]
  
  # === Forecasts ===
  forecast_1 <- alpha * (1 - rho_i) + rho_i * Y_lag + V
  forecast_2 <- rho_i * Y_lag
  
  return(list(
    Y = Y_outcome,
    forecast_1 = forecast_1,
    forecast_2 = forecast_2
  ))
}