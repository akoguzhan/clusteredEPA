#' Generate Forecast Simulation Data (Original Design, True Parameters)
#'
#' @description Design of the original Monte Carlo study: the forecasters use the true AR(1)
#' parameters; Forecaster 1 adds a noise with cluster-specific variance.
#'
#' @param N Integer; number of units.
#' @param Tobs Integer; number of time periods (after burn-in).
#' @param burn Integer; burn-in length.
#' @param alpha Means of the AR(1) process by cluster (length 3, or a scalar recycled; default 1).
#' @param rho_vec Numeric vector of length 3; AR coefficients for 3 clusters.
#' @param gamma Integer vector of length N; cluster membership in 1, 2, 3.
#' @param phi AR(1) coefficient of the forecast noise \eqn{V_{it}} of Forecaster 1.
#' @param lambda Cross-sectional loading of the common shock.
#' @param delta_vec Numeric vector of length 3; deviations from optimal forecast variance.
#'
#' @return A list with Y, forecast_1, and forecast_2 matrices (N x (Tobs - 1)).
#' @export
generate_forecast_simulation_data <- function(N, Tobs, burn,
                                              alpha = c(1, 1, 1),
                                              rho_vec = c(0.2, 0.3, 0.4),
                                              gamma,
                                              phi = 0.2,
                                              lambda = 0.2,
                                              delta_vec = c(0, -0.1, -0.2)) {
  if (length(alpha) == 1) alpha <- rep(alpha, 3)
  if (length(alpha) != 3) stop("alpha must be a vector of length 3.")
  if (length(rho_vec) != 3) stop("rho_vec must be a vector of length 3.")
  if (length(gamma) != N) stop("gamma must be a vector of length N with values in {1,2,3}.")
  if (length(delta_vec) != 3) stop("delta_vec must be a vector of length 3.")

  T_total <- burn + Tobs
  alpha_i <- alpha[gamma]
  rho_i <- rho_vec[gamma]
  delta_i <- delta_vec[gamma]

  sigma2_V_i <- alpha_i^2 * (1 - rho_i)^2 + delta_i
  if (any(sigma2_V_i <= 0)) stop("Resulting forecast variance is non-positive. Adjust delta_vec.")

  sigma2_xi_i <- sigma2_V_i * (1 - phi^2) - lambda^2
  if (any(sigma2_xi_i < 0)) stop("Choose smaller lambda or phi; some cluster variances invalid.")

  Y <- matrix(0, nrow = N, ncol = T_total)
  U <- matrix(rnorm(N * T_total), N, T_total)
  for (t in 2:T_total) {
    Y[, t] <- alpha_i * (1 - rho_i) + rho_i * Y[, t - 1] + U[, t]
  }
  Y <- Y[, (burn + 1):T_total]
  Y_lag <- Y[, 1:(Tobs - 1)]
  Y_outcome <- Y[, 2:Tobs]

  F_common <- rnorm(burn + Tobs - 1)
  V <- matrix(0, nrow = N, ncol = burn + Tobs - 1)
  V[, 1] <- rnorm(N, mean = 0, sd = sqrt(sigma2_V_i))
  for (t in 2:(burn + Tobs - 1)) {
    xi_t <- rnorm(N, 0, sqrt(sigma2_xi_i))
    V[, t] <- phi * V[, t - 1] + lambda * F_common[t] + xi_t
  }
  V <- V[, (burn + 1):(burn + Tobs - 1)]

  forecast_1 <- alpha_i * (1 - rho_i) + rho_i * Y_lag + V
  forecast_2 <- rho_i * Y_lag

  return(list(Y = Y_outcome, forecast_1 = forecast_1, forecast_2 = forecast_2))
}

# ---- Extended design (Section 5 of the revised paper) ------------------------

# Simulate the target: Y_it = alpha_i (1 - rho_i) + rho_i Y_i,t-1 + g_i F_t + U_it,
# F_t AR(1) with coefficient fac_ar. Returns an N x n matrix (after burn-in).
# Optional C++ versions of .sim_target and .rolling_errors (inst/cpp/gen.cpp), compiled once with
# Rcpp and cached; the R versions below are used if Rcpp or a compiler is not available.
.gen_env <- new.env()
#' Compiled Simulation Routine by Name
#'
#' @description Compiles \code{inst/cpp/gen.cpp} with Rcpp once per session (cached on disk) and
#' returns the routine \code{name}. Returns \code{NULL} if the option \code{clusteredEPA.use_cpp}
#' is FALSE or if Rcpp or a compiler is not available.
#'
#' @param name Name of the compiled routine (e.g. \code{"ar_target_cpp"}).
#'
#' @return A function, or \code{NULL}.
#' @keywords internal
.get_cpp_gen <- function(name) {
  if (!isTRUE(getOption("clusteredEPA.use_cpp", TRUE))) return(NULL)
  if (is.null(.gen_env$e)) {
    .gen_env$e <- tryCatch({
      src <- system.file("cpp", "gen.cpp", package = "clusteredEPA")
      if (!nzchar(src) || !requireNamespace("Rcpp", quietly = TRUE)) stop("Rcpp not available")
      cache <- file.path(tools::R_user_dir("clusteredEPA", "cache"), "cpp")
      dir.create(cache, recursive = TRUE, showWarnings = FALSE)
      e <- new.env()
      ok <- tryCatch({
        suppressWarnings(suppressMessages(Rcpp::sourceCpp(src, env = e, cacheDir = cache, rebuild = FALSE, showOutput = FALSE)))
        TRUE
      }, error = function(err) FALSE)
      if (!ok) {   # e.g. concurrent compilation by several processes: compile privately
        e <- new.env()
        suppressWarnings(suppressMessages(Rcpp::sourceCpp(src, env = e, cacheDir = tempfile("cepa_cpp"), showOutput = FALSE)))
      }
      e
    }, error = function(err) {
      message("clusteredEPA: C++ simulation helpers not available (", conditionMessage(err),
              "); using the slower R versions.")
      FALSE
    })
  }
  if (isFALSE(.gen_env$e)) NULL else .gen_env$e[[name]]
}

#' Simulate the Target Variable of the Extended Design
#'
#' @description Simulates \eqn{Y_{it} = \alpha_i(1-\rho_i) + \rho_iY_{i,t-1} + g_iF_t + U_{it}}, with
#' \eqn{F_t} an AR(1) factor with coefficient \code{fac_ar} and standard normal innovations, starting
#' from \eqn{Y_{i1} = \alpha_i}, and drops the first \code{burn} periods. Uses the compiled
#' \code{ar_target_cpp} if available.
#'
#' @param n Number of periods returned.
#' @param alpha_i,rho_i,load_i Unit means, AR coefficients and factor loadings (length N).
#' @param fac_ar AR coefficient of the factor.
#' @param burn Burn-in length.
#'
#' @return An N x n matrix.
#' @keywords internal
.sim_target <- function(n, alpha_i, rho_i, load_i, fac_ar, burn) {
  N <- length(alpha_i)
  n_tot <- n + burn
  Fy <- as.numeric(stats::filter(rnorm(n_tot), fac_ar, method = "recursive"))
  Y <- matrix(0, N, n_tot)
  Y[, 1] <- alpha_i
  U <- matrix(rnorm(N * n_tot), N, n_tot)
  f <- .get_cpp_gen("ar_target_cpp")
  if (!is.null(f)) return(f(alpha_i, rho_i, load_i, Fy, U)[, (burn + 1):n_tot, drop = FALSE])
  for (t in 2:n_tot) Y[, t] <- alpha_i * (1 - rho_i) + rho_i * Y[, t - 1] + load_i * Fy[t] + U[, t]
  Y[, (burn + 1):n_tot, drop = FALSE]
}

# Row-wise rolling sums over windows of length w ending at column t (columns w..n)
#' Row-Wise Rolling Sums
#'
#' @description Computes, for each row of \code{M}, the sums over windows of \code{w} consecutive
#' columns ending at columns w, ..., n.
#'
#' @param M Numeric matrix with n columns.
#' @param w Window length.
#'
#' @return A matrix with \code{nrow(M)} rows and n - w + 1 columns.
#' @keywords internal
.roll_sum <- function(M, w) {
  cs <- cbind(0, t(apply(M, 1, cumsum)))
  n <- ncol(M)
  cs[, (w + 1):(n + 1), drop = FALSE] - cs[, 1:(n - w + 1), drop = FALSE]
}

# Rolling-window OLS forecasts of Y_t from AR(1) with and without intercept, window of R
# observations Y_{t-R},...,Y_{t-1} (R - 1 regression pairs). Returns errors for t = R+1..n
# and the lagged target.
#' Rolling-Window AR(1) Forecast Errors
#'
#' @description For each unit, estimates by least squares an AR(1) with intercept (Forecaster 1)
#' and without intercept (Forecaster 2) on rolling windows of R observations and computes the
#' one-step-ahead forecast errors. Uses the compiled \code{rolling_errors_cpp} if available.
#'
#' @param Y N x n matrix of the target (rows are units).
#' @param R Rolling window length.
#'
#' @return A list of N x (n - R) matrices: \code{e1}, \code{e2} (forecast errors), \code{ylag}
#'   (lagged target) and \code{y} (target).
#' @keywords internal
.rolling_errors <- function(Y, R) {
  f <- .get_cpp_gen("rolling_errors_cpp")
  if (!is.null(f)) return(f(Y, as.integer(R)))
  n <- ncol(Y)
  a <- Y[, 1:(n - 1), drop = FALSE]   # y_{s-1}
  b <- Y[, 2:n, drop = FALSE]         # y_s
  w <- R - 1
  Sa <- .roll_sum(a, w); Sb <- .roll_sum(b, w)
  Saa <- .roll_sum(a^2, w); Sab <- .roll_sum(a * b, w)
  # column j of the rolling sums uses pairs s = j+1..j+w, i.e. observations y_j..y_{j+w};
  # the forecast target is y_{j+w+1} = y_{j+R}, with lag y_{j+R-1}
  m <- ncol(Sa) - 1                    # number of forecasts (last window has no target)
  Sa <- Sa[, 1:m, drop = FALSE]; Sb <- Sb[, 1:m, drop = FALSE]
  Saa <- Saa[, 1:m, drop = FALSE]; Sab <- Sab[, 1:m, drop = FALSE]
  ylag <- Y[, R:(R + m - 1), drop = FALSE]
  y <- Y[, (R + 1):(R + m), drop = FALSE]
  beta1 <- (w * Sab - Sa * Sb) / (w * Saa - Sa^2)
  c1 <- (Sb - beta1 * Sa) / w
  beta2 <- Sab / Saa
  list(e1 = y - (c1 + beta1 * ylag), e2 = y - beta2 * ylag, ylag = ylag, y = y)
}

# Unit-variance AR(1) noise with a common factor: e_t = phi e_{t-1} + lambda F_t + sqrt(1-phi^2-lambda^2) xi_t
#' Simulate Unit-Variance AR(1) Noise with a Common Factor
#'
#' @description Simulates \eqn{e_{it} = \phi e_{i,t-1} + \lambda F_t + \sqrt{1-\phi^2-\lambda^2}\,\xi_{it}}
#' with standard normal \eqn{F_t} and \eqn{\xi_{it}}, and drops the first \code{burn} periods.
#'
#' @param N Number of units.
#' @param n Number of periods returned.
#' @param phi AR coefficient.
#' @param lambda Loading of the common factor.
#' @param burn Burn-in length (default 50).
#'
#' @return An N x n matrix.
#' @keywords internal
.sim_noise <- function(N, n, phi, lambda, burn = 50) {
  n_tot <- n + burn
  Fc <- rnorm(n_tot)
  E <- matrix(0, N, n_tot)
  E[, 1] <- rnorm(N)
  s <- sqrt(1 - phi^2 - lambda^2)
  for (t in 2:n_tot) E[, t] <- phi * E[, t - 1] + lambda * Fc[t] + s * rnorm(N)
  E[, (burn + 1):n_tot, drop = FALSE]
}

#' Calibrate the Forecast Noise of the Extended Design
#'
#' @description For each unit, simulates a long path (length \code{Tcal}) with the unit's
#' parameters and computes the moments of the loss differential without noise,
#' \eqn{d_t = e_{1t}^2 - e_{2t}^2}. The net noise variance \eqn{v_1(y)-v_2(y) = c_i + s_i w_i(y)}, with
#' \eqn{w_i(y) = \Phi((y-\alpha_i)/sd_i)}, solves
#' \eqn{E(\Delta L_{it}) = \psi_i} and \eqn{E(Y_{i,t-1}\Delta L_{it}) = \alpha_i\psi_i}
#' (\code{"both"}), or only the first equation with \eqn{s_i = 0} (\code{"unconditional"}).
#' Positive parts are added to Forecaster 1, negative parts to Forecaster 2.
#'
#' @param psi_i Target expected loss differential of each unit (length N).
#' @param alpha_i,rho_i,load_i Unit means, AR coefficients and factor loadings (length N).
#' @param fac_ar AR coefficient of the target factor.
#' @param R Rolling window length.
#' @param Tcal Length of the calibration path.
#' @param calibrate \code{"both"} or \code{"unconditional"}.
#' @param burn Burn-in length.
#'
#' @return A list with the constants \code{c}, the slopes \code{s} and the standard deviations
#'   \code{sd} of the lagged target (each of length N).
#' @keywords internal
.calibrate_noise <- function(psi_i, alpha_i, rho_i, load_i, fac_ar, R, Tcal, calibrate, burn) {
  Yc <- .sim_target(Tcal + R, alpha_i, rho_i, load_i, fac_ar, burn)
  er <- .rolling_errors(Yc, R)
  d <- er$e1^2 - er$e2^2
  yl <- er$ylag
  sd_i <- sqrt(rowSums((yl - rowMeans(yl))^2) / (ncol(yl) - 1))
  wv <- pnorm((yl - alpha_i) / sd_i)
  m0 <- rowMeans(d); m1 <- rowMeans(yl * d)
  Ew <- rowMeans(wv); Ey <- rowMeans(yl); Eyw <- rowMeans(yl * wv)
  N <- length(psi_i)
  if (calibrate == "unconditional") {
    cc <- psi_i - m0; ss <- rep(0, N)
  } else {
    det <- Eyw - Ew * Ey
    ss <- ((psi_i * alpha_i - m1) - Ey * (psi_i - m0)) / det
    cc <- (psi_i - m0) - ss * Ew
  }
  list(c = cc, s = ss, sd = sd_i)
}

#' Generate Forecast Simulation Data (Extended Design)
#'
#' @description Extended design of the revised Monte Carlo study (Section 5 of the paper):
#' \eqn{Y_{it} = \alpha_i(1-\rho_{k_i}) + \rho_{k_i}Y_{i,t-1} + \gamma_iF^y_t + U_{it}}, with
#' \eqn{\alpha_i\sim N(1, 0.25)}, \eqn{\gamma_i\sim U(0,1)} drawn once per call and \eqn{F^y_t} an AR(1)
#' common factor with coefficient 0.5. Both forecasters estimate an AR(1) by least squares on rolling
#' windows of length R: Forecaster 1 with an intercept, Forecaster 2 without. The forecast noise
#' \eqn{\sqrt{v_a(Y_{i,t-1})}\,\tilde\varepsilon^{(a)}_{it}}, with \eqn{\tilde\varepsilon^{(a)}} a
#' unit-variance AR(1) process with coefficient \eqn{\phi} and common factor with loading
#' \eqn{\lambda}, is calibrated unit by unit (\code{.calibrate_noise}) so that
#' \eqn{E(H_{i,t-1}\Delta L_{it}) = \psi_{k_i}(1,\alpha_i)'}. Under \eqn{\psi = 0}, EPA holds for every
#' unit, unconditionally and conditionally on \eqn{Y_{i,t-1}}.
#'
#' @param N Number of units.
#' @param Tobs Number of evaluation periods.
#' @param gamma Cluster membership (length N).
#' @param psi Expected loss differential for each cluster (length max(gamma)).
#' @param psi_break Optional; values of \code{psi} for the second half of the sample (breaks).
#' @param rho_vec AR coefficients by cluster.
#' @param alpha_mean,alpha_var Mean and variance of the individual effects.
#' @param fac_ar AR coefficient of the target factor.
#' @param R Rolling window length.
#' @param phi,lambda Parameters of the forecast noise.
#' @param Tcal Length of the calibration path.
#' @param calibrate "both" (default) or "unconditional".
#' @param burn Burn-in.
#'
#' @return A list with \code{DL} and \code{Ylag} (Tobs x N matrices, columns are units),
#'   \code{alpha}, \code{load} and \code{net_var} (Tobs x N matrix of net noise variances).
#' @export
generate_epa_simulation_data <- function(N, Tobs, gamma, psi, psi_break = NULL,
                                         rho_vec = c(0.1, 0.2, 0.3),
                                         alpha_mean = 1, alpha_var = 0.25, fac_ar = 0.5,
                                         R = 60, phi = 0.2, lambda = 0.2,
                                         Tcal = 20000, calibrate = c("both", "unconditional"),
                                         burn = 50) {
  calibrate <- match.arg(calibrate)
  if (length(gamma) != N) stop("gamma must have length N.")
  alpha_i <- rnorm(N, alpha_mean, sqrt(alpha_var))
  load_i <- runif(N)
  rho_i <- rho_vec[gamma]
  psi_i <- psi[gamma]

  cal1 <- .calibrate_noise(psi_i, alpha_i, rho_i, load_i, fac_ar, R, Tcal, calibrate, burn)
  cal2 <- if (is.null(psi_break)) cal1 else
    .calibrate_noise(psi_break[gamma], alpha_i, rho_i, load_i, fac_ar, R, Tcal, calibrate, burn)

  Y <- .sim_target(Tobs + R + 1, alpha_i, rho_i, load_i, fac_ar, burn)
  er <- .rolling_errors(Y, R)
  m <- ncol(er$e1)
  idx <- (m - Tobs + 1):m
  e1 <- er$e1[, idx, drop = FALSE]; e2 <- er$e2[, idx, drop = FALSE]
  yl <- er$ylag[, idx, drop = FALSE]

  regime2 <- matrix(rep(seq_len(Tobs) > Tobs / 2 & !is.null(psi_break), each = N), N)
  cc <- ifelse(regime2, cal2$c, cal1$c)
  ss <- ifelse(regime2, cal2$s, cal1$s)
  wv <- pnorm((yl - alpha_i) / cal1$sd)
  net <- cc + ss * wv
  # net variance can change sign with y only through s; split each part by its sign
  v1 <- pmax(cc, 0) + pmax(ss, 0) * wv
  v2 <- pmax(-cc, 0) + pmax(-ss, 0) * wv
  n1 <- sqrt(v1) * .sim_noise(N, Tobs, phi, lambda)
  n2 <- sqrt(v2) * .sim_noise(N, Tobs, phi, lambda)
  DL <- (e1 - n1)^2 - (e2 - n2)^2

  list(DL = t(DL), Ylag = t(yl), alpha = alpha_i, load = load_i, net_var = t(net))
}
