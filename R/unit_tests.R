# =============================================================================
# Unit-level tests of joint EPA (EPA for every unit) and power enhancement.
# Used only as benchmarks and robustness checks: the paper's tests are the selective C-EPA tests.
# =============================================================================

#' Unit-Level Wald Statistics with the OS Long-Run Variance
#'
#' @description For each unit i, computes \eqn{W_i = T\,\bar a_i'\hat\Omega_i^{-1}\bar a_i}, where
#' \eqn{\bar a_i} is an estimate (by default the time average of unit i's series) and
#' \eqn{\hat\Omega_i} is the orthonormal series (OS) estimator with B cosine terms of the long-run
#' variance of the T x P series \eqn{A_i}. Under fixed-B asymptotics,
#' \eqn{F_i = (B - P + 1) W_i / (P B)} is \eqn{F(P, B - P + 1)}.
#'
#' @param A List of N matrices (T x P), the series whose long-run variance is estimated.
#' @param est List of N vectors of length P, the estimates; default: the column means of A.
#' @param B Number of cosine terms.
#'
#' @return A list with \code{W}, \code{F}, unit p-values \code{p}, and the degrees of freedom
#'   \code{d1}, \code{d2} of the reference F law.
#' @keywords internal
.unit_wald_os <- function(A, est = NULL, B) {
  N <- length(A)
  P <- ncol(A[[1]])
  Tobs <- nrow(A[[1]])
  if (is.null(est)) est <- lapply(A, colMeans)
  Cb <- sqrt(2 / Tobs) * cos(outer(seq_len(Tobs) - 0.5, seq_len(B)) * pi / Tobs)  # T x B
  W <- vapply(seq_len(N), function(i) {
    Lam <- crossprod(Cb, sweep(A[[i]], 2, colMeans(A[[i]])))                      # B x P
    Tobs * drop(t(est[[i]]) %*% solve(crossprod(Lam) / B, est[[i]]))
  }, 0)
  Fi <- (B - P + 1) / (P * B) * W
  list(W = W, F = Fi, p = stats::pf(Fi, P, B - P + 1, lower.tail = FALSE), d1 = P, d2 = B - P + 1)
}

#' Mean and Variance of the F Distribution
#'
#' @description Closed-form mean \eqn{d_2/(d_2-2)} and variance of the \eqn{F(d_1, d_2)} distribution.
#'
#' @param d1,d2 Degrees of freedom (\code{d2 > 4}).
#' @return A list with the mean \code{m} and the variance \code{v} of \eqn{F(d_1, d_2)}.
#' @keywords internal
.F_moments <- function(d1, d2) {
  list(m = d2 / (d2 - 2), v = 2 * d2^2 * (d1 + d2 - 2) / (d1 * (d2 - 2)^2 * (d2 - 4)))
}

#' Sum of Squared Cross-Sectional Correlations
#'
#' @description Computes \eqn{\sum_{i\neq j}\tilde\rho_{ij}^2}, where \eqn{\tilde\rho_{ij}^2} is the
#' squared contemporaneous correlation between the standardized series of units i and j, bias
#' corrected by \eqn{-1/T}. For P > 1 the P x P cross-correlation blocks \eqn{C_{ij}} enter through
#' \eqn{\mathrm{tr}(C_{ij}C_{ij}')/P}.
#'
#' @param Zl List of N matrices (T x P).
#' @return A scalar, bounded below by \eqn{-N/2}.
#' @keywords internal
.cross_rho2 <- function(Zl) {
  N <- length(Zl)
  P <- ncol(Zl[[1]])
  Tobs <- nrow(Zl[[1]])
  U <- do.call(cbind, lapply(Zl, function(Zi) {
    Zc <- sweep(Zi, 2, colMeans(Zi))
    e <- eigen(crossprod(Zc) / Tobs, symmetric = TRUE)
    Zc %*% e$vectors %*% diag(1 / sqrt(e$values), P) %*% t(e$vectors)
  }))
  R2 <- (crossprod(U) / Tobs)^2 - 1 / Tobs
  blk <- rep(seq_len(N), each = P)
  R2[outer(blk, blk, "==")] <- 0
  max(sum(R2) / P, -N / 2)
}

#' Unit-Level Tests of Joint EPA
#'
#' @description Tests of the hypothesis that EPA holds for every unit, \eqn{\mu_i = 0} for all i,
#' based on unit-level Wald statistics with the OS long-run variance (fixed B):
#' \describe{
#'   \item{\code{p_maxt}}{Bonferroni test over units, \eqn{\min(1, N\min_i p_i)}; valid under any
#'     cross-sectional dependence.}
#'   \item{\code{p_delta}}{Standardized average of the unit statistics, as in the \eqn{\Delta} test of
#'     Pesaran and Yamagata (2008), one-sided, with the exact fixed-B moments of the F law; assumes
#'     cross-sectional independence.}
#'   \item{\code{p_delta_cc}}{The same with the variance of the sum corrected for cross-sectional
#'     correlation, \eqn{\mathrm{var}(F)[N + \sum_{i\neq j}\tilde\rho_{ij}^2]}, as in the alpha test of
#'     Pesaran and Yamagata (2024), with bias-corrected rather than thresholded correlations.}
#'   \item{\code{p_choi}}{Inverse normal combination of the unit p-values (Choi, 2001); assumes
#'     cross-sectional independence.}
#'   \item{\code{p_cce_delta}, \code{p_cce_max}}{Defactored tests: each unit's series is regressed
#'     on an intercept and the cross-sectional averages; since \eqn{\mu_i = \alpha_i + \beta_i'\bar\mu},
#'     joint EPA holds if and only if all \eqn{\alpha_i = 0} and \eqn{\bar\mu = 0}. The intercept
#'     statistics (OS variance of the influence series of the intercept) are combined by \eqn{\Delta}
#'     or by Bonferroni, and the result is combined with the O-EPA p-value by Bonferroni.}
#' }
#'
#' @param Z NT x P matrix of loss differentials or moments.
#' @param id,time NT vectors of unit and time identifiers (balanced panel).
#' @param B Number of cosine terms of the OS estimator (default \code{epa_default_B}).
#' @param p_oepa O-EPA p-value with the same B (default: computed with \code{overall_EPA_test}).
#'
#' @return A list with the p-values above, the unit p-values \code{unit_p}, the unit statistics
#'   \code{unit_W} and \code{B}.
#' @references
#' Choi, I. (2001). Unit root tests for panel data. Journal of International Money and Finance.
#'
#' Pesaran, M. H. and Yamagata, T. (2008). Testing slope homogeneity in large panels. Journal of
#' Econometrics.
#'
#' Pesaran, M. H. and Yamagata, T. (2024). Testing for alpha in linear factor pricing models with a
#' large number of securities. Journal of Financial Econometrics.
#' @export
epa_unit_tests <- function(Z, id, time, B = NULL, p_oepa = NULL) {
  Z <- as.matrix(Z)
  o <- order(id, time)
  Z <- Z[o, , drop = FALSE]; id <- id[o]; time <- time[o]
  ids <- unique(id)
  N <- length(ids)
  Tobs <- length(unique(time))
  if (length(id) != N * Tobs) stop("epa_unit_tests requires a balanced panel.")
  if (is.null(B)) B <- epa_default_B(Z, id, time)
  if (is.null(p_oepa)) p_oepa <- overall_EPA_test(Z, id, time, lrv_par = B)$pval
  Zl <- lapply(seq_len(N), function(i) Z[((i - 1) * Tobs + 1):(i * Tobs), , drop = FALSE])

  us <- .unit_wald_os(Zl, B = B)
  mv <- .F_moments(us$d1, us$d2)
  Delta <- sum(us$F - mv$m) / sqrt(N * mv$v)
  Delta_cc <- sum(us$F - mv$m) / sqrt(mv$v * (N + .cross_rho2(Zl)))
  Zc <- sum(stats::qnorm(pmin(pmax(us$p, 1e-300), 1 - 1e-16))) / sqrt(N)

  # defactored intercepts: regression on (1, cross-sectional averages)
  Zbar <- Reduce(`+`, Zl) / N
  X <- cbind(1, Zbar)
  xb <- colMeans(Zbar)
  Xc <- sweep(Zbar, 2, xb)
  wt <- drop(1 - Xc %*% solve(crossprod(Xc) / Tobs, xb))   # influence weights of the intercept
  XtXi <- solve(crossprod(X))
  cce <- lapply(Zl, function(Zi) {
    bh <- XtXi %*% crossprod(X, Zi)
    list(a = bh[1, ], inf = wt * (Zi - X %*% bh))
  })
  ua <- .unit_wald_os(lapply(cce, `[[`, "inf"), lapply(cce, `[[`, "a"), B)
  p_alpha_delta <- stats::pnorm(sum(ua$F - mv$m) / sqrt(N * mv$v), lower.tail = FALSE)
  p_alpha_max <- min(1, N * min(ua$p))

  list(p_maxt = min(1, N * min(us$p)),
       p_delta = stats::pnorm(Delta, lower.tail = FALSE),
       p_delta_cc = stats::pnorm(Delta_cc, lower.tail = FALSE),
       p_choi = stats::pnorm(Zc),
       p_cce_delta = min(1, 2 * min(p_alpha_delta, p_oepa)),
       p_cce_max = min(1, 2 * min(p_alpha_max, p_oepa)),
       unit_p = us$p, unit_W = us$W, B = B)
}

#' Power Enhancement of a Selective p-Value by Unit-Level Screening
#'
#' @description Weighted Bonferroni combination of a selective p-value \eqn{p_0} and the unit-level
#' max test \eqn{p_{max}} (Fan, Liao and Yao, 2015, with a finite-sample screening budget):
#' at level \eqn{\alpha}, the test rejects if \eqn{p_0 \le \alpha - \gamma} or some unit has
#' \eqn{p_i \le \gamma / N}. The size is at most \eqn{\alpha} whenever both components are valid.
#' Valid for joint EPA only: under C-EPA, units within a cluster may have non-zero means.
#'
#' @param p0 Selective p-value.
#' @param p_maxt Bonferroni p-value over units (\code{epa_unit_tests()$p_maxt}).
#' @param alpha Level at which the combination is calibrated.
#' @param gam Part of \code{alpha} spent on the screening.
#' @return The combined p-value.
#' @references Fan, J., Liao, Y. and Yao, J. (2015). Power enhancement in high-dimensional
#'   cross-sectional tests. Econometrica.
#' @export
power_enhance <- function(p0, p_maxt, alpha = 0.05, gam = 0.01) {
  w <- 1 - gam / alpha
  min(1, p0 / w, p_maxt / (1 - w))
}
