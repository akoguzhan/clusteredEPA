#' Default Number of Cosines for all EPA Tests on a Sample
#'
#' @description \eqn{B = \lceil c_BT^{2/3}\rceil} with \eqn{c_B} from \code{PanelKmeansInference::choose_B}
#' applied to the cross-sectional average \eqn{\bar Z_t}, bounded below by \eqn{K_{max}P + 4}. Using the
#' same B for all tests on a sample follows Section 4.3 of the paper.
#'
#' @param Z NT x P matrix.
#' @param id,time NT vectors.
#' @param Kmax Largest number of clusters used by any test.
#' @param cB Optional constant \eqn{c_B}.
#' @return Integer B.
#' @export
epa_default_B <- function(Z, id, time, Kmax = 5, cB = NULL) {
  Z <- as.matrix(Z)
  Zbar_t <- aggregate_matrix(Z, time)
  choose_B(Zbar_t, cB = cB, Bmin = Kmax * ncol(Z) + 4)
}

#' Overall Equal Predictive Ability (O-EPA) Test for Panel Data
#'
#' @description Wald test of a zero mean of the cross-sectional average \eqn{\bar Z_t}, with the OS
#' (EWC) long-run variance and F critical values, or Newey-West and chi-squared critical values.
#'
#' @param Z NT x P matrix.
#' @param id,time NT vectors.
#' @param lrv "EWC" (default) or "NeweyWest".
#' @param lrv_par Number of cosines B (EWC) or lag (NW). Default for EWC: \code{choose_B}.
#' @param cB Optional constant \eqn{c_B}.
#' @param demean "global" or "local" demeaning for the long-run variance (EWC).
#'
#' @return A list with S, W, pval, lrv_par.
#' @export
overall_EPA_test <- function(Z, id, time, lrv = "EWC", lrv_par = NULL, cB = NULL, demean = "global") {
  if (!is.matrix(Z)) stop("Z must be a matrix.")
  NT <- nrow(Z)
  if (length(id) != NT || length(time) != NT) stop("Length of 'id' and 'time' must match number of rows in Z.")
  Tobs <- length(unique(time))
  P <- ncol(Z)
  if (Tobs < P) stop("Number of time periods (Tobs) must be at least as large as the number of moments (P).")

  Zbar_Nt <- aggregate_matrix(Z, time)
  mu_hat_oepa <- colMeans(Zbar_Nt)

  if (identical(lrv, "NeweyWest")) {
    S <- NeweyWest(Zbar_Nt, lrv_par)
    S_inv <- chol2inv(chol(S))
    W <- as.numeric(Tobs * t(mu_hat_oepa) %*% S_inv %*% mu_hat_oepa)
    pval <- pchisq(W, df = P, lower.tail = FALSE)
  } else if (identical(lrv, "EWC")) {
    if (is.null(lrv_par)) lrv_par <- choose_B(Zbar_Nt, cB = cB, Bmin = P + 4)
    S <- EWC(Zbar_Nt, lrv_par, demean = demean)$S
    adj_term <- ((lrv_par - P + 1) / lrv_par) / P
    S_inv <- chol2inv(chol(S))
    W <- as.numeric(adj_term * Tobs * t(mu_hat_oepa) %*% S_inv %*% mu_hat_oepa)
    pval <- pf(W, df1 = P, df2 = lrv_par - P + 1, lower.tail = FALSE)
  } else {
    stop("Invalid 'lrv' value. Use 'NeweyWest' or 'EWC'.")
  }
  return(list(S = S, W = W, pval = pval, lrv_par = lrv_par))
}

#' Clustered EPA (C-EPA) Test with Known Group Memberships
#'
#' @description Wald test of a zero mean of the stacked cluster averages
#' \eqn{(\bar Z_{1t}',\dots,\bar Z_{Kt}')'} for given clusters, with the OS (EWC) long-run variance
#' and F critical values, or Newey-West and chi-squared critical values.
#'
#' @param Z Matrix (NT x P) of loss differentials or moment conditions.
#' @param id Vector of length NT for unit identifiers.
#' @param time Vector of length NT for time identifiers.
#' @param gamma Vector of length N (units in increasing id order), cluster assignments in 1..K.
#' @param lrv "EWC" (default) or "NeweyWest".
#' @param lrv_par Number of cosines (EWC) or lag (NW). Default for EWC: \code{choose_B} on
#'   \eqn{\bar Z_t} with lower bound \eqn{KP+4}.
#' @param cB Optional constant \eqn{c_B}.
#' @param demean "global" or "local" demeaning for the long-run variance (EWC).
#'
#' @return A list with W, S, pval, lrv_par.
#' @export
epa_clustered_known <- function(Z, id, time, gamma, lrv = "EWC", lrv_par = NULL, cB = NULL,
                                demean = "global") {
  if (!is.matrix(Z)) stop("Z must be a matrix.")
  P <- ncol(Z)
  N <- length(unique(id))
  Tobs <- length(unique(time))
  if (length(gamma) != N) stop("Length of gamma must match number of unique units.")

  df <- data.frame(Z, id = id, time = time)
  df <- df[order(df$time, df$id), ]
  Z <- as.matrix(df[, 1:P])
  id <- df$id
  time <- df$time

  id_unique <- sort(unique(id))
  gamma_long <- gamma[match(id, id_unique)]
  K <- length(unique(gamma))
  KP <- K * P
  if (Tobs < KP) stop("Number of time periods (Tobs) must be at least as large as K*P.")
  Z_by_group <- matrix(NA, nrow = Tobs, ncol = KP)
  for (k in 1:K) {
    idx_k <- gamma_long == k
    Z_by_group[, ((k - 1) * P + 1):(k * P)] <- aggregate_matrix(Z[idx_k, , drop = FALSE], time[idx_k])
  }
  mu_hat <- colMeans(Z_by_group)

  if (lrv == "NeweyWest") {
    S <- NeweyWest(Z_by_group, lrv_par = lrv_par)
    S_inv <- tryCatch(chol2inv(chol(S)), error = function(e) stop("S is not positive definite."))
    W <- as.numeric(Tobs * t(mu_hat) %*% S_inv %*% mu_hat)
    pval <- pchisq(W, df = KP, lower.tail = FALSE)
  } else if (lrv == "EWC") {
    if (is.null(lrv_par)) lrv_par <- choose_B(aggregate_matrix(Z, time), cB = cB, Bmin = KP + 4)
    S <- EWC(Z_by_group, lrv_par = lrv_par, demean = demean)$S
    adj_term <- ((lrv_par - KP + 1) / lrv_par) / KP
    S_inv <- tryCatch(chol2inv(chol(S)), error = function(e) stop("S is not positive definite."))
    W <- as.numeric(adj_term * Tobs * t(mu_hat) %*% S_inv %*% mu_hat)
    pval <- pf(W, df1 = KP, df2 = lrv_par - KP + 1, lower.tail = FALSE)
  } else {
    stop("Unknown lrv. Use 'EWC' or 'NeweyWest'.")
  }
  return(list(W = W, S = S, pval = pval, lrv_par = lrv_par))
}

#' Clustered EPA (C-EPA) Test with Sample Splitting
#'
#' @description Clusters are estimated on the first \code{prop} share of the periods (training
#' sample, K chosen by the IC if \code{Kmax} is given), the first \code{maxlag} periods of the
#' remaining sample are dropped (gap), and the C-EPA test with known clusters is computed on the
#' testing sample. B is chosen on the testing sample unless \code{lrv_par} is given. Columns of
#' Z that are identically zero in the training sample (e.g. the product with a dummy that is zero
#' there) are not used to estimate the clusters.
#'
#' @inheritParams epa_clustered_known
#' @param K,Kmax Number of clusters or largest number for the IC.
#' @param prop Training proportion \eqn{\gamma}.
#' @param maxlag Gap between the samples (default \eqn{\lfloor\sqrt{prop\,T}\rfloor}).
#' @param Ninit,iter.max,kappa,init_type Panel Kmeans settings.
#'
#' @return A list with clustering, W, S, pval, lrv_par.
#' @export
epa_clustered_split <- function(
  Z, id, time,
  K = NULL,
  Kmax = NULL,
  prop = 0.5,
  lrv = "EWC",
  lrv_par = NULL,
  Ninit = 10,
  iter.max = 100,
  kappa = 1.5,
  maxlag = NULL,
  cB = NULL,
  demean = "global",
  init_type = "forgy") {

  split <- split_panel_matrix(Z, id = id, time = time, prop = prop, maxlag = maxlag)
  # Columns that are identically zero in the training sample (e.g. d x dummy when the dummy is 0
  # there) carry no information on the clusters and make the IC determinant zero: cluster on the
  # other columns.
  ZTr <- as.matrix(split$Z_Tr)
  nz <- colSums(abs(ZTr)) > 0
  if (!any(nz)) stop("The training sample is identically zero.")
  km_res <- panel_kmeans_estimation(Z = ZTr[, nz, drop = FALSE], id = split$id_Tr, time = split$time_Tr,
                                    K = K, Kmax = Kmax, Ninit = Ninit, iter.max = iter.max,
                                    kappa = kappa, init_type = init_type)
  gamma_Tr <- km_res$final_cluster
  gamma_Te <- gamma_Tr[as.character(sort(unique(split$id_Te)))]
  if (any(is.na(gamma_Te))) stop("Some units in the testing sample were not found in the training clustering.")

  test_res <- epa_clustered_known(Z = split$Z_Te, id = split$id_Te, time = split$time_Te,
                                  gamma = as.integer(gamma_Te), lrv = lrv, lrv_par = lrv_par,
                                  cB = cB, demean = demean)
  return(append(list(clustering = gamma_Tr), test_res))
}

#' Selective Inference for Clustered EPA
#'
#' @description Selective C-EPA tests with unknown clusters (Section 4 of the paper):
#' \enumerate{
#'   \item Panel Kmeans with \code{Ninit} random initializations; K fixed or chosen by the IC
#'     from 2 to \code{Kmax}.
#'   \item Selective p-values of all pairwise equality hypotheses, conditional on the selected
#'     partition and \eqn{\hat K_{IC}} (given the data-independent initializations), with the
#'     decomposition in the metric of the long-run covariance matrix.
#'   \item Homogeneity: merged pairwise p-values (\eqn{F_{homo,r}}) and the joint selective test.
#'   \item O-EPA test with the same B.
#'   \item C-EPA: \eqn{F_{SI,r}} merging the pairwise and O-EPA p-values (main test, Theorem 2),
#'     the joint selective C-EPA test (\eqn{\mathsf R = I_{KP}}, Remark 7) and the weighted
#'     Bonferroni combination of the homogeneity and O-EPA p-values (Remark 9).
#'   \item Optionally (\code{cluster_specific = TRUE}), the selective tests that each selected
#'     cluster mean is zero, with Holm-adjusted p-values (Section 5.3), which identify the clusters
#'     for which EPA fails after a rejection.
#' }
#'
#' @param Z NT x P matrix of loss differentials or moments.
#' @param id,time NT vectors of unit and time identifiers.
#' @param K Number of clusters (alternative to \code{Kmax}).
#' @param Kmax Largest number of clusters for the IC.
#' @param pairs Optional matrix of cluster pairs.
#' @param pcombine_fun Merging function (default "Genmean_neq").
#' @param r Order of the generalized mean (default -20).
#' @param lrv Only "EWC" is supported.
#' @param lrv_par Number of cosines B (default: \code{choose_B} with lower bound \eqn{K_{max}P+4}).
#' @param cB Optional constant \eqn{c_B}.
#' @param Ninit Number of random initializations.
#' @param iter.max Maximum number of Lloyd iterations.
#' @param kappa IC tuning constant.
#' @param metric "xi" (default), "shrink" or "euclidean" (see \code{selective_test}).
#' @param conditioning "partition" (default: condition on the selected partition and K_IC only)
#'   or "paths" (see \code{selective_test}).
#' @param condition_K,condition_starts Options for \code{conditioning = "paths"}.
#' @param ref Reference law of the truncated statistic: "chi" (default, \eqn{\chi_q}) or "F"
#'   (Hotelling law with B degrees of freedom; see \code{selective_test}).
#' @param init_type "forgy" (default) or "partition" (see \code{pk_cluster}).
#' @param w_bonf Weight of the homogeneity p-value in the weighted Bonferroni combination.
#' @param demean "global" or "local" demeaning for the long-run variance.
#' @param fit Optional precomputed \code{pk_fit} object for the same Z.
#' @param cluster_specific Also compute the cluster-specific selective tests (default FALSE; K
#'   additional truncation sets).
#'
#' @return A list with pairwise_pvalues, pairs, overall_pval, hom_test_pval (merged),
#'   joint_hom_pval, joint_cepa_pval, wbonf_pval, pval (\eqn{F_{SI,r}}), clustering, K, B, in_set
#'   (diagnostic), fit; with \code{cluster_specific = TRUE} also cluster_pvals, cluster_holm
#'   (Holm-adjusted), cluster_holm_pval (smallest adjusted p-value), cluster_centers (K x P) and
#'   cluster_sizes (see \code{panel_bycluster_tests}).
#' @export
epa_clustered_selective <- function(
  Z, id, time,
  K = NULL,
  Kmax = NULL,
  pairs = NULL,
  pcombine_fun = "Genmean_neq",
  r = -20,
  lrv = "EWC",
  lrv_par = NULL,
  cB = NULL,
  Ninit = 10,
  iter.max = 100,
  kappa = 1.5,
  metric = "xi",
  conditioning = "partition",
  condition_K = TRUE,
  condition_starts = TRUE,
  ref = "chi",
  w_bonf = 0.5,
  demean = "global",
  init_type = "forgy",
  fit = NULL,
  cluster_specific = FALSE
) {
  if (lrv != "EWC") stop("Selective tests are implemented with the OS (EWC) estimator only.")
  if (is.null(fit)) {
    fit <- panel_selective_fit(Z, id, time, K = K, Kmax = Kmax, Ninit = Ninit, iter.max = iter.max,
                               kappa = kappa, lrv_par = lrv_par, cB = cB, demean = demean,
                               init_type = init_type)
  }
  st <- panel_selective_tests(fit, pairs = pairs, pcombine_fun = pcombine_fun, r = r, joint = TRUE,
                              metric = metric, conditioning = conditioning, condition_K = condition_K,
                              condition_starts = condition_starts, ref = ref)
  oepa <- overall_EPA_test(Z = Z, id = id, time = time, lrv = "EWC", lrv_par = fit$B, demean = demean)
  overall_pval <- oepa$pval

  all_pvals <- c(st$pairwise_pvalues, overall = overall_pval)
  selective_pval <- merge_pvalues(all_pvals, pcombine_fun, r)

  out <- list(
    pairwise_pvalues = st$pairwise_pvalues,
    pairs = st$pairs,
    overall_pval = overall_pval,
    hom_test_pval = st$hom_pval,
    joint_hom_pval = st$joint_hom_pval,
    joint_cepa_pval = st$joint_cepa_pval,
    wbonf_pval = weighted_bonferroni(st$hom_pval, overall_pval, w_bonf),
    pval = selective_pval,
    clustering = fit$clustering,
    K = fit$K,
    B = fit$B,
    in_set = st$in_set,
    fit = fit
  )
  if (cluster_specific) {
    cs <- panel_bycluster_tests(fit = fit, metric = metric, conditioning = conditioning,
                                condition_K = condition_K, condition_starts = condition_starts,
                                ref = ref)
    out$cluster_pvals <- cs$bycluster_pvals
    out$cluster_holm <- cs$bycluster_holm
    out$cluster_holm_pval <- cs$holm_pval
    out$cluster_centers <- cs$centers
    out$cluster_sizes <- cs$sizes
    out$in_set <- c(out$in_set, cs$in_set)
  }
  out
}

#' Selective Tests that Each Estimated Cluster Has Zero Mean
#'
#' @description Merges the selective p-values of \eqn{\theta^0_k(\hat{\mathcal C}) = 0}, k = 1..K.
#' @inheritParams epa_clustered_selective
#' @return A list with pval, bycluster_pvals, clustering.
#' @export
epa_clustered_selective_mode2 <- function(
  Z, id, time,
  K = NULL,
  Kmax = NULL,
  pcombine_fun = "Genmean_neq",
  r = -20,
  lrv = "EWC",
  lrv_par = NULL,
  cB = NULL,
  Ninit = 10,
  iter.max = 100,
  kappa = 1.5,
  metric = "xi",
  conditioning = "partition",
  init_type = "forgy",
  fit = NULL
) {
  res <- panel_bycluster_tests(Z = Z, id = id, time = time, K = K, Kmax = Kmax, Ninit = Ninit,
                               iter.max = iter.max, lrv = lrv, lrv_par = lrv_par, kappa = kappa,
                               cB = cB, metric = metric, conditioning = conditioning,
                               init_type = init_type, fit = fit)
  list(pval = merge_pvalues(res$bycluster_pvals, pcombine_fun, r),
       bycluster_pvals = res$bycluster_pvals, clustering = res$clustering)
}

#' Clark-West Test
#'
#' @description Clark-West test of equal predictive accuracy for nested models: the adjusted loss
#' differential \eqn{f_t = (y_t - y_{1t})^2 - [(y_t - y_{2t})^2 - (y_{1t} - y_{2t})^2]} is tested with
#' a t statistic using the Newey-West variance, and the p-value is \eqn{1 - \Phi(T)}.
#'
#' @param y1,y2 Forecasts.
#' @param y Realizations.
#' @param lrv_par NW lag.
#' @return A list with T and pval.
#' @export
CWtest <- function(y1, y2, y, lrv_par = 0) {
  if (!(length(y1) == length(y2) && length(y2) == length(y))) {
    stop("Input vectors y1, y2, and y must have the same length.")
  }
  P <- length(y1)
  e1 <- (y - y1)^2
  e2 <- (y - y2)^2
  f <- e1 - (e2 - (y1 - y2)^2)
  mu <- sum(f) / P
  sigma2 <- NeweyWest(f, lrv_par = lrv_par)
  T_stat <- mu / sqrt(sigma2 / P)
  pval <- 1 - pnorm(T_stat)
  return(list(T = T_stat, pval = pval))
}

#' Giacomini-White Test
#'
#' @description Regresses the loss differential on the instruments by least squares and computes
#' the Wald statistic \eqn{\hat\beta'\hat S^{-1}\hat\beta}, with \eqn{\hat S} the Newey-West variance
#' of \eqn{x_te_t}, and chi-squared p-value with \code{ncol(x_t)} degrees of freedom.
#'
#' @param d_t Loss differential.
#' @param x_t Instruments (default: a constant).
#' @param max_lag NW lag.
#' @return A list with statistic, pval, S_hat.
#' @export
GWtest <- function(d_t, x_t = NULL, max_lag = 0) {
  if (!is.vector(d_t)) stop("d_t must be a vector")
  if (anyNA(d_t)) stop("d_t contains missing values")
  P <- length(d_t)
  if (is.null(x_t)) {
    x_t <- matrix(1, nrow = P, ncol = 1)
  } else {
    x_t <- as.matrix(x_t)
    if (nrow(x_t) != P) stop("x_t must have same number of rows as d_t")
    if (anyNA(x_t)) stop("x_t contains missing values")
  }
  X <- x_t
  Y <- d_t
  beta_hat <- qr.solve(t(X) %*% X, t(X) %*% Y)
  e_t <- Y - X %*% beta_hat
  S_hat <- NeweyWest(X * e_t, lrv_par = max_lag)
  stat <- t(beta_hat) %*% solve(S_hat) %*% beta_hat
  df <- ncol(X)
  pval <- 1 - pchisq(stat, df)
  return(list(statistic = as.numeric(stat), pval = as.numeric(pval), S_hat = S_hat))
}
