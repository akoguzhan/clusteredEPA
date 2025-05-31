#' @title Overall Equal Predictive Ability (O-EPA) Test
#'
#' @description Computes the test statistic and p-value for the Overall Equal Predictive Ability (O-EPA) hypothesis 
#' using cross-sectional averages of transformed loss differentials. The test can be unconditional or conditional 
#' on testing functions. Long-run variance estimation can be done using either Newey-West or EWC estimators.
#'
#' @param df A data frame containing the panel dataset with unit and time identifiers, loss differentials, and (optionally) testing functions.
#' @param DL A character string naming the column in `df` that contains the loss differential \eqn{\Delta L_{it}}.
#' @param H A character vector of column names corresponding to testing functions \eqn{H_{it}}, or \code{NULL} for the unconditional case.
#' @param id A character string naming the column in `df` that contains the unit (cross-sectional) identifier.
#' @param time A character string naming the column in `df` that contains the time index.
#' @param lrv A character string indicating which long-run variance estimator to use. Either \code{"NeweyWest"} or \code{"EWC"}.
#' @param lrv_par A numeric value specifying the tuning parameter: lag truncation for Newey-West or degrees of freedom for EWC.
#'
#' @return A list with the following elements:
#' \describe{
#'   \item{\code{S.oepa}}{Estimated long-run variance matrix of the test moment vector.}
#'   \item{\code{W.oepa}}{Test statistic value.}
#'   \item{\code{p.oepa}}{P-value of the O-EPA test.}
#' }
#'
#' @examples
#' # Suppose 'df' contains columns 'DL', 'H1', 'id', and 'time':
#' overall_EPA_test(df, DL = "DL", H = c("H1"), id = "id", time = "time", lrv = "NeweyWest", lrv_par = 4)
#'
#' @export
#' 
overall_EPA_test = function(df,DL,H,id,time,lrv,lrv_par) {
  # Preliminaries
  Tobs = length(unique(df[,time]))
  Z_names = c(DL,H)
  K = length(Z_names)
  .data = df[,c(Z_names,id,time)]
  if (!is.null(H)) {
    for (k in 1:length(K-1)) {
      .data[,H[k]] = .data[,DL]*.data[,H[k]]
    }
  }
  
  # Taking cross-sectional averages
  if (K==1) {
    Z = as.matrix(c(.data[,Z_names]))
  } else {
    Z = as.matrix(.data[,Z_names])
  }
  Zbar_Nt = aggregate_matrix(Z,.data[,time])
  
  # Test statistic
  mu.hat.oepa = colMeans(Zbar_Nt)
  if (identical(lrv,"NeweyWest")) {
    S.oepa = NeweyWest(Zbar_Nt,lrv_par)
    W.oepa = Tobs*t(mu.hat.oepa)%*%solve(S.oepa)%*%(mu.hat.oepa)
    p.oepa = pchisq(q=W.oepa,df=K,lower.tail=FALSE)
  } else if (identical(lrv,"EWC")) {
    adj.term_overall = ((lrv_par-K+1)/lrv_par)/K
    S.oepa = EWC(Zbar_Nt,lrv_par)
    W.oepa = adj.term_overall*Tobs*t(mu.hat.oepa)%*%solve(S.oepa)%*%(mu.hat.oepa)
    p.oepa = pf(W.oepa,K,lrv_par-K+1,lower.tail=FALSE)
  }
  
  # Return results
  result_list = list("S.oepa" = S.oepa,
                     "W.oepa" = W.oepa,
                     "p.oepa" = p.oepa)
  
  return(result_list)
}
