% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/test_functions.R
\name{GWtest}
\alias{GWtest}
\title{Giacomini and White (2006) Predictive Ability Test}
\usage{
GWtest(d_t, x_t = NULL, max_lag = 0)
}
\arguments{
\item{d_t}{A numeric vector of out-of-sample loss differentials (L1 - L2) of length P.}

\item{x_t}{Optional matrix or vector of instruments (predictors) to condition on. Defaults to 1 (unconditional test).}

\item{max_lag}{Integer; lag length for Newey-West LRV estimation.}
}
\value{
A list with the test statistic, p-value, and long-run variance estimate.
}
\description{
Implements the conditional predictive ability (CPA) test for two competing forecasts
using the time series of out-of-sample loss differentials and a Newey-West long-run variance estimate.
}
