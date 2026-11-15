// Simulation helpers of generate_epa_simulation_data (see R/generate_panels.r): AR(1) recursion
// of the target and rolling-window least squares forecast errors. The random numbers are drawn in
// R, so the C++ and R versions use the same draws.
#include <Rcpp.h>
using namespace Rcpp;

// Y[, 1] = alpha; Y[, t] = alpha (1 - rho) + rho Y[, t - 1] + load Fy[t] + U[, t]
// [[Rcpp::export]]
NumericMatrix ar_target_cpp(NumericVector alpha, NumericVector rho, NumericVector load,
                            NumericVector Fy, NumericMatrix U) {
  int N = U.nrow(), n = U.ncol();
  NumericMatrix Y(N, n);
  for (int i = 0; i < N; ++i) Y(i, 0) = alpha[i];
  for (int t = 1; t < n; ++t)
    for (int i = 0; i < N; ++i)
      Y(i, t) = alpha[i] * (1 - rho[i]) + rho[i] * Y(i, t - 1) + load[i] * Fy[t] + U(i, t);
  return Y;
}

// Rolling sums of w consecutive values, from cumulative sums as in R's cumsum
static void roll(const std::vector<double>& x, int w, std::vector<double>& out) {
  int n = x.size();
  std::vector<double> cs(n + 1);
  long double s = 0;
  cs[0] = 0;
  for (int t = 0; t < n; ++t) { s += x[t]; cs[t + 1] = (double)s; }
  out.resize(n - w + 1);
  for (int j = 0; j <= n - w; ++j) out[j] = cs[j + w] - cs[j];
}

// Forecast errors of the AR(1) with intercept (e1) and without intercept (e2), rolling window R
// [[Rcpp::export]]
List rolling_errors_cpp(NumericMatrix Y, int R) {
  int N = Y.nrow(), n = Y.ncol();
  int w = R - 1;
  int m = (n - 1) - w;          // number of forecasts
  NumericMatrix e1(N, m), e2(N, m), yl(N, m), yy(N, m);
  std::vector<double> a(n - 1), b(n - 1), aa(n - 1), ab(n - 1), Sa, Sb, Saa, Sab;
  for (int i = 0; i < N; ++i) {
    for (int t = 0; t < n - 1; ++t) {
      a[t] = Y(i, t); b[t] = Y(i, t + 1); aa[t] = a[t] * a[t]; ab[t] = a[t] * b[t];
    }
    roll(a, w, Sa); roll(b, w, Sb); roll(aa, w, Saa); roll(ab, w, Sab);
    for (int j = 0; j < m; ++j) {
      double ylag = Y(i, R - 1 + j), y = Y(i, R + j);
      double beta1 = (w * Sab[j] - Sa[j] * Sb[j]) / (w * Saa[j] - Sa[j] * Sa[j]);
      double c1 = (Sb[j] - beta1 * Sa[j]) / w;
      double beta2 = Sab[j] / Saa[j];
      e1(i, j) = y - (c1 + beta1 * ylag);
      e2(i, j) = y - beta2 * ylag;
      yl(i, j) = ylag; yy(i, j) = y;
    }
  }
  return List::create(Named("e1") = e1, Named("e2") = e2, Named("ylag") = yl, Named("y") = yy);
}
