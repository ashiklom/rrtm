#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
arma::mat solvearray(const arma::cube& a, const arma::mat& b) {
  unsigned int nwl = a.n_slices, s = a.n_cols;
  unsigned int w;
  arma::mat result(s, nwl);
  // arma::sp_mat M;

  for (w = 0; w < nwl; w++) {
    result.col(w) = arma::solve(a.slice(w), b.col(w), arma::solve_opts::fast);
  }

  return result;
}
