#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]

//' Solve the matrix equation AX = B for X for an array of matrices A
//'
//' @param a 3D array. For dimensions x,y,z, iteration is over dimension z ("slices")
//' @param b Matrix with dimensions y,x
//' @return Matrix of solutions (x,z)
//' @import Rcpp
// [[Rcpp::export]]
arma::mat solvearray(const arma::cube& a, const arma::mat& b) {
  unsigned int nwl = a.n_slices, s = a.n_cols;
  unsigned int w;
  arma::mat result(s, nwl);

  for (w = 0; w < nwl; w++) {
    result.col(w) = arma::solve(a.slice(w), b.col(w), arma::solve_opts::fast);
  }

  return result;
}
