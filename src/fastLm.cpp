// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

using namespace Rcpp;

// [[Rcpp::export]]
List fastLm(Rcpp::NumericVector y, Rcpp::NumericMatrix x)
{
  int n = x.nrow(), k = x.ncol();
  
  arma::mat xMatrix{ x.begin(), n, k, false };
  arma::colvec yVector{ y.begin(), y.size(), false };
  
  arma::colvec coef = arma::solve(xMatrix, yVector);
  
  return List::create(coef);
}