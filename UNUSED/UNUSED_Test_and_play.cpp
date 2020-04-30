// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

using namespace Rcpp;

arma::mat xMatrix(arma::colvec xVector, int polyOrder);

// test function to play around with
// [[Rcpp::export]]
List test(NumericVector x, double q)
{
  arma::colvec xVector{ x.begin(), x.size(), false };
  
  arma::colvec out = pow(3, 2) * xVector;
  
  return List::create(out + 2);
}

int test0(int x, int (*testFcn)(int))
{
  return testFcn(x) - 1;
}

// [[Rcpp::export]]
arma::vec test2(double q)
{
  arma::vec qVec{ arma::vec(5).fill(q) };
  return pow(qVec, 2) - 1 ;
}

// [[Rcpp::export]]
arma::mat test3(int colIndex, int bndw)
{
  arma::colvec  uBound(arma::regspace(-colIndex, bndw) / bndw);
  arma::mat     xMatBound{ xMatrix(uBound, 3) };
  return xMatBound;
}


