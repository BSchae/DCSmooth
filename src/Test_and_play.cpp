// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

using namespace Rcpp;

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

int test2(int x)
{
  return x + 1;
}

// [[Rcpp::export]]
int test3(int x)
{
  return test0(x, test2);
}


