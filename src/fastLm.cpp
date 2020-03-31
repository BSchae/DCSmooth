// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

using namespace Rcpp;

// y is the vector of endogenous variables
// x is the matrix of exogenous variables

// [[Rcpp::export]]
List fastLm(NumericVector y, NumericMatrix x, NumericMatrix w)
{
  int n = x.nrow(), k = x.ncol();
  
  arma::mat xMatrix{ x.begin(), n, k, false };
  arma::colvec yVector{ y.begin(), y.size(), false };
  
  arma::colvec coef = arma::solve(xMatrix, yVector);
  
  return List::create(coef);
<<<<<<< HEAD
}


arma::colvec fastLm2(arma::colvec y, arma::mat x)
{
  return  arma::solve(x, y);;
}


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

=======
}
>>>>>>> b5df29349f4b8889e2d5ea5910a895694b4e9768
