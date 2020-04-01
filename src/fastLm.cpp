// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <cmath>

using namespace Rcpp;


// Multiple linear regression, where y is the vector of endogenous variables
// and x is the matrix of exogenous variables. w is the diagonal of the
// weighting matrix. The function is progrtammed in c++/armadillo.
arma::mat fastLm(arma::colvec y, arma::mat x, arma::colvec w)
{
  // put vector of weights into diagonal matrix
  arma::mat weights{ arma::diagmat(w) };
  
  // calculate weighted variables
  arma::mat xWeighted{ weights * x };
  arma::colvec yWeighted{ weights * y };
  
  // do and return coefficients of mult. linear regression
  return arma::solve(xWeighted, yWeighted);
}



// // [[Rcpp::export]]
// NumericVector LocalLinear(NumericVector y, NumericMatrix x, NumericVector kernWeights)
// {
//   // get bandwidth (kernWeights should have defined length 2*bndw + 1)
//   int bndw{ (kernWeights.size() - 1)/2 };
// 
//   NumericVector ySmooth(y.length(), 0);       // declare empty vector for results
// 
//   for (int index{ 0 }; index < y.size(); ++index)
//   {
//     double x0{ x(index, 1) }; // point of evaluation
//     NumericMatrix xSubset{ x(abs(x - x0)/bndw, _) }
//     NumericVector ySubset{ y(abs(x - x0)/bndw) };
//   //}
// 
//   return ;
// 
// 
// }

arma::mat xMatrix(arma::colvec xVector, int order)
{
  arma::mat returnMatrix{ arma::ones(xVector.n_rows, order + 1) };

  for (int indexOrder{ 0 }; indexOrder < order; ++indexOrder)
    returnMatrix.col(indexOrder + 1) = pow(xVector, indexOrder + 1);
  
  return returnMatrix;

}

// [[Rcpp::export]]
List outputTest(arma::colvec y, arma::colvec x, arma::colvec w, int order)
{
  arma::mat xMat{ xMatrix(x, order) };
  
  return List::create(fastLm(y, xMat, w));
}