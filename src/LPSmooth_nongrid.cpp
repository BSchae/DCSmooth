// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <cmath>

using namespace Rcpp;

//---------------------------------------------------------------------------//

// Multiple linear regression, where y is the vector of endogenous variables
// and x is the matrix of exogenous variables. w is the diagonal of the
// weighting matrix. The function is progrtammed in c++/armadillo.
arma::mat lmFunction(arma::colvec y, arma::mat x, arma::colvec w)
{
  // put vector of weights into diagonal matrix
  arma::mat weights{ arma::diagmat(w) };

  // calculate weighted variables
  arma::mat xWeighted{ weights * x };
  arma::colvec yWeighted{ weights * y };

  // do and return coefficients of mult. linear regression
  return arma::solve(xWeighted, yWeighted);
}

//---------------------------------------------------------------------------//

// rewrite x-Vector as x-Matrix for lm model, x can be a vector before this function
arma::mat xMatrix(arma::colvec xVector, int polyOrder)
{
  arma::mat returnMatrix{ arma::ones(xVector.n_rows, polyOrder + 1) };

  // put powers of x into matrix depending on order (local linear: order = 1)
  for (int indexOrder{ 0 }; indexOrder < polyOrder; ++indexOrder)
    returnMatrix.col(indexOrder + 1) = pow(xVector, indexOrder + 1);

  return returnMatrix;
}

//---------------------------------------------------------------------------//

// [[Rcpp::export]]
arma::colvec LPSmooth_nongrid(const arma::colvec y, const arma::colvec x,
                                  const double h, const int polyOrder)
{
  // vector for results
  arma::colvec yOut(y.n_rows);

  for (int index{ 0 }; index < y.n_rows; ++index)
  {
    // define subvectors around evaluation point x0
    double x0{ x(index) };
    arma::uvec   subIndex{ find(abs(x - x0) < h) };
    arma::colvec xSub{ x(subIndex) - x0 };
    arma::colvec ySub{ y(subIndex) };
    arma::colvec u{ xSub/h };
    arma::colvec w{ pow(1 - pow(u, 2), 2) };
    
    arma::mat xMat{ xMatrix(xSub, polyOrder) };
    arma::mat coefMat{ lmFunction(ySub, xMat, w) };
    yOut(index) = coefMat(0, 0);
  }

  return yOut;
}