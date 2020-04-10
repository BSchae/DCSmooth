// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <cmath>
#include "DCSmooth.h"

using namespace Rcpp;

//---------------------------------------------------------------------------//

int factorialFunction(int value);

//---------------------------------------------------------------------------//

// [[Rcpp::export]]
arma::mat weightMatrix(arma::colvec weights, arma::mat matrix)
{
  arma::mat matrixOut{ arma::mat(matrix.n_rows, matrix.n_cols) };
  for (int j{ 0 }; j < matrix.n_cols; ++j)
  {
    matrixOut.col(j) = weights % matrix.col(j);
  }

  return matrixOut;
}

//---------------------------------------------------------------------------//

double solveCoefs(arma::colvec y, arma::mat xMat, int drv)
{
  arma::mat returnMatrix{ arma::inv(xMat.t() * xMat) * xMat.t()
          * y };
  return returnMatrix(drv, 0);
}

//---------------------------------------------------------------------------//

// [[Rcpp::export]]
arma::colvec LPSmooth_grid(const arma::colvec y,
                      const double h, const int polyOrder, const int drv)
{
  int n{ y.n_rows };                      // number of observations
  int bndw{ static_cast<int>(h * n) };    // calculate absolute bandwidth, decimals will be dumped
  int windowWidth{ 2*bndw + 1 };          // width of estimation window

  arma::colvec yOut(n);               // vector for results

// smoothing over interior values
  arma::colvec  uInterior{ arma::linspace(-bndw, bndw, windowWidth)/bndw}; //
  arma::colvec  weightsInterior{ sqrt(pow(1 - pow(uInterior, 2), 2)) };
  arma::colvec  yInterior{ arma::zeros(windowWidth) };                         // empty vector for use in loop
  arma::mat     xMatInterior{ xMatrix(arma::linspace(-bndw, bndw, 
                windowWidth)/bndw, polyOrder) };
  xMatInterior = weightMatrix(weightsInterior, xMatInterior);

  for (int index{ bndw }; index < (n - bndw); ++index)
  {
    yInterior   = weightsInterior % y.subvec(index - bndw, index + bndw);
    yOut(index) = solveCoefs(yInterior, xMatInterior, drv);                    // speed up code by inserting 'solveCoefs' directly
  }


// smoothing over boundaries
  // initialise empty variables for data inside the loop
  arma::colvec  xBound(arma::linspace(0, windowWidth - 1, windowWidth));
  arma::colvec  uBound(windowWidth);
  arma::mat     xMatBound(windowWidth, polyOrder + 1);
  arma::colvec  weightsBound(windowWidth);
  arma::colvec  yLeft(windowWidth);
  arma::colvec  yRight(windowWidth);

  for (int index{ 0 }; index < bndw; ++index)
  {
    uBound        = (index - xBound)/(2 * bndw - index);
    weightsBound  = (pow(1 - pow(uBound, 2), 2));
    yLeft         = weightsBound % y.subvec(0, windowWidth - 1);
    yRight        = weightsBound % reverse(y.subvec(n - windowWidth, n - 1));
    xMatBound     = xMatrix((xBound - index)/bndw, polyOrder);
    xMatBound     = weightMatrix(weightsBound, xMatBound);

    yOut(index)   = solveCoefs(yLeft, xMatBound, drv);
    yOut(n - index - 1) = pow(-1, drv) * solveCoefs(yRight, xMatBound, drv); // pow... ensures the correct sign
  }

  return factorialFunction(drv) * yOut / pow(h, drv);
}