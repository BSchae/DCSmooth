// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <cmath>
#include "DCSmooth.h"

using namespace Rcpp;

arma::mat weightMatrix(arma::colvec weights, arma::mat matrix);

//---------------------------------------------------------------------------//

// [[Rcpp::export]]
int factorialFunction(int value)
{
  int outValue{ 1 };

  for(int count{ 1 }; count <= value; ++count)
    outValue *= count;

  return outValue;
}

//---------------------------------------------------------------------------//

// function smoothes over the rows of a matrix yMat, conditional on columns

// [[Rcpp::export]]
arma::mat LPSmooth_matrix(const arma::mat yMat, const double h,
                          const int polyOrder, const int drv)
{
  int nRow{ yMat.n_rows };                // number of conditional Time-Series
  int nCol{ yMat.n_cols };                // number of observations per Time-Series
  int bndw{ static_cast<int>(h * nCol) }; // calculate absolute bandwidth, decimals will be dumped
  int windowWidth{ 2*bndw + 1 };          // width of estimation window
  
  arma::mat yMatOut(nRow, nCol);          // matrix for results
  
// smoothing over interior values
  arma::colvec  uInterior{ arma::linspace(-bndw, bndw, windowWidth)/(h * nCol) };   // vector from -1 to 1 to compute weights
  arma::colvec  weightsInterior{ (0.75*pow(1 - pow(uInterior, 2), 2)) };            // computation of weights (put in kernel-function later)
  arma::colvec  yRowInterior{ arma::zeros(windowWidth) };                           // empty vector for use inside loop
  arma::mat     coefMat(polyOrder, 1);                                              // empty matrix for lm results
  arma::mat     xMatInterior{ xMatrix(arma::linspace(-bndw, bndw,
                                    windowWidth)/bndw, polyOrder) };                // compute x-matrix for lm-regression
  xMatInterior  = weightMatrix(weightsInterior, xMatInterior);
  arma::mat     xMatSolved{ arma::inv(xMatInterior.t() * xMatInterior)
                          * xMatInterior.t() };                                     // compute inv(X*X)*X only once here

  // Loops smooth over the columns, conditional on rows. That is, every row is
  // consiedered to be an individual time series. To speed up computation, the
  // smoothing order is inverted (computing of weights only once per column, as
  // the weights are the same on a grid)
  for (int colIndex{ bndw }; colIndex < (nCol - bndw); ++colIndex)                  // outer loop over columns
  {
    for (int rowIndex{ 0 }; rowIndex < nRow; ++rowIndex)                            // inner loop over rows
    {
      yRowInterior  = yMat(rowIndex, arma::span(colIndex - bndw, colIndex + bndw)).t();
      coefMat = xMatSolved * (weightsInterior % yRowInterior);
      yMatOut(rowIndex, colIndex) = coefMat(drv, 0);
    }
  }
 
// smoothing over boundaries
  arma::colvec  xBound(arma::linspace(0, windowWidth - 1, windowWidth));
  arma::colvec  uBound(windowWidth);
  arma::mat     xMatBound(windowWidth, polyOrder + 1);
  arma::colvec  weightsBound(windowWidth);
  arma::colvec  yLeft(windowWidth);
  arma::colvec  yRight(windowWidth);
  arma::mat     coefMatrix(polyOrder, 1);
  
  for (int colIndex{ 0 }; colIndex < bndw; ++colIndex)
  {
    uBound        = (colIndex - xBound)/(2 * bndw - colIndex);
    weightsBound  = (pow(1 - pow(uBound, 2), 2));
    xMatBound     = xMatrix((xBound - colIndex)/bndw, polyOrder);
    xMatBound     = weightMatrix(weightsBound, xMatBound);
    xMatSolved    = arma::inv(xMatBound.t() * xMatBound)
                    * xMatBound.t();

    for (int rowIndex{ 0 }; rowIndex < nRow; ++rowIndex)
    {
      yLeft       = weightsBound % yMat(rowIndex, arma::span(0, 
                                                windowWidth - 1)).t();
      yRight      = weightsBound % reverse(yMat(rowIndex, arma::span(nCol 
                                              - windowWidth, nCol - 1))).t();
      
      coefMatrix  = xMatSolved * yLeft;
      yMatOut(rowIndex, colIndex) = coefMatrix(drv, 0);
      
      coefMatrix  = xMatSolved * yRight;
      yMatOut(rowIndex, nCol - colIndex - 1) = pow(-1, drv) * coefMatrix(drv, 0); // pow(-1, drv) ensures the correct sign
    }
  }
  
  return factorialFunction(drv) * yMatOut / pow(h, drv);
}

//---------------------------------------------------------------------------//

// [[Rcpp::export]]
arma::mat FastDoubleSmooth(arma::mat yMat, arma::colvec hVec,
                       arma::icolvec polyOrderVec, arma::icolvec drvVec)
{
  arma::mat mMatTemp{ LPSmooth_matrix(yMat, hVec(0),
                                      polyOrderVec(0), drvVec(0)) };
  arma::mat yMatOut{ LPSmooth_matrix(mMatTemp.t(), hVec(1),
                                      polyOrderVec(1), drvVec(1)) };
  
  return yMatOut.t();
}



