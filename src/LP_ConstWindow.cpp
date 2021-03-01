// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <cmath>
#include "DCSmooth.h"
#include "DCSmooth_types.h"

using namespace Rcpp;

//---------------------------------------------------------------------------//

// function smoothes over the rows of a matrix yMat, conditional on columns

// [[Rcpp::export]]
arma::mat LPSmooth_matrix(const arma::mat yMat, const double h,
                          const int polyOrder, const int drv, SEXP kernFcnPtr)
{
  int nRow{ yMat.n_rows };    // number of conditional Time-Series 
  int nCol{ yMat.n_cols };    // number of observations per Time-Series
  int bndw{ std::max(static_cast<int>(h * nCol), polyOrder + 1) };
                              // calculate absolute bandwidth
  int windowWidth{ std::min(2*bndw + 1, nCol) };  // width of estimation window
  arma::mat yMatOut(nRow, nCol);  // matrix for results
  
  // enable Kernel function
  XPtr<funcPtr> xpfun(kernFcnPtr);
  funcPtr kernFcn = *xpfun;

  // calculate weights for interior smoothing
  arma::colvec  uVec{ arma::regspace(-bndw, bndw)/std::max(h * nCol + 1,
                                     static_cast<double>(polyOrder + 1)) }; // vector from -1 to 1 to compute weights
  arma::colvec  xVec{ arma::regspace(-bndw, bndw)/std::max(h * nCol,
                                     static_cast<double>(polyOrder + 1)) }; // vector from -1 to 1 to compute weights
  arma::colvec  weightsVec{ kernFcn(uVec, 1) };       // computation of weights

  if (h < 0.5)
  {
    // smothing over interior values
    arma::mat yInterior(nRow, windowWidth); // empty vector for use inside loop
    arma::mat coefMat(polyOrder, windowWidth);    // empty matrix for lm results
    arma::mat xMatInterior{ xMatrix(xVec, polyOrder) };  // compute x-matrix for lm-regression
    arma::mat xMatWeight{ weightMatrix(weightsVec, xMatInterior) };
    arma::mat xMatSolved{ arma::inv(xMatWeight.t() * xMatInterior)
      * xMatWeight.t() };         // compute inv(X*W*X)^(-1)*X*W only once here
    arma::rowvec xWeightsVec{ factorialFunction(drv)
      * xMatSolved.row(drv) / pow(h, drv) };

    // Loops smooth over the columns, conditional on rows. That is, every row is
    // consiedered to be an individual time series. To speed up computation, the
    // smoothing order is inverted (computing of weights only once per column, as
    // the weights are the same on a grid)
    for (int colIndex{ bndw }; colIndex < (nCol - bndw); ++colIndex)                  // outer loop over columns
    {
      yInterior  = yMat.cols(colIndex - bndw, colIndex + bndw);
      yMatOut.col(colIndex) = yInterior * xWeightsVec.t();
    }
  }

  // smoothing over boundaries
  arma::colvec  xBound(arma::regspace(0, windowWidth - 1));

  //declare empty matrices/vectors for linear regression in loop
  arma::colvec  uBound(windowWidth);
  arma::mat     xMatBound(windowWidth, polyOrder + 1);
  arma::mat     xWeightBound(windowWidth, polyOrder + 1);
  arma::mat     xMatSolved(windowWidth, polyOrder + 1);
  arma::colvec  weightsBound(windowWidth);
  arma::colvec  yLeft(windowWidth);
  arma::colvec  yRight(windowWidth);
  arma::mat     coefMatrix(polyOrder, 1);

  for (int colIndex{ 0 }; colIndex < bndw; ++colIndex)
  {
    // calculate weights for linear regression
    uBound        = (colIndex - xBound)/(windowWidth - colIndex);
    weightsBound  = kernFcn(uBound, 1);
    xMatBound     = xMatrix((xBound - colIndex)/bndw, polyOrder);
    xWeightBound  = weightMatrix(weightsBound, xMatBound);
    xMatSolved    = arma::inv(xWeightBound.t() * xMatBound)
                     * xWeightBound.t();
    arma::rowvec xWeightsLeft{ factorialFunction(drv)
      * xMatSolved.row(drv) / pow(h, drv) };
    // drv^(-1) ensures the correct sign
    arma::rowvec xWeightsRight{ pow(-1, drv) * xWeightsLeft };

    // calculation of estimates (complete column)
    yMatOut.col(colIndex) = yMat.cols(0, xBound.n_rows - 1) * xWeightsLeft.t();
    yMatOut.col(nCol - colIndex - 1) = arma::reverse(yMat.cols(nCol -
      xBound.n_rows, nCol - 1), 1) * xWeightsRight.t();
  }
  return yMatOut;
}

//---------------------------------------------------------------------------//

// [[Rcpp::export]]
arma::mat LP_DoubleSmooth(arma::mat yMat, arma::colvec hVec,
                       arma::icolvec polyOrderVec, arma::icolvec drvVec,
                       SEXP kernFcnPtr)
{
  // Smoothing over cond. on rows first (e.g. over single days).
  // Thus, drv and order is (1) instead of (0) here (depending on t)
  arma::mat mMatTemp{ LPSmooth_matrix(yMat, hVec(1),
                                      polyOrderVec(1), drvVec(1), kernFcnPtr) };
  // Smoothing over cols, drv and order is (0) (depending on x)
  arma::mat yMatOut{ LPSmooth_matrix(mMatTemp.t(), hVec(0),
                                      polyOrderVec(0), drvVec(0), kernFcnPtr) };
  
  return yMatOut.t();
}