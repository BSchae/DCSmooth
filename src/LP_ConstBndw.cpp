// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <cmath>
#include "DCSmooth.h"

using namespace Rcpp;

//---------------------------------------------------------------------------//

arma::mat weightMatrix(arma::colvec weights, arma::mat matrix);

//---------------------------------------------------------------------------//

// rewrite x-Vector as x-Matrix for lm model, x can be a vector before this function
arma::mat xMatrix(arma::colvec xVector, int polyOrder);

//---------------------------------------------------------------------------//

int factorialFunction(int value);

//---------------------------------------------------------------------------//

// [[Rcpp::export]]
arma::mat LPSmooth_matrix3(const arma::mat yMat, const double h,
                           const int polyOrder, const int drv)
{
  int nRow{ yMat.n_rows };                // number of conditional Time-Series
  int nCol{ yMat.n_cols };                // number of observations per Time-Series
  int bndw{ std::max(static_cast<int>(h * nCol), polyOrder + 1) }; // calculate absolute bandwidth, decimals will be dumped
  int windowWidth{ 2*bndw + 1 };  // width of estimation window
  
  //return bndw;
  
  arma::mat yMatOut(nRow, nCol);          // matrix for results
  
  // calculate weights
  arma::colvec  uVec{ arma::regspace(-bndw, bndw)/std::max(h * nCol, 
                                     static_cast<double>(polyOrder + 1)) }; // vector from -1 to 1 to compute weights
  arma::colvec  weightsVec{ (0.9375*pow(1 - pow(uVec, 2), 2)) };         // computation of weights (put in kernel-function later)
  
  // smoothing over boundaries
  for (int colIndex{ 0 }; colIndex < bndw; ++colIndex)
  {
    if (colIndex == (nCol/2) + 1)     // break condition, if bndw > 0.5*nCol
    {
      break;
    }
    
    int startIndex{ (windowWidth - 1)/2 - colIndex };
    int stopIndex{ startIndex + std::min(bndw, nCol - 1) };
    
    arma::colvec weightsBound{ weightsVec.subvec(startIndex, stopIndex) };
    arma::colvec uBound{ uVec.subvec(startIndex, stopIndex) };
    arma::mat     xMatBound{ xMatrix(uBound, polyOrder) }; // put this two lines in one function???
    xMatBound     = weightMatrix(weightsBound, xMatBound);
    //arma::colvec  yLeft(uBound.n_rows);
    //arma::colvec  yRight(uBound.n_rows);
    arma::mat     coefMatrix(polyOrder, 1);
    arma::mat     xMatSolved{ arma::inv(xMatBound.t() * xMatBound)
      * xMatBound.t() };
    
    // ??? get rid of the .t() maybe?
    arma::mat coefLeft  = xMatSolved * diagmat(weightsBound) 
                  * yMat.cols(0, uBound.n_rows - 1).t();
    arma::mat coefRight = xMatSolved * diagmat(weightsBound) 
                  * reverse(yMat.cols(nCol - uBound.n_rows, nCol - 1), 1).t();
    yMatOut.col(colIndex) = coefLeft.row(drv).t();
    yMatOut.col(nCol - colIndex - 1) = pow(-1, drv) * coefRight.row(drv).t();
                                      // pow(-1, drv) ensures the correct sign
  }
  
  if (h < 0.5)
  {
    // smothing over interior values
    arma::mat     yInterior(nRow, windowWidth);                           // empty vector for use inside loop
    arma::mat     coefMat(polyOrder, windowWidth);                                              // empty matrix for lm results
    arma::mat     xMatInterior{ xMatrix(uVec, polyOrder) };                // compute x-matrix for lm-regression
    xMatInterior  = weightMatrix(weightsVec, xMatInterior);
    arma::mat     xMatSolved{ arma::inv(xMatInterior.t() * xMatInterior)
      * xMatInterior.t() };                                     // compute inv(X*X)*X only once here
    
    // Loops smooth over the columns, conditional on rows. That is, every row is
    // consiedered to be an individual time series. To speed up computation, the
    // smoothing order is inverted (computing of weights only once per column, as
    // the weights are the same on a grid)
    for (int colIndex{ bndw }; colIndex < (nCol - bndw); ++colIndex)                  // outer loop over columns
    {
      yInterior  = yMat.cols(colIndex - bndw, colIndex + bndw);
      coefMat = (yInterior * diagmat(weightsVec)) * xMatSolved.t();
      yMatOut.col(colIndex) = coefMat.col(drv);
    }
  }
  
  return factorialFunction(drv) * yMatOut / pow(h, drv);
}


// function smoothes over the rows of a matrix yMat, conditional on columns

// this is only for testing purposes.

// [[Rcpp::export]]
arma::mat LPSmooth_matrix2(const arma::mat yMat, const double h,
                          const int polyOrder, const int drv)
{
  int nRow{ yMat.n_rows };                // number of conditional Time-Series
  int nCol{ yMat.n_cols };                // number of observations per Time-Series
  int bndw{ std::max(static_cast<int>(h * nCol), polyOrder + 1) }; // calculate absolute bandwidth, decimals will be dumped
  int windowWidth{ 2*bndw + 1 };  // width of estimation window
  
  //return bndw;
  
  arma::mat yMatOut(nRow, nCol);          // matrix for results

  // calculate weights
  arma::colvec  uVec{ arma::regspace(-bndw, bndw)/std::max(h * nCol, 
                                  static_cast<double>(polyOrder + 1)) }; // vector from -1 to 1 to compute weights
  arma::colvec  weightsVec{ (0.9375*pow(1 - pow(uVec, 2), 2)) };         // computation of weights (put in kernel-function later)

  // smoothing over boundaries
  for (int colIndex{ 0 }; colIndex < bndw; ++colIndex)
    {
      if (colIndex == (nCol/2) + 1)     // break condition, if bndw > 0.5*nCol
      {
        break;
      }

      int startIndex{ (windowWidth - 1)/2 - colIndex };
      int stopIndex{ startIndex + std::min(bndw, nCol - 1) };

      arma::colvec weightsBound{ weightsVec.subvec(startIndex, stopIndex) };
      arma::colvec uBound{ uVec.subvec(startIndex, stopIndex) };
      arma::mat     xMatBound{ xMatrix(uBound, polyOrder) }; // put this two lines in one function???
      xMatBound     = weightMatrix(weightsBound, xMatBound);
      arma::colvec  yLeft(uBound.n_rows);
      arma::colvec  yRight(uBound.n_rows);
      arma::mat     coefMatrix(polyOrder, 1);
      arma::mat     xMatSolved{ arma::inv(xMatBound.t() * xMatBound)
                                * xMatBound.t() };

      for (int rowIndex{ 0 }; rowIndex < nRow; ++rowIndex)
      {
        yLeft       = weightsBound % yMat(rowIndex, arma::span(0,
                                                  uBound.n_rows - 1)).t();
        yRight      = weightsBound % reverse(yMat(rowIndex, arma::span(nCol
                                                - uBound.n_rows, nCol - 1))).t();

        coefMatrix  = xMatSolved * yLeft;
        yMatOut(rowIndex, colIndex) = coefMatrix(drv, 0);

        coefMatrix  = xMatSolved * yRight;
        yMatOut(rowIndex, nCol - colIndex - 1) = pow(-1, drv) * coefMatrix(drv, 0); // pow(-1, drv) ensures the correct sign
      }
  }

  if (h < 0.5)
  {
    // smothing over interior values
    arma::colvec  yRowInterior{ arma::zeros(windowWidth) };                           // empty vector for use inside loop
    arma::mat     coefMat(polyOrder, 1);                                              // empty matrix for lm results
    arma::mat     xMatInterior{ xMatrix(uVec, polyOrder) };                // compute x-matrix for lm-regression
    xMatInterior  = weightMatrix(weightsVec, xMatInterior);
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
        coefMat = xMatSolved * (weightsVec % yRowInterior);
        yMatOut(rowIndex, colIndex) = coefMat(drv, 0);
      }
    }
  }

  return factorialFunction(drv) * yMatOut / pow(h, drv);
}

//---------------------------------------------------------------------------//

// [[Rcpp::export]]
arma::mat LP_DoubleSmooth2(arma::mat yMat, arma::colvec hVec,
                       arma::icolvec polyOrderVec, arma::icolvec drvVec)
{
  arma::mat mMatTemp{ LPSmooth_matrix2(yMat, hVec(1),
                                      polyOrderVec(1), drvVec(1)) };
  arma::mat yMatOut{ LPSmooth_matrix2(mMatTemp.t(), hVec(0),
                                      polyOrderVec(0), drvVec(0)) };

  return yMatOut.t();
}