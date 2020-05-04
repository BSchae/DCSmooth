// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include "DCSmooth_types.h"

using namespace Rcpp;

arma::vec kernFkt_MW200(arma::vec&, double);
arma::vec kernFkt_MW210(arma::vec&, double);
arma::vec kernFkt_MW220(arma::vec&, double);
arma::vec kernFkt_MW320(arma::vec&, double);
arma::vec kernFkt_MW420(arma::vec&, double);
arma::vec kernFkt_MW422(arma::vec&, double);


// Test
// [[Rcpp::export]]
arma::mat KRTest(arma::mat yMat, double h,
  SEXP kernFcnPtr)
{
  int nRow{ yMat.n_rows };                // number of conditional Time-Series
  int nCol{ yMat.n_cols };                // number of observations per Time-Series
  int bndw{ std::max(static_cast<int>(h * nCol), 2) }; // calculate absolute bandwidth, decimals will be dumped
  int windowWidth{ 2*bndw + 1 };          // width of estimation window

  arma::mat yMatOut(nRow, nCol);          // matrix for results

  //enable Kernel function
  XPtr<funcPtr> xpfun(kernFcnPtr);
  funcPtr kernFcn = *xpfun;
  
  // smoothing over boundaries
  for (int colIndex{ 0 }; colIndex < bndw; ++colIndex)
  {
    double q = static_cast<double>(colIndex)/bndw;
    arma::colvec  uBound(arma::regspace(colIndex, -bndw) / bndw);
    arma::colvec  weightsBound{ kernFcn(uBound, q)*(1 + q)/uBound.n_rows };
    
    arma::mat     yLeftMat{ yMat.cols(0, uBound.n_rows - 1) };
    arma::mat     yRightMat{ yMat.cols(nCol - uBound.n_rows, nCol - 1) };
    
    yMatOut.col(colIndex) = yLeftMat * weightsBound;
    yMatOut.col(nCol - colIndex - 1) = yRightMat * reverse(weightsBound);
    
  }
  return yMatOut;
}




















//---------------------------------------------------------------------------//
//                    Kernel Regression Functions                            //
//---------------------------------------------------------------------------//

// function smoothes over the rows of a matrix yMat, conditional on columns

// [[Rcpp::export]]
arma::mat KRSmooth_matrix2(arma::mat yMat, double h,
                         SEXP kernFcnPtr) //arma::vec (*kernFktPtr)(const arma::vec&, double))
{
  int nRow{ yMat.n_rows };                // number of conditional Time-Series
  int nCol{ yMat.n_cols };                // number of observations per Time-Series
  int bndw{ std::max(static_cast<int>(h * nCol), 2) }; // calculate absolute bandwidth, decimals will be dumped
  int windowWidth{ 2*bndw + 1 };          // width of estimation window

  arma::mat yMatOut(nRow, nCol);          // matrix for results

  // enable Kernel function
  XPtr<funcPtr> xpfun(kernFcnPtr);
  funcPtr kernFcn = *xpfun;

  // smoothing over interior values
  arma::colvec  uInterior{ arma::linspace(-bndw, bndw, windowWidth)/bndw }; //(h * nCol) };   // vector from -1 to 1 to compute weights
  arma::colvec  weightsInterior{ kernFcn(uInterior, 1)/bndw };            // computation of weights (put in kernel-function later)
  arma::colvec  yRowInterior{ arma::zeros(windowWidth) };                           // empty vector for use inside loop

  // Loops smooth over the columns, conditional on rows. That is, every row is
  // consiedered to be an individual time series. To speed up computation, the
  // smoothing order is inverted (computing of weights only once per column, as
  // the weights are the same on a grid)
  for (int colIndex{ bndw }; colIndex < (nCol - bndw); ++colIndex)                  // outer loop over columns
  {
    for (int rowIndex{ 0 }; rowIndex < nRow; ++rowIndex)                            // inner loop over rows
    {
      yRowInterior  = yMat(rowIndex, arma::span(colIndex - bndw, colIndex + bndw)).t();
      yMatOut(rowIndex, colIndex) = sum(weightsInterior % yRowInterior);
    }
  }

  // smoothing over boundaries
  for (int colIndex{ 0 }; colIndex < bndw; ++colIndex)
  {
    double q = static_cast<double>(colIndex)/bndw;
    arma::colvec  uBound(arma::regspace(colIndex, -bndw) / bndw);
    arma::colvec  weightsBound{ kernFcn(uBound, q)/bndw };
    arma::colvec  yLeft(uBound.n_rows);
    arma::colvec  yRight(uBound.n_rows);

    for (int rowIndex{ 0 }; rowIndex < nRow; ++rowIndex)
    {
      yLeft       = weightsBound % yMat(rowIndex, arma::span(0,
                            uBound.n_rows - 1)).t();
      yRight      = weightsBound % reverse(yMat(rowIndex, arma::span(nCol
                            - uBound.n_rows, nCol - 1))).t();

      yMatOut(rowIndex, colIndex) = sum(yLeft);
      yMatOut(rowIndex, nCol - colIndex - 1) = sum(yRight);
    }
  }
  return yMatOut;
}

//---------------------------------------------------------------------------//

// [[Rcpp::export]]
arma::mat KR_DoubleSmooth2(arma::mat yMat, arma::colvec hVec,
                    arma::colvec drvVec, SEXP kernFcnPtrX, SEXP kernFcnPtrT)
{
  // Smoothing over cond. on rows first (e.g. over single days).
  // Thus, drv and order is (1) instead of (0) here (depending on t)
  arma::mat mMatTemp{ KRSmooth_matrix2(yMat, hVec(1),
                        kernFcnPtrT)/pow(hVec(1), drvVec(1)) };
  // Smoothing over cols, drv and order is (0) (depending on x)
  arma::mat yMatOut{ KRSmooth_matrix2(mMatTemp.t(), hVec(0),
                        kernFcnPtrX)/pow(hVec(0), drvVec(0)) };

  return yMatOut.t();
}