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

//---------------------------------------------------------------------------//
//                    Kernel Regression Functions                            //
//---------------------------------------------------------------------------//

// function smoothes over the rows of a matrix yMat, conditional on columns

// [[Rcpp::export]]
arma::mat KRSmooth_matrix(arma::mat yMat, double h, int drv,
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
  arma::colvec  uInterior{ arma::linspace(-bndw, bndw, windowWidth)/(h * nCol) };   // vector from -1 to 1 to compute weights
  arma::colvec  weightsInterior{ kernFcn(uInterior, 1)/bndw };            // computation of weights (put in kernel-function later)
  arma::mat     yMatInterior(nRow, windowWidth + 1);                           // empty matrix for use inside loop

  // Loops smooth over the columns, conditional on rows. That is, every row is
  // consiedered to be an individual time series. To speed up computation, the
  // smoothing order is inverted (computing of weights only once per column, as
  // the weights are the same on a grid)
  for (int colIndex{ bndw }; colIndex < (nCol - bndw); ++colIndex)                  // outer loop over columns
  {
    yMatInterior = yMat.cols(colIndex - bndw, colIndex + bndw);
    yMatOut.col(colIndex) = yMatInterior * weightsInterior;
  }

  // smoothing over boundaries
  arma::mat    yLeftMat{ yMat.cols(0, windowWidth - 1) };
  arma::mat    yRightMat{ yMat.cols(nCol - windowWidth, nCol - 1) };
  arma::colvec uBound(windowWidth);
  arma::colvec weightsBound(windowWidth);

  for (int colIndex{ 0 }; colIndex < bndw; ++colIndex)
  {
    double q = (colIndex - 1.0)/(2*bndw - colIndex + 1.0);
    uBound = arma::linspace(q, -1, windowWidth);

    weightsBound = kernFcn(uBound, q);
    // if (drv ==  0)
    // {
    //   weightsBound = weightsBound/sum(weightsBound);
    // } else {
      weightsBound = weightsBound * (1.0 + q)/static_cast<double>(windowWidth - 1);
    // }
    yMatOut.col(colIndex) = yLeftMat * weightsBound;
    yMatOut.col(nCol - colIndex - 1) = yRightMat * reverse(weightsBound);
  }
  
  return yMatOut;
}

//---------------------------------------------------------------------------//

// [[Rcpp::export]]
arma::mat KR_DoubleSmooth(arma::mat yMat, arma::colvec hVec,
              arma::icolvec drvVec, SEXP kernFcnPtrX, SEXP kernFcnPtrT)
{
  // Smoothing over cond. on rows first (e.g. over single days).
  // Thus, drv and order is (1) instead of (0) here (depending on t)
  arma::mat mMatTemp{ KRSmooth_matrix(yMat, hVec(1), drvVec(1), kernFcnPtrT) };
  // Smoothing over cols, drv and order is (0) (depending on x)
  arma::mat yMatOut{ KRSmooth_matrix(mMatTemp.t(), hVec(0), drvVec(1),
                                        kernFcnPtrX) };

  return yMatOut.t();
}