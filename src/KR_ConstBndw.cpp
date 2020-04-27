// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include "DCSmooth_types.h"

using namespace Rcpp;

arma::vec kernFkt_MW200(const arma::vec&, double);
arma::vec kernFkt_MW210(const arma::vec&, double);
arma::vec kernFkt_MW220(const arma::vec&, double);
arma::vec kernFkt_MW320(const arma::vec&, double);
arma::vec kernFkt_MW420(const arma::vec&, double);
arma::vec kernFkt_MW422(const arma::vec&, double);

//---------------------------------------------------------------------------//
//                    Kernel Regression Functions                            //
//---------------------------------------------------------------------------//

// function smoothes over the rows of a matrix yMat, conditional on columns

// [[Rcpp::export]]
arma::mat KRSmooth_matrix2(const arma::mat yMat, const double h,
                         SEXP kernFcnPtr) //arma::vec (*kernFktPtr)(const arma::vec&, double))
{
  int nRow{ yMat.n_rows };                // number of conditional Time-Series
  int nCol{ yMat.n_cols };                // number of observations per Time-Series
  int bndw{ static_cast<int>(h * nCol) }; // calculate absolute bandwidth, decimals will be dumped
  int windowWidth{ 2*bndw + 1 };          // width of estimation window

  arma::mat yMatOut(nRow, nCol);          // matrix for results

  // enable Kernel function
  XPtr<funcPtr> xpfun(kernFcnPtr);
  funcPtr kernFcn = *xpfun;

  // smoothing over interior values
  arma::colvec  uInterior{ arma::linspace(-bndw, bndw, windowWidth)/(h * nCol) };   // vector from -1 to 1 to compute weights
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
              SEXP kernFcnPtrX, SEXP kernFcnPtrT)
{
  arma::mat mMatTemp{ KRSmooth_matrix2(yMat, hVec(0),
                        kernFcnPtrX) };
  arma::mat yMatOut{ KRSmooth_matrix2(mMatTemp.t(), hVec(1),
                        kernFcnPtrT) };

  return yMatOut.t();
}