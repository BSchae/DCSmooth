// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <cmath>
#include "DCSmooth.h"
#include "DCSmooth_types.h"

using namespace Rcpp;

//---------------------------------------------------------------------------//

// [[Rcpp::export]]
arma::mat LPSmooth_matrix2_BMod(const arma::mat yMat, const double h,
                           const int polyOrder, const int drv, SEXP kernFcnPtr)
{
  // get additional information on nX, nT, bndw etc.
  int nRow{ yMat.n_rows };
  int nCol{ yMat.n_cols };
  int bndw{ std::max(static_cast<int>(h * nCol), polyOrder + 1) };    // calculate absolute bandwidth, decimals will be dumped
  arma::mat yMatOut(nRow, nCol);   // result matrix
  
  // enable Kernel function
  XPtr<funcPtr> xpfun(kernFcnPtr);
  funcPtr kernFcn = *xpfun;
  
  // smoothing over boundaries
  for (int colIndex{ 0 }; colIndex < bndw + 1; ++colIndex)
  {
    if (colIndex == (nCol/2) + 1)     // break condition, if bndw > 0.5*nCol
    {
      break;
    }
    
    // calculate kernel weights
    double q{ static_cast<double>(colIndex)/(bndw - 1) };
    arma::colvec xBound{ arma::regspace(-colIndex, bndw) / (nCol - 1) }; // vector for exogenous variables. is [q, -1]
    arma::colvec uBound{ - xBound / h };
    arma::colvec wBound{ (kernFcn(uBound, q)) };           // computation of weights

    // calculate regression weights for linear regression
    arma::mat    xMatBound{ xMatrix(xBound, polyOrder) };
    arma::mat    xMatWeight{ weightMatrix(wBound, xMatBound) };
    arma::mat    weightsMat{ arma::inv(xMatWeight.t() * xMatBound) *
                             xMatWeight.t() };
    
    arma::rowvec weightsLeft{ factorialFunction(drv) * weightsMat.row(drv) };
    arma::rowvec weightsRight{ pow(-1, drv) * weightsLeft }; // pow(-1, drv) ensures the correct sign

    // calculation of estimates (complete column)
    yMatOut.col(colIndex) = yMat.cols(0, xBound.n_rows - 1) * weightsLeft.t();
    yMatOut.col(nCol - colIndex - 1) = arma::reverse(yMat.cols(nCol -
                              uBound.n_rows, nCol - 1), 1) * weightsRight.t();
  }
  
  if (h < 0.5)
  {
    // calculate weights for interior smoothing
    arma::colvec xVec{ arma::regspace(-bndw, bndw) / (nCol - 1) };
    arma::colvec uVec{ - xVec / h };
    arma::colvec wVec{ kernFcn(uVec, 1) };
    
    arma::mat xMatInterior{ xMatrix(xVec, polyOrder) };
    arma::mat xWeightsInterior{ weightMatrix(wVec, xMatInterior) };
    arma::mat weightsMat{ arma::inv(xWeightsInterior.t() * xMatInterior) *
      xWeightsInterior.t() * factorialFunction(drv) };
    arma::rowvec weightsInterior{ weightsMat.row(drv) };
    
    // Loops smooth over the columns, conditional on rows. That is, every row is
    // considered to be an individual time series. To speed up computation, the
    // smoothing order is inverted (computing of weights only once per column, as
    // the weights are the same on a grid)
    for (int colIndex{ bndw + 1 }; colIndex < (nCol - bndw - 1); ++colIndex)
    {
      arma::mat yInterior{ yMat.cols(colIndex - bndw, colIndex + bndw) };
      yMatOut.col(colIndex) = yInterior * weightsInterior.t();
    }
  }
  
  return yMatOut;
}

//---------------------------------------------------------------------------//

// [[Rcpp::export]]
arma::mat LP_dcs_const0_BMod(arma::mat yMat, arma::colvec hVec,
                       arma::colvec polyOrderVec, arma::icolvec drvVec,
                       SEXP kernFcnPtr_x, SEXP kernFcnPtr_t)
{
  arma::mat mMatTemp{ LPSmooth_matrix2_BMod(yMat, hVec(1),
                                    polyOrderVec(1), drvVec(1), kernFcnPtr_t) };
  arma::mat yMatOut{ LPSmooth_matrix2_BMod(mMatTemp.t(), hVec(0),
                                    polyOrderVec(0), drvVec(0), kernFcnPtr_x) };

  return yMatOut.t();
}