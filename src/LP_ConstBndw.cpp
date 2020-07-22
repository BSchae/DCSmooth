// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <cmath>
#include "DCSmooth.h"
#include "DCSmooth_types.h"

using namespace Rcpp;

arma::vec kernFkt_MW200(arma::vec&, double);
arma::vec kernFkt_MW210(arma::vec&, double);
arma::vec kernFkt_MW220(arma::vec&, double);
arma::vec kernFkt_MW320(arma::vec&, double);
arma::vec kernFkt_MW420(arma::vec&, double);
arma::vec kernFkt_MW422(arma::vec&, double);

//---------------------------------------------------------------------------//

arma::mat weightMatrix(arma::colvec weights, arma::mat matrix);

//---------------------------------------------------------------------------//

// rewrite x-Vector as x-Matrix for lm model, x can be a vector before this function
arma::mat xMatrix(arma::colvec xVector, int polyOrder);

//---------------------------------------------------------------------------//

int factorialFunction(int value);

//---------------------------------------------------------------------------//

// [[Rcpp::export]]
arma::mat LPSmooth_matrix2(const arma::mat yMat, const double h,
                           const int polyOrder, const int drv, SEXP kernFcnPtr)
{
  int nRow{ yMat.n_rows };                // number of conditional Time-Series
  int nCol{ yMat.n_cols };                // number of observations per Time-Series
  int bndw{ std::max(static_cast<int>(h * nCol), polyOrder + 1) }; // calculate absolute bandwidth, decimals will be dumped
  int windowWidth{ 2*bndw + 1 };  // width of estimation window
  arma::mat yMatOut(nRow, nCol);          // matrix for results
  
  // enable Kernel function
  XPtr<funcPtr> xpfun(kernFcnPtr);
  funcPtr kernFcn = *xpfun;
  
  // calculate weights
  arma::colvec  uVec{ arma::regspace(-bndw, bndw)/std::max(h * nCol, 
                                     static_cast<double>(polyOrder + 1)) }; // vector from -1 to 1 to compute weights
  arma::colvec  weightsVec{ kernFcn(uVec, 1) };         // computation of weights (put in kernel-function later)
  
  // smoothing over boundaries
  for (int colIndex{ 0 }; colIndex < bndw; ++colIndex)
  {
    if (colIndex == (nCol/2) + 1)     // break condition, if bndw > 0.5*nCol
    {
      break;
    }
    
    int startIndex{ bndw - colIndex };
    int stopIndex{ std::min(windowWidth - 1, startIndex + nCol - 1) };

    arma::colvec weightsBound{ weightsVec.subvec(startIndex, stopIndex) };
    arma::colvec uBound{ uVec.subvec(startIndex, stopIndex) };
    arma::mat     xMatBound{ xMatrix(uBound, polyOrder) }; // put this two lines in one function???
    xMatBound     = weightMatrix(weightsBound, xMatBound);
    //arma::colvec  yLeft(uBound.n_rows);
    //arma::colvec  yRight(uBound.n_rows);
    //arma::mat     coefMatrix(polyOrder, 1);
    arma::mat     xMatSolved{ arma::inv(xMatBound.t() * xMatBound)
      * xMatBound.t() };
    
    // ??? get rid of the .t() maybe?
    arma::mat coefLeft  = xMatSolved * diagmat(weightsBound) 
                  * yMat.cols(0, uBound.n_rows - 1).t();
    arma::mat coefRight = xMatSolved * diagmat(weightsBound) 
                  * reverse(yMat.cols(nCol - uBound.n_rows, nCol - 1), 1).t();
    yMatOut.col(colIndex) = pow(-1, drv) * coefLeft.row(drv).t();
    yMatOut.col(nCol - colIndex - 1) = pow(-1, drv) *   coefRight.row(drv).t();
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

//---------------------------------------------------------------------------//

// [[Rcpp::export]]
arma::mat LP_DoubleSmooth2(arma::mat yMat, arma::colvec hVec,
                       arma::icolvec polyOrderVec, arma::icolvec drvVec,
                       SEXP kernFcnPtr)
{
  arma::mat mMatTemp{ LPSmooth_matrix2(yMat, hVec(1),
                                    polyOrderVec(1), drvVec(1), kernFcnPtr) };
  arma::mat yMatOut{ LPSmooth_matrix2(mMatTemp.t(), hVec(0),
                                    polyOrderVec(0), drvVec(0), kernFcnPtr) };

  return yMatOut.t();
}