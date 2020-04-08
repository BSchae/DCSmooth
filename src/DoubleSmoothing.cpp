// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <cmath>
#include "DCSmooth.h"

using namespace Rcpp;

arma::colvec LPSmooth_grid(const arma::colvec y,
  const double h, const int polyOrder, const int drv);

//---------------------------------------------------------------------------//

// [[Rcpp::export]]
arma::mat DoubleSmooth(arma::mat yMat, arma::colvec hVec, 
          arma::colvec polyOrderS, arma::colvec drv)
{
  int nRows{ static_cast<int>(yMat.n_rows) };
  int nCols{ static_cast<int>(yMat.n_cols) };
  
  arma::mat mMatTemp(nRows, nCols);
  arma::mat yMatOut(nRows, nCols);
  
// smoothing over rows
  for (int rowIndex{ 0 }; rowIndex < nRows; ++rowIndex)
  {
    mMatTemp.row(rowIndex) = (LPSmooth_grid(yMat.row(rowIndex).t(), hVec(0), polyOrderS(0), drv(0))).t();
  }
  
// smoothing over colums
  for (int colIndex{ 0 }; colIndex < nCols; ++colIndex)
  {
    yMatOut.col(colIndex) = LPSmooth_grid(mMatTemp.col(colIndex), hVec(1), polyOrderS(1), drv(1));
  }
  
  return yMatOut;
}

