// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <cmath>
#include "DCSmooth.h"

using namespace Rcpp;

arma::colvec LPSmooth_grid(const arma::colvec y,
  const double h, const int polyOrder);

//---------------------------------------------------------------------------//

arma::mat DoubleSmooth(arma::mat yMat, arma::colvec hVec, 
          arma:colvec polyOrderS, arma::colvec drv)
{
  int nRows{ yMat.n_rows };
  int nCols{ yMat.n_cols };
  
  arma::mat mMatTemp(nRows, nCols);
  arma::mat yMatOut(nRows, nCols);
  
// smoothing over rows
  for (int rowIndex{ 0 }; rowIndex < nRows; ++rowIndex)
  {
    mMatTemp.row(rowIndex) = LPSmooth(yMat.row(rowIndex), hVec(0), polyOrderS(0));
  }
  
// smoothing over colums
  for (int colIndex{ 0 }; colIndex < nCols; ++colIndex)
  {
    yMatOut.col(colIndex) = LPSmooth(mMatTemp.col(colIndex), hVec(1), polyOrderS(1));
  }
  
  return yMatOut;
  
}

