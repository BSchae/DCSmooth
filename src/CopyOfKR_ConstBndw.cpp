// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <cmath>
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
arma::rowvec KRSlow_single(arma::rowvec yVec, double h, SEXP kernFcnPtr)
{
  int n{ yVec.n_cols };                // number of conditional Time-Series

  arma::rowvec  xVec{ arma::linspace(0, 1, n).t() };
  arma::rowvec  yEst(n);          // vector for results
  
  // enable Kernel function
  XPtr<funcPtr> xpfun(kernFcnPtr);
  funcPtr kernFcn = *xpfun;

  int leftIndex{ static_cast<int>(n*(1 - h)) };

  for (int i{ 0 }; i < leftIndex; ++i)
  {
    double x0{ xVec(i) };
    arma::uvec index = find(abs(xVec - x0) < h);
    arma::vec  xIndex = x0 - xVec(index);
    arma::vec  yIndex = yVec(index);
    arma::vec  uVec{ xIndex/h };
    double     q{ (uVec.n_rows - 1)/(h*n - 1) - 1 };
    arma::vec  weights{ kernFcn(uVec, q)/(n*h) };
    yEst(i) =  sum(weights % yIndex);
  }

  for (int i{ leftIndex }; i < n; ++i)
  {
    double x0{ xVec(i) };
    arma::uvec index = find(abs(xVec - x0) < h);
    arma::vec  xIndex = xVec(index) - x0;
    arma::vec  yIndex = yVec(index);
    arma::vec  uVec{ xIndex/h };
    double     q{ (uVec.n_rows - 1)/(h*n - 1) - 1 };
    arma::vec  weights{ kernFcn(uVec, q)/(n*h) };
    yEst(i) =  sum(weights % yIndex);
  }

  return yEst;
}

// KRSlow_single(Y, 0.1, kernFcn0)

//---------------------------------------------------------------------------//

// [[Rcpp::export]]
arma::mat KR_DSSlow(arma::mat yMat, arma::colvec hVec, arma::colvec drvVec,
                    SEXP kernFcnPtrX, SEXP kernFcnPtrT)
{
  int nRows = yMat.n_rows;
  int nCols = yMat.n_cols;
  arma::mat matTemp(nRows, nCols);
  arma::mat yMatOut(nCols, nRows);
  
  for(int i{ 0 }; i < nRows; ++i)
  {
    matTemp.row(i) = KRSlow_single(yMat.row(i), hVec(1), kernFcnPtrT)/
                      pow(hVec(1), drvVec(1));
  }

  arma::mat matTempT{ matTemp.t() };

  for(int i{ 0 }; i < nCols; ++i)
  {
    yMatOut.row(i) = KRSlow_single(matTempT.row(i), hVec(0), kernFcnPtrX)/
                           pow(hVec(0), drvVec(0));
  }

  return yMatOut.t();
}