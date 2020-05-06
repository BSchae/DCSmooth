// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <cstdlib>
#include "DCSmooth_types.h"

using namespace Rcpp;

// [[Rcpp::export]]
int cppSample(arma::vec x)
{
  int randomIndex{ std::rand() % x.n_rows };
  int out = static_cast<int>(x(randomIndex));
  return out;
}

// [[Rcpp::export]]
arma::mat thinnedMat(arma::mat yMat, int seed) //, int sRow, int sCol)
{
  std::srand(seed);

  int sRow{ yMat.n_rows/1000 + 1};
  int sCol{ yMat.n_cols/1000 + 1};

  int nRow{ yMat.n_rows/sRow };
  int nCol{ yMat.n_cols/sCol };

  arma::vec rowRange = arma::regspace(0, sRow - 1);
  arma::vec colRange = arma::regspace(0, sCol - 1);

  int rowMark{ 0 };
  int colMark{ 0 };

  arma::mat yMatOut(nRow, nCol);


  for (int rowIndex{ 0 }; rowIndex < nRow; ++rowIndex)
  {
    for (int colIndex{ 0 }; colIndex < nCol; ++colIndex)
    {
      rowMark = cppSample(rowRange);
      colMark = cppSample(colRange);
      yMatOut(rowIndex, colIndex) = yMat(rowMark, colMark);
      colRange += sCol;
    }
    rowRange += sRow;
    colRange = arma::regspace(0, sCol - 1);
  }

  return yMatOut;
}