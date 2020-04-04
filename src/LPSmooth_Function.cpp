// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <cmath>

using namespace Rcpp;

//---------------------------------------------------------------------------//

// declare xMatrix here
arma::mat xMatrix(arma::colvec xVector, int order);

//---------------------------------------------------------------------------//

arma::mat weightMatrix(arma::colvec weights, arma::mat matrix)
{
  arma::mat matrixOut{ arma::mat(matrix.n_rows, matrix.n_cols) };
  for (int j{ 0 }; j < matrix.n_cols; ++j)
  {
    matrixOut.col(j) = weights % matrix.col(j);
  }
  
  return matrixOut;
}

//---------------------------------------------------------------------------//

// [[Rcpp::export]]
arma::colvec LPSmooth(const arma::colvec y,
                      const double h, const int polyOrder)
{

  arma::colvec yOut{ y*0 };               // vector for results (??? edit me please)
  int n{ y.n_rows };                      // number of observations
  int bndw{ static_cast<int>(h * n) };    // calculate absolute bandwidth, decimals will be dumped
  int windowWidth{ 2*bndw + 1 };          // width of estimation window
  
  
// smoothing over interior values
  arma::colvec  uInterior{ arma::linspace(-bndw, bndw, windowWidth)/(h * n) };          //
  arma::colvec  weightsInterior{ (0.75*pow(1 - pow(uInterior, 2), 2)) };
  arma::colvec  yInterior{ arma::zeros(windowWidth) };                                  // empty vector for use in loop
  arma::mat     xMatInterior{ xMatrix(arma::linspace(-bndw, bndw, windowWidth), polyOrder) }; //
  xMatInterior = weightMatrix(weightsInterior, xMatInterior);

  for (int index{ bndw }; index < (n - bndw); ++index)
  {
    yInterior = weightsInterior % y.subvec(index - bndw, index + bndw);
    arma::mat coefMatrix{ arma::solve(xMatInterior, yInterior) };
    yOut(index) = coefMatrix(0,0);
  }

  
// smoothing over left boundary
  arma::colvec  xSubLeft{ x.subvec(0, 2 * bndw) }; // constant window width at boundaries
  arma::colvec  ySubLeft{ y.subvec(0, 2 * bndw) }; // constant window width at boundaries

  // initialise empty variables for data inside the loop
  arma::colvec  uLeft(xSubLeft.n_rows);                        // vector for holding normalized x values around x0
  arma::mat     xWeightedMat(xSubLeft.n_rows, polyOrder + 1);  // holds weighted x-matrix for lm
  arma::mat     weights(xSubLeft.n_rows, xSubLeft.n_rows);     // holds computed weights
  arma::colvec  yWeighted(ySubLeft.n_rows);                    // holds weigted y-values

  for (int index{ 0 }; index < bndw; ++index)
  {
    double x0{ x(index) };
    uLeft   = (xSubLeft - x0)/h;
    weights = arma::diagmat(pow(1 - pow(uLeft, 2), 2));
    // compute weighted values for x and Y
    xWeightedMat = xMatrix(xSubLeft, polyOrder);
    yWeighted    = weights * ySubLeft;

    arma::mat coefMatrix{ arma::solve(xWeightedMat, yWeighted) };
    yOut(index) = coefMatrix(0,0);
  }

  return yOut;
}