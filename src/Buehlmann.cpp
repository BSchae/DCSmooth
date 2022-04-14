// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
// #include <cmath>
#include "DCSmooth.h"
#include "DCSmooth_types.h"

using namespace Rcpp;

arma::mat acfMatrix_cpp(const arma::mat& y_Mat);
arma::mat acfMatrix_quarter(const arma::mat& y0_Mat);
arma::mat bartlett_mat(int Lx, int Lt, arma::vec drv1, arma::vec drv2,
                       int power);
double specInt00_cpp(const arma::mat& acfMat);
double specIntDrv_cpp(const arma::mat& acfMat, arma::vec drv1, arma::vec drv2,
                      arma::vec hLag);
arma::cx_mat fourier_mat(int Lx, int Lt, arma::vec omega);
double specDensEst_cpp(const arma::mat& acfMat, arma::vec drv, arma::vec hLag,
                       arma::vec omega);
arma::vec localBndw_cpp(const arma::mat& acfMat, arma::vec hLag,
                        arma::vec omega);
arma::vec globalBndw_cpp(const arma::mat& acfMat, arma::vec hLag);

//-------------------------------Main Function--------------------------------//

//[[Rcpp::export]]
arma::vec specDens_cpp(const arma::mat& y_mat, arma::vec omega)
{
  int nX{ static_cast<int>(y_mat.n_rows) };
  int nT{ static_cast<int>(y_mat.n_cols) };

  arma::mat acfMat{ acfMatrix_cpp(y_mat) };

  // initial bandwidth values
  arma::vec hVec(2);
  arma::vec hVecInfl(2);
  arma::vec hVecTemp(2);
  arma::vec inflFct(2);
  hVec(0) = trunc(nX/2); hVec(1) = trunc(nT/2);
  inflFct(0) = 1/pow(nX, 2.0/21.0); inflFct(1) = 1/pow(nT, 2.0/21.0);

  // global step
  for (int g{ 0 }; g < 5; g++)
  {
    hVecTemp = hVec;
    hVecInfl = trunc(hVec % inflFct) + 1;
    hVec = globalBndw_cpp(acfMat, hVecInfl);
    if (arma::approx_equal(hVec, hVecTemp, "absdiff", 0.001))
    {
      break;
    }
  }

  // local step
  hVecInfl = trunc(hVec % inflFct) + 1;

  arma::vec hOpt{ localBndw_cpp(acfMat, hVecInfl, omega) };

  arma::vec drv0 = { 0, 0 };

  double specDensOut{ specDensEst_cpp(acfMat, drv0, hOpt, omega) };

  arma::vec returnVec(3);
  returnVec(0) = specDensOut * pow(2 * M_PI, 2);
  returnVec.subvec(1, 2) = hOpt;

  return returnVec;
}

//---------------------------Calculation of ACF-------------------------------//

// estimation of complete acf matrix
//[[Rcpp::export]]
arma::mat acfMatrix_cpp(const arma::mat& y0_Mat)
{
  int nRow{ static_cast<int>(y0_Mat.n_rows) };
  int nCol{ static_cast<int>(y0_Mat.n_cols) };

  // matrix should already be centered
  // arma::mat y0_Mat{ y_Mat - arma::accu(y_Mat)/(nRow * nCol) };
  // top right submatrix
  arma::mat y1_Mat{ reverse(y0_Mat, 0) }; // a bit slow here

  arma::mat acf0_Mat{ acfMatrix_quarter(y0_Mat) };
  arma::mat acf1_Mat{ acfMatrix_quarter(y1_Mat) };

  arma::mat acfMat_out(2*nRow - 1, 2*nCol - 1, arma::fill::zeros);
  acfMat_out.submat(0, 0, nRow - 1, nCol - 1) =
    reverse(reverse(acf0_Mat, 1), 0);
  acfMat_out.submat(nRow - 1, nCol - 1, 2*nRow - 2, 2*nCol - 2) = acf0_Mat;
  acfMat_out.submat(nRow - 1, 0, 2*nRow - 2, nCol - 1) =
    reverse(acf1_Mat, 1);
  acfMat_out.submat(0, nCol - 1, nRow - 1, 2*nCol - 2) =
    reverse(acf1_Mat, 0);

  return acfMat_out;
}

// estimation of acf matrix for s = 0,...,nX; u = 0,...,nT

//[[Rcpp::export]]
arma::mat acfMatrix_quarter(const arma::mat& y_Mat)
{
  int nRow{ static_cast<int>(y_Mat.n_rows) };
  int nCol{ static_cast<int>(y_Mat.n_cols) };
  arma::mat y0;
  arma::mat y1;
  arma::mat R_out;

  R_out.zeros(nRow, nCol);

  for (arma::uword j{ 0 }; j < nCol; j++)
  {
    for (arma::uword i{ 0 }; i < nRow; i++)
    {
      y0 = y_Mat.submat(0, 0, nRow - i - 1, nCol - j - 1);
      y1 = y_Mat.submat(i, j, nRow - 1, nCol - 1);
      R_out.at(i, j) = arma::accu(y0 % y1);
    }
  }

  return R_out/(nRow * nCol);
}

//---------------------Estimation of Global Bandwidth-------------------------//

//[[Rcpp::export]]
arma::vec globalBndw_cpp(const arma::mat& acfMat, arma::vec hLag)
{
  int nX{ static_cast<int>(acfMat.n_rows)/2 + 1 };
  int nT{ static_cast<int>(acfMat.n_cols)/2 + 1 };
  double mu_2w { 4.0/9.0 };

  double F00;
  double F10;
  double F01;
  double F11;

  double Fx;
  double Ft;

  int Lx;
  int Lt;

  arma::vec drv10{ 1, 0 };
  arma::vec drv01{ 0, 1 };

  if (hLag(0) > 1 && hLag(1) > 1)
  {
    F00 = specInt00_cpp(acfMat);
    F10 = specIntDrv_cpp(acfMat, drv10, drv10, hLag);
    F01 = specIntDrv_cpp(acfMat, drv01, drv01, hLag);
    F11 = specIntDrv_cpp(acfMat, drv10, drv01, hLag);

    Fx = F10 * (sqrt(F10/F01) + F11/F01) / F00;
    Ft = F01 * (sqrt(F01/F10) + F11/F10) / F00;

    Lx = static_cast<int>(pow(2 * Fx * nX*nT/mu_2w, 0.25));
    Lt = static_cast<int>(pow(2 * Ft * nX*nT/mu_2w, 0.25));
  }
  else if (hLag(0) == 1 && hLag(1) > 1)
  {
    F00 = specInt00_cpp(acfMat);
    F10 = specIntDrv_cpp(acfMat, drv10, drv10, hLag);

    Fx = F10/F00;

    Lx = static_cast<int>(pow(2 * Fx * nX*nT/mu_2w, 1.0/3.0));
    Lt = 0;
  }
  else if (hLag(0) > 1 && hLag(1) == 1)
  {
    F00 = specInt00_cpp(acfMat);
    F01 = specIntDrv_cpp(acfMat, drv01, drv01, hLag);

    Ft =  F01/F00;

    Lx = 0;
    Lt = static_cast<int>(pow(2 * Ft * nX*nT/mu_2w, 1.0/3.0));
  }
  else
  {
    Lx = 0;
    Lt = 0;
  }

  arma::vec hOut = { static_cast<double>(std::min(Lx + 1, nX - 1)),
                     static_cast<double>(std::min(Lt + 1, nT - 1)) };

  return hOut;
}

//[[Rcpp::export]]
double specInt00_cpp(const arma::mat& acfMat)
{
  double out{ accu(acfMat % acfMat) };
  return pow(2 * M_PI, - 2) * out * 0.5;
}

//[[Rcpp::export]]
double specIntDrv_cpp(const arma::mat& acfMat, arma::vec drv1, arma::vec drv2,
                      arma::vec hLag)
{
  int origin_x{ static_cast<int>(acfMat.n_rows)/2 };  // as first element has index 0,
  int origin_t{ static_cast<int>(acfMat.n_cols)/2 };  // origin is at floor(n/2)
  int Lx = hLag(0);
  int Lt = hLag(1);

  arma::mat w_mat(2*(Lx - 1) + 1, 2*(Lt - 1) + 1);

  // compute matrix of bartlett weights without "fourier factor"
  arma::mat weights = bartlett_mat(Lx, Lt, drv1, drv2, 2);

  if (Lx > 1 && Lt > 1)
  {
    w_mat.submat(0, 0, Lx - 1, Lt - 1) =
      reverse(reverse(weights.submat(0, 0, Lx - 1, Lt - 1), 1), 0);
    w_mat.submat(Lx, 0, 2*(Lx - 1), Lt - 2) =
      reverse(weights.submat(1, 1, Lx - 1, Lt - 1), 1);
    w_mat.submat(0, Lt, Lx - 2, 2*(Lt - 1)) =
      reverse(weights.submat(1, 1, Lx - 1, Lt - 1), 0);
    w_mat.submat(Lx - 1, Lt - 1, 2*(Lx - 1), 2*(Lt - 1)) =
      weights.submat(0, 0, Lx - 1, Lt - 1);
    // are the "reverse" necessary, as w_mat is probably symmetric?
  } else if (Lx == 1 && Lt == 1) {
    w_mat.submat(0, 0, 0, 0) = 1;
  }

  arma::mat acf_submat{ acfMat.submat(origin_x - Lx + 1, origin_t - Lt + 1,
                                      origin_x + Lx - 1, origin_t + Lt - 1) };

  double F_out = accu(w_mat % acf_submat % acf_submat);

  return pow(2 * M_PI, - 2) * F_out;
}

//----------------------Estimation of Local Bandwidth-------------------------//

// [[Rcpp::export]]
arma::vec localBndw_cpp(const arma::mat& acfMat, arma::vec hLag, arma::vec omega)
{
  int nX{ static_cast<int>(acfMat.n_rows)/2 + 1 };
  int nT{ static_cast<int>(acfMat.n_cols)/2 + 1 };
  double mu_2w { 4.0/9.0 };

  double f00;
  double f10;
  double f01;
  double f11;

  double fx;
  double ft;

  int Lx;
  int Lt;

  arma::vec drv00{ 0, 0 };
  arma::vec drv10{ 1, 0 };
  arma::vec drv01{ 0, 1 };

  if (hLag(0) > 1 && hLag(1) > 1)
  {
    f00 = specDensEst_cpp(acfMat, drv00, hLag, omega);
    f10 = specDensEst_cpp(acfMat, drv10, hLag, omega);
    f01 = specDensEst_cpp(acfMat, drv01, hLag, omega);

    Lx = static_cast<int>(pow(4 * std::abs(pow(f10, 3)/f01) * 1/pow(f00, 2) *
                          nX*nT/mu_2w, 0.25));
    Lt = static_cast<int>(pow(4 * std::abs(pow(f01, 3)/f10) * 1/pow(f00, 2) *
                          nX*nT/mu_2w, 0.25));
  }
  else if (hLag(0) == 1 && hLag(1) > 1)
  {
    f00 = specDensEst_cpp(acfMat, drv00, hLag, omega);
    f10 = specDensEst_cpp(acfMat, drv10, hLag, omega);

    Lx = static_cast<int>(pow(2 * pow(f10/f00, 2) *
                          nX*nT/sqrt(mu_2w), 1.0/3.0));
    Lt = 0;
  }
  else if (hLag(0) > 1 && hLag(1) == 1)
  {
    f00 = specDensEst_cpp(acfMat, drv00, hLag, omega);
    f01 = specDensEst_cpp(acfMat, drv01, hLag, omega);

    Lx = 0;
    Lt = static_cast<int>(pow(2 * pow(f01/f00, 2) *
                          nX*nT/sqrt(mu_2w), 1.0/3.0)); // power 1/3 is correct?
  }
  else
  {
    Lx = 0;
    Lt = 0;
  }

  arma::vec hOut = { static_cast<double>(std::min(Lx + 1, nX - 1)),
                     static_cast<double>(std::min(Lt + 1, nT - 1)) };

  return hOut;
}

//[[Rcpp::export]]
double specDensEst_cpp(const arma::mat& acfMat, arma::vec drv, arma::vec hLag,
                       arma::vec omega)
{
  int origin_x{ static_cast<int>(acfMat.n_rows)/2 };  // as first element has index 0,
  int origin_t{ static_cast<int>(acfMat.n_cols)/2 };  // origin is at floor(n/2)
  int Lx{ static_cast<int>(hLag(0)) };
  int Lt{ static_cast<int>(hLag(1)) };
  arma::vec drv0 = { 0, 0 };

  arma::cx_mat w_mat(2*(Lx - 1) + 1, 2*(Lt - 1) + 1);

  // compute matrix of bartlett weights including "fourier factor"
  arma::cx_mat weights{ bartlett_mat(Lx, Lt, drv, drv0, 1) %
                fourier_mat(Lx, Lt, omega) };

  if (Lx > 1 && Lt > 1)
  {
    w_mat.submat(0, 0, Lx - 1, Lt - 1) =
      reverse(reverse(weights.submat(0, 0, Lx - 1, Lt - 1), 1), 0);
    w_mat.submat(Lx, 0, 2*(Lx - 1), Lt - 2) =
      reverse(weights.submat(1, 1, Lx - 1, Lt - 1), 1);
    w_mat.submat(0, Lt, Lx - 2, 2*(Lt - 1)) =
      reverse(weights.submat(1, 1, Lx - 1, Lt - 1), 0);
    w_mat.submat(Lx - 1, Lt - 1, 2*(Lx - 1), 2*(Lt - 1)) =
      weights.submat(0, 0, Lx - 1, Lt - 1);
    // are the "reverse" necessary, as w_mat is probably symmetric?
  } else if (Lx == 1 && Lt == 1) {
    w_mat.submat(0, 0, 0, 0) = 1;
  }

  arma::mat acf_submat{ acfMat.submat(origin_x - Lx + 1, origin_t - Lt + 1,
                                      origin_x + Lx - 1, origin_t + Lt - 1) };

  // calculate f. Note that f(-kX, - kT) = f(kX, kT), off diagonal sub matrices
  // are different etc.
  // note that w is (Lx + 1)x(Lt + 1) but the outer values are all zeroes,
  // thus summation over them is not needed.
  double fOut{ arma::accu(arma::real(w_mat) % acf_submat) };

  return pow(2 * M_PI, -2) * fOut;
}

//---------------------------Additional Functions-----------------------------//

// [[Rcpp::export]]
arma::cx_mat fourier_mat(int Lx, int Lt, arma::vec omega)
{
  arma::vec arg_x{ -arma::regspace(0, Lx) * omega(0) };
  arma::vec arg_t{ -arma::regspace(0, Lt) * omega(1) };

  arma::cx_mat mat_out{ arma::cx_vec(cos(arg_x), sin(arg_x)) *
          arma::cx_vec(cos(arg_t), sin(arg_t)).st() };

  return mat_out;
}

//[[Rcpp::export]]
arma::mat bartlett_mat(int Lx, int Lt, arma::vec drv1, arma::vec drv2,
                       int power)
{
  arma::mat lx_mat(Lx + 1, Lt + 1);
  arma::mat lt_mat(Lt + 1, Lx + 1); // needs to be transposed later
  arma::vec lx_col{ pow(arma::regspace(0, Lx), drv1(0) + drv2(0)) };
  arma::vec lt_row{ pow(arma::regspace(0, Lt), drv1(1) + drv2(1)) };

  for (int i{ 0 }; i <= Lt; i++)
  {
    lx_mat.col(i) = lx_col;
  }

  for (int j{ 0 }; j <= Lx; j++)
  {
    lt_mat.col(j) = lt_row; // cols are rows after transposing
  }

  arma::mat weights = pow((1 - arma::regspace(0, Lx)/std::max(Lx, 1)) *
    (1 - arma::regspace(0, Lt).t()/std::max(Lt, 1)), power) %
    lx_mat % lt_mat.t();

  return weights;
}