// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include "DCSmooth_types.h"

using namespace Rcpp;

// Müller-Wang kernels
// Form is: Type_k_mu_nu

arma::vec kernFkt_MW200(const arma::vec& u, double q = 1)
{
  arma::vec qVec{ arma::vec(u.n_rows).fill(q) };
  arma::vec q2{ pow(qVec, 2) };
  double n0{ 2.0/pow(q + 1, 3) };

  return n0 * (2*(q2 - qVec + 1) + 3*u%(1 - qVec));
}

arma::vec kernFkt_MW210(const arma::vec& u, double q = 1)
{
  arma::vec qVec{ arma::vec(u.n_rows).fill(q) };
  arma::vec wq{ (1 + u) };
  arma::vec q2{ pow(qVec, 2) };
  double n0{ 6.0/pow(q + 1, 4) };

  return n0 * ((3*q2 - 2*qVec + 1) + 2*u % (1 - 2*qVec)) % wq;
}

arma::vec kernFkt_MW220(const arma::vec& u, double q = 1)
{
  arma::vec qVec{ arma::vec(u.n_rows).fill(q) };
  arma::vec wq{ (1 + u) % (1 + u) % (qVec - u) };
  arma::vec q2{ pow(qVec, 2) };
  double n0{ 60.0/pow(q + 1, 6) };

  return n0 * ((2*q2 - 2*qVec + 1) - (3*qVec - 2) % u) % wq;
}

// [[Rcpp::export]]
arma::vec kernFkt_MW320(const arma::vec& u, double q = 1)
{
  arma::vec qVec{ arma::vec(u.n_rows).fill(q) };
  arma::vec wq{ (1 + u) % (1 + u) % (qVec - u) };
  arma::vec q4{ pow(qVec, 4) };
  arma::vec q3{ pow(qVec, 3) };
  arma::vec q2{ pow(qVec, 2) };
  double n0{ 60.0/pow(q + 1, 8) };
  arma::vec mid{ (10*q4 - 30*q3 + 39*q2 - 16*qVec + 3)
                  - 7*(5*qVec - 2) % (qVec - 1) % (qVec - 1) % u
                  + 14*(2*q2 - 4*qVec + 1) % u % u };

  return (n0 * mid % wq);
}

// [[Rcpp::export]]
arma::vec kernFkt_MW420(const arma::vec& u, double q = 1)
{
  arma::vec qVec{ arma::vec(u.n_rows).fill(q) };
  arma::vec wq{ (1 + u) % (1 + u) % (qVec - u) };
  arma::vec q6{ pow(qVec, 6) };
  arma::vec q5{ pow(qVec, 5) };
  arma::vec q4{ pow(qVec, 4) };
  arma::vec q3{ pow(qVec, 3) };
  arma::vec q2{ pow(qVec, 2) };
  double n0{ 420.0/pow(q + 1, 10) };
  arma::vec mid{ (5*q6 - 30*q5%u - 30*q5 + 54*q4%u%u + 152*q4%u + 77*q4
                    - 30*q3%u%u%u - 240*q3%u%u - 280*q3%u - 84*q3
                    + 120*q2%u%u%u + 300*q2%u%u + 216*q2%u + 45*q2
                    - 90*qVec%u%u%u - 144*qVec%u%u - 70*qVec%u - 10*qVec
                    + 12*u%u%u + 18*u%u + 8*u + 1) };
    return (n0 * mid % wq);
}

// [[Rcpp::export]]
arma::vec kernFkt_MW422(const arma::vec& u, double q = 1)
{
  arma::vec qVec{ arma::vec(u.n_rows).fill(q) };
  arma::vec wq{ (1 + u) % (1 + u) % (qVec - u) };
  arma::vec q4{ pow(qVec, 4) };
  arma::vec q3{ pow(qVec, 3) };
  arma::vec q2{ pow(qVec, 2) };
  double n0{ 5040.0/pow(q + 1, 10) };
  arma::vec mid{ (9*q4 - 61*q3%u - 40*q3 + 119*q2%u%u + 160*q2%u + 50*q2
                 - 70*qVec%u%u%u - 182*qVec%u%u - 127*qVec%u - 24*qVec
                 + 56*u%u%u + 77*u%u + 30*u + 3) };
  return (n0 * mid % wq);
}

// [[Rcpp::export]]
XPtr<funcPtr> kernelFcn_assign(std::string fstr) {
  if (fstr == "MW200")
    return(XPtr<funcPtr>(new funcPtr(&kernFkt_MW200)));
  else if (fstr == "MW210")
    return(XPtr<funcPtr>(new funcPtr(&kernFkt_MW210)));
  else if (fstr == "MW220")
    return(XPtr<funcPtr>(new funcPtr(&kernFkt_MW220)));
  else if (fstr == "MW320")
    return(XPtr<funcPtr>(new funcPtr(&kernFkt_MW320)));
  else if (fstr == "MW420")
    return(XPtr<funcPtr>(new funcPtr(&kernFkt_MW420)));
  else if (fstr == "MW422")
    return(XPtr<funcPtr>(new funcPtr(&kernFkt_MW422)));
  else
    return XPtr<funcPtr>(R_NilValue); // runtime error as NULL no XPtr
}

// [[Rcpp::export]]
arma::vec kernelFcn_use(const arma::vec x, double q, SEXP xpsexp) {
  XPtr<funcPtr> xpfun(xpsexp);
  funcPtr fun = *xpfun;
  arma::vec y = fun(x, q);
  return y;
}