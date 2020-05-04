// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include "DCSmooth_types.h"

using namespace Rcpp;

// Müller-Wang kernels
// Form is: Type_k_mu_nu

arma::vec kernFkt_MW200(arma::vec& u, double q = 1)
{
  arma::vec qVec{ arma::vec(u.n_rows).fill(q) };
  arma::vec q2{ pow(qVec, 2) };
  double n0{ 2.0/pow(q + 1, 3) };

  return n0 * (2*(q2 - qVec + 1) + 3*u%(1 - qVec));
}

arma::vec kernFkt_MW210(arma::vec& u, double q = 1)
{
  arma::vec qVec{ arma::vec(u.n_rows).fill(q) };
  arma::vec wq{ (1 + u) };
  arma::vec q2{ pow(qVec, 2) };
  double n0{ 6.0/pow(q + 1, 4) };

  return n0 * ((3*q2 - 2*qVec + 1) + 2*u % (1 - 2*qVec)) % wq;
}

// [[Rcpp::export]]
arma::vec kernFkt_MW220(arma::vec& uVec, double q)
{
  int nBound{ uVec.size() };
  
  arma::vec uOut(nBound);
  
  double q2{ q * q };
  double q2_1{ pow(1 + q, 6) };
  double out{ };
  double u{ };
  
  for (int i{ 0 }; i < nBound; ++i)
  {
    u   = uVec(i);
    out = (60.0/q2_1) * (u*(2 - 3*q) + (1 - 2*q + 2*q2))
                      * (1 + u)*(1 + u)*(q - u);
    uOut(i) = out;
  }
  return uOut;
}

// [[Rcpp::export]]
arma::vec kernFkt_MW220_old(arma::vec& u, double q = 1)
{
  arma::vec qVec{ arma::vec(u.n_rows).fill(q) };
  arma::vec wq{ (1 + u) % (1 + u) % (qVec - u) };
  arma::vec q2{ pow(qVec, 2) };
  double n0{ 60.0/pow(q + 1, 6) };

  return n0 * ((2*q2 - 2*qVec + 1) - (3*qVec - 2) % u) % wq;
}

// [[Rcpp::export]]
arma::vec kernFkt_MW320(arma::vec& u, double q = 1)
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
arma::vec kernFkt_MW420(arma::vec& uVec, double q)
{
  int nBound{ uVec.size() };
  
  arma::vec uOut(nBound);
  
  double q2{ q * q };
  double q3{ q2 * q };
  double q4{ q3 * q };
  double q5{ q4 * q };
  double q6{ q5 * q };
  double q10_1{ pow(1 + q, 10) };
  double out{ };
  double u{ };
  
  for (int i{ 0 }; i < nBound; ++i)
  {
    u   = uVec(i);
    out = (420.0/q10_1) * ( u*u*u*(12 - 90*q + 120*q2 - 30*q3) 
                  + u*u*(18 - 144*q + 300*q2 - 240*q3 + 54*q4)
                  + u*(8 - 70*q + 216*q2 - 280*q3 + 152*q4 - 30*q5)
                  + (1 - 10*q + 45*q2 - 84*q3 + 77*q4 - 30*q5 + 5*q6))
                * (1 + u)*(1 + u)*(q - u);
    uOut(i) = out;
  }
  return uOut;
}

// [[Rcpp::export]]
arma::vec kernFkt_MW422(arma::vec& uVec, double q)
{
  int nBound{ uVec.size() };
  
  arma::vec uOut(nBound);
  
  double q2{ q * q };
  double q3{ q2 * q };
  double q4{ q3 * q };
  double q10_1{ pow(1 + q, 10) };
  double out{ };
  double u{ };
  
  for (int i{ 0 }; i < nBound; ++i)
  {
    u   = uVec(i);
    out = (5040.0/q10_1) * ( u*u*u*(56 - 70*q) + u*u*(77 - 182*q + 119*q2)
                         + u*(30 - 127*q + 160*q2 - 61*q3)
                         + (3 - 24*q + 50*q2 - 40*q3 + 9*q4))
                       * (1 + u)*(1 + u)*(q - u);
    uOut(i) = out;
  }
  return uOut;
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
arma::vec kernelFcn_use(arma::vec x, double q, SEXP xpsexp) {
  XPtr<funcPtr> xpfun(xpsexp);
  funcPtr fun = *xpfun;
  arma::vec y = fun(x, q);
  return y;
}





/// Test

// [[Rcpp::export]]
arma::vec MWTest1(arma::vec uVec, double q)
{
  int nBound{ uVec.size() };
  
  arma::vec uOut(nBound);
  
  double q2{ q * q };
  double q3{ q2 * q };
  double q4{ q3 * q };
  double q5{ q4 * q };
  double q6{ q5 * q };
  double q10_1{ pow(1 + q, 10) };
  double out{ };
  double u{ };
  
  for (int i{ 0 }; i < nBound; ++i)
  {
    u   = uVec(i);
    out = (420.0/q10_1) * ( u*u*u*(12 - 90*q + 120*q2 - 30*q3) 
      + u*u*(18 - 144*q + 300*q2 - 240*q3 + 54*q4)
      + u*(8 - 70*q + 216*q2 - 280*q3 + 152*q4 - 30*q5)
      + (1 - 10*q + 45*q2 - 84*q3 + 77*q4 - 30*q5 + 5*q6))
      * (1 + u)*(1 + u)*(q - u);
      uOut(i) = out;
  }
  return uOut;
}

// [[Rcpp::export]]
double MWTest2(double u, double q)
{
  double q2{ q * q };
  double q3{ q2 * q };
  double q4{ q3 * q };
  double q5{ q4 * q };
  double q6{ q5 * q };
  double q10_1{ pow(1 + q, 10) };
  
  double out = (420.0/q10_1) * ( u*u*u*(12 - 90*q + 120*q2 - 30*q3) 
    + u*u*(18 - 144*q + 300*q2 - 240*q3 + 54*q4)
    + u*(8 - 70*q + 216*q2 - 280*q3 + 152*q4 - 30*q5)
    + (1 - 10*q + 45*q2 - 84*q3 + 77*q4 - 30*q5 + 5*q6))
    * (1 + u)*(1 + u)*(q - u);
  
  return out;
}