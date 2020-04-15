// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include "DCSmooth_types.h"

using namespace Rcpp;

// Müller-Wang kernels
// Form is: Type_k_mu_nu

arma::vec kernFkt_MW200(const arma::vec& u, double q)
{
  arma::vec qVec{ arma::vec(u.n_rows).fill(q) };
  arma::vec q2{ pow(qVec, 2) };
  double n0{ 2.0/pow(q + 1, 3) };

  return n0 * (2*(q2 - qVec + 1) + 3*u%(1 - qVec));
}

arma::vec kernFkt_MW210(const arma::vec& u, double q)
{
  arma::vec qVec{ arma::vec(u.n_rows).fill(q) };
  arma::vec wq{ (1 + u) };
  arma::vec q2{ pow(qVec, 2) };
  double n0{ 6.0/pow(q + 1, 4) };

  return n0 * ((3*q2 - 2*qVec + 1) + 2*u % (1 - 2*qVec)) % wq;
}

arma::vec kernFkt_MW220(const arma::vec& u, double q)
{
  arma::vec qVec{ arma::vec(u.n_rows).fill(q) };
  arma::vec wq{ (1 + u) % (1 + u) % (qVec - u) };
  arma::vec q2{ pow(qVec, 2) };
  double n0{ 60.0/pow(q + 1, 6) };

  return n0 * ((2*q2 - 2*qVec + 1) - (3*qVec - 2) % u) % wq;
}

// [[Rcpp::export]]
XPtr<funcPtr> kernelFcn_assign(std::string fstr) {
  if (fstr == "MW200")
    return(XPtr<funcPtr>(new funcPtr(&kernFkt_MW200)));
  else if (fstr == "MW210")
    return(XPtr<funcPtr>(new funcPtr(&kernFkt_MW210)));
  else if (fstr == "MW220")
    return(XPtr<funcPtr>(new funcPtr(&kernFkt_MW220)));
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