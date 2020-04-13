// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include "DCSmooth_types.h"

using namespace arma; 
using namespace Rcpp;


vec fun1_cpp(const vec& x) {	// a first function 
  vec y = x + x;
  return (y);
}

vec fun2_cpp(const vec& x) {	// and a second function
  vec y = 10*x;
  return (y);
}


// Müller-Wang kernels
// Form is: Type_k_mu_nu

arma::vec kernFkt_MW200(const arma::vec& x, double q)
{
  return (6*x - 4*q - 6*q*x + 4*pow(q, 2) + 4)/pow(q + 1, 3);
}

arma::vec kernFkt_MW210(const arma::vec& x, double q)
{
  return (6*(x + 1) % (2*x - 2*q - 4*q*x + 3*pow(q, 2) + 1))/pow(q + 1, 4);
}

arma::vec kernFkt_MW220(const arma::vec& x, double q)
{
  return (60*(q - x)*pow(x + 1, 2)*(2*x - 2*q - 3*q*x + 2*pow(q, 2) + 1))/pow(q + 1, 6);
}

// [[Rcpp::export]]
XPtr<funcPtr> kernelFkt_assign(std::string fstr) {
  if (fstr == "MW200")
    return(XPtr<funcPtr>(new funcPtr(&kernFkt_MW200)));
  else if (fstr == "MW210")
    return(XPtr<funcPtr>(new funcPtr(&kernFkt_MW210)));
  else
    return XPtr<funcPtr>(R_NilValue); // runtime error as NULL no XPtr
}

// [[Rcpp::export]]
arma::vec kernelFkt_use(const arma::vec x, double q, SEXP xpsexp) {
  XPtr<funcPtr> xpfun(xpsexp);
  funcPtr fun = *xpfun;
  vec y = fun(x, q);
  return (y);
}