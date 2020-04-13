// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include "KernelFunctions.h"

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

// [[Rcpp::export]]
XPtr<funcPtr> putFunPtrInXPtr(std::string fstr) {
  if (fstr == "fun1")
    return(XPtr<funcPtr>(new funcPtr(&fun1_cpp)));
  else if (fstr == "fun2")
    return(XPtr<funcPtr>(new funcPtr(&fun2_cpp)));
  else
    return XPtr<funcPtr>(R_NilValue); // runtime error as NULL no XPtr
}