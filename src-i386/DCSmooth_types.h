#include <RcppArmadillo.h>
#include <RcppEigen.h>

typedef arma::vec (*funcPtr)(arma::vec&, double);