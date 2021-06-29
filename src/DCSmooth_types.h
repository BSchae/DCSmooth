#include <RcppArmadillo.h>
#include <RcppEigen.h>

typedef arma::vec (*funcPtr)(arma::vec&, double);
typedef arma::vec (*weightPtr)(arma::vec&, double, int);