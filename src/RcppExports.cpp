// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// fastLm
List fastLm(NumericVector y, NumericMatrix x, NumericMatrix w);
RcppExport SEXP _DCSmooth_fastLm(SEXP ySEXP, SEXP xSEXP, SEXP wSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type y(ySEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type x(xSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type w(wSEXP);
    rcpp_result_gen = Rcpp::wrap(fastLm(y, x, w));
    return rcpp_result_gen;
END_RCPP
}
// test
List test(NumericVector x, double q);
RcppExport SEXP _DCSmooth_test(SEXP xSEXP, SEXP qSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< double >::type q(qSEXP);
    rcpp_result_gen = Rcpp::wrap(test(x, q));
    return rcpp_result_gen;
END_RCPP
}
// test3
int test3(int x);
RcppExport SEXP _DCSmooth_test3(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(test3(x));
    return rcpp_result_gen;
END_RCPP
}
// rcpparma_hello_world
arma::mat rcpparma_hello_world();
RcppExport SEXP _DCSmooth_rcpparma_hello_world() {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    rcpp_result_gen = Rcpp::wrap(rcpparma_hello_world());
    return rcpp_result_gen;
END_RCPP
}
// rcpparma_outerproduct
arma::mat rcpparma_outerproduct(const arma::colvec& x);
RcppExport SEXP _DCSmooth_rcpparma_outerproduct(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::colvec& >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(rcpparma_outerproduct(x));
    return rcpp_result_gen;
END_RCPP
}
// rcpparma_innerproduct
double rcpparma_innerproduct(const arma::colvec& x);
RcppExport SEXP _DCSmooth_rcpparma_innerproduct(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::colvec& >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(rcpparma_innerproduct(x));
    return rcpp_result_gen;
END_RCPP
}
// rcpparma_bothproducts
Rcpp::List rcpparma_bothproducts(const arma::colvec& x);
RcppExport SEXP _DCSmooth_rcpparma_bothproducts(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::colvec& >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(rcpparma_bothproducts(x));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_DCSmooth_fastLm", (DL_FUNC) &_DCSmooth_fastLm, 3},
    {"_DCSmooth_test", (DL_FUNC) &_DCSmooth_test, 2},
    {"_DCSmooth_test3", (DL_FUNC) &_DCSmooth_test3, 1},
    {"_DCSmooth_rcpparma_hello_world", (DL_FUNC) &_DCSmooth_rcpparma_hello_world, 0},
    {"_DCSmooth_rcpparma_outerproduct", (DL_FUNC) &_DCSmooth_rcpparma_outerproduct, 1},
    {"_DCSmooth_rcpparma_innerproduct", (DL_FUNC) &_DCSmooth_rcpparma_innerproduct, 1},
    {"_DCSmooth_rcpparma_bothproducts", (DL_FUNC) &_DCSmooth_rcpparma_bothproducts, 1},
    {NULL, NULL, 0}
};

RcppExport void R_init_DCSmooth(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
