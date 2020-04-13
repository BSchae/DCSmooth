// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include "DCSmooth_types.h"
#include <RcppArmadillo.h>
#include <RcppEigen.h>
#include <Rcpp.h>

using namespace Rcpp;

// DoubleSmooth
arma::mat DoubleSmooth(arma::mat yMat, arma::colvec hVec, arma::colvec polyOrderS, arma::colvec drv);
RcppExport SEXP _DCSmooth_DoubleSmooth(SEXP yMatSEXP, SEXP hVecSEXP, SEXP polyOrderSSEXP, SEXP drvSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type yMat(yMatSEXP);
    Rcpp::traits::input_parameter< arma::colvec >::type hVec(hVecSEXP);
    Rcpp::traits::input_parameter< arma::colvec >::type polyOrderS(polyOrderSSEXP);
    Rcpp::traits::input_parameter< arma::colvec >::type drv(drvSEXP);
    rcpp_result_gen = Rcpp::wrap(DoubleSmooth(yMat, hVec, polyOrderS, drv));
    return rcpp_result_gen;
END_RCPP
}
// factorialFunction
int factorialFunction(int value);
RcppExport SEXP _DCSmooth_factorialFunction(SEXP valueSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type value(valueSEXP);
    rcpp_result_gen = Rcpp::wrap(factorialFunction(value));
    return rcpp_result_gen;
END_RCPP
}
// LPSmooth_matrix
arma::mat LPSmooth_matrix(const arma::mat yMat, const double h, const int polyOrder, const int drv);
RcppExport SEXP _DCSmooth_LPSmooth_matrix(SEXP yMatSEXP, SEXP hSEXP, SEXP polyOrderSEXP, SEXP drvSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat >::type yMat(yMatSEXP);
    Rcpp::traits::input_parameter< const double >::type h(hSEXP);
    Rcpp::traits::input_parameter< const int >::type polyOrder(polyOrderSEXP);
    Rcpp::traits::input_parameter< const int >::type drv(drvSEXP);
    rcpp_result_gen = Rcpp::wrap(LPSmooth_matrix(yMat, h, polyOrder, drv));
    return rcpp_result_gen;
END_RCPP
}
// FastDoubleSmooth
arma::mat FastDoubleSmooth(arma::mat yMat, arma::colvec hVec, arma::icolvec polyOrderVec, arma::icolvec drvVec);
RcppExport SEXP _DCSmooth_FastDoubleSmooth(SEXP yMatSEXP, SEXP hVecSEXP, SEXP polyOrderVecSEXP, SEXP drvVecSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type yMat(yMatSEXP);
    Rcpp::traits::input_parameter< arma::colvec >::type hVec(hVecSEXP);
    Rcpp::traits::input_parameter< arma::icolvec >::type polyOrderVec(polyOrderVecSEXP);
    Rcpp::traits::input_parameter< arma::icolvec >::type drvVec(drvVecSEXP);
    rcpp_result_gen = Rcpp::wrap(FastDoubleSmooth(yMat, hVec, polyOrderVec, drvVec));
    return rcpp_result_gen;
END_RCPP
}
// kernelFkt_assign
XPtr<funcPtr> kernelFkt_assign(std::string fstr);
RcppExport SEXP _DCSmooth_kernelFkt_assign(SEXP fstrSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string >::type fstr(fstrSEXP);
    rcpp_result_gen = Rcpp::wrap(kernelFkt_assign(fstr));
    return rcpp_result_gen;
END_RCPP
}
// kernelFkt_use
arma::vec kernelFkt_use(const arma::vec x, double q, SEXP xpsexp);
RcppExport SEXP _DCSmooth_kernelFkt_use(SEXP xSEXP, SEXP qSEXP, SEXP xpsexpSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec >::type x(xSEXP);
    Rcpp::traits::input_parameter< double >::type q(qSEXP);
    Rcpp::traits::input_parameter< SEXP >::type xpsexp(xpsexpSEXP);
    rcpp_result_gen = Rcpp::wrap(kernelFkt_use(x, q, xpsexp));
    return rcpp_result_gen;
END_RCPP
}
// weightMatrix
arma::mat weightMatrix(arma::colvec weights, arma::mat matrix);
RcppExport SEXP _DCSmooth_weightMatrix(SEXP weightsSEXP, SEXP matrixSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::colvec >::type weights(weightsSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type matrix(matrixSEXP);
    rcpp_result_gen = Rcpp::wrap(weightMatrix(weights, matrix));
    return rcpp_result_gen;
END_RCPP
}
// LPSmooth_grid
arma::colvec LPSmooth_grid(const arma::colvec y, const double h, const int polyOrder, const int drv);
RcppExport SEXP _DCSmooth_LPSmooth_grid(SEXP ySEXP, SEXP hSEXP, SEXP polyOrderSEXP, SEXP drvSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::colvec >::type y(ySEXP);
    Rcpp::traits::input_parameter< const double >::type h(hSEXP);
    Rcpp::traits::input_parameter< const int >::type polyOrder(polyOrderSEXP);
    Rcpp::traits::input_parameter< const int >::type drv(drvSEXP);
    rcpp_result_gen = Rcpp::wrap(LPSmooth_grid(y, h, polyOrder, drv));
    return rcpp_result_gen;
END_RCPP
}
// LPSmooth_nongrid
arma::colvec LPSmooth_nongrid(const arma::colvec y, const arma::colvec x, const double h, const int polyOrder);
RcppExport SEXP _DCSmooth_LPSmooth_nongrid(SEXP ySEXP, SEXP xSEXP, SEXP hSEXP, SEXP polyOrderSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::colvec >::type y(ySEXP);
    Rcpp::traits::input_parameter< const arma::colvec >::type x(xSEXP);
    Rcpp::traits::input_parameter< const double >::type h(hSEXP);
    Rcpp::traits::input_parameter< const int >::type polyOrder(polyOrderSEXP);
    rcpp_result_gen = Rcpp::wrap(LPSmooth_nongrid(y, x, h, polyOrder));
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
arma::colvec test3(arma::colvec a, arma::colvec b);
RcppExport SEXP _DCSmooth_test3(SEXP aSEXP, SEXP bSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::colvec >::type a(aSEXP);
    Rcpp::traits::input_parameter< arma::colvec >::type b(bSEXP);
    rcpp_result_gen = Rcpp::wrap(test3(a, b));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_DCSmooth_DoubleSmooth", (DL_FUNC) &_DCSmooth_DoubleSmooth, 4},
    {"_DCSmooth_factorialFunction", (DL_FUNC) &_DCSmooth_factorialFunction, 1},
    {"_DCSmooth_LPSmooth_matrix", (DL_FUNC) &_DCSmooth_LPSmooth_matrix, 4},
    {"_DCSmooth_FastDoubleSmooth", (DL_FUNC) &_DCSmooth_FastDoubleSmooth, 4},
    {"_DCSmooth_kernelFkt_assign", (DL_FUNC) &_DCSmooth_kernelFkt_assign, 1},
    {"_DCSmooth_kernelFkt_use", (DL_FUNC) &_DCSmooth_kernelFkt_use, 3},
    {"_DCSmooth_weightMatrix", (DL_FUNC) &_DCSmooth_weightMatrix, 2},
    {"_DCSmooth_LPSmooth_grid", (DL_FUNC) &_DCSmooth_LPSmooth_grid, 4},
    {"_DCSmooth_LPSmooth_nongrid", (DL_FUNC) &_DCSmooth_LPSmooth_nongrid, 4},
    {"_DCSmooth_test", (DL_FUNC) &_DCSmooth_test, 2},
    {"_DCSmooth_test3", (DL_FUNC) &_DCSmooth_test3, 2},
    {NULL, NULL, 0}
};

RcppExport void R_init_DCSmooth(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
