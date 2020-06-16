// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include "DCSmooth_types.h"
#include <RcppArmadillo.h>
#include <RcppEigen.h>
#include <Rcpp.h>

using namespace Rcpp;

// KRSlow_single
arma::rowvec KRSlow_single(arma::rowvec yVec, double h, SEXP kernFcnPtr);
RcppExport SEXP _DCSmooth_KRSlow_single(SEXP yVecSEXP, SEXP hSEXP, SEXP kernFcnPtrSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::rowvec >::type yVec(yVecSEXP);
    Rcpp::traits::input_parameter< double >::type h(hSEXP);
    Rcpp::traits::input_parameter< SEXP >::type kernFcnPtr(kernFcnPtrSEXP);
    rcpp_result_gen = Rcpp::wrap(KRSlow_single(yVec, h, kernFcnPtr));
    return rcpp_result_gen;
END_RCPP
}
// KR_DSSlow
arma::mat KR_DSSlow(arma::mat yMat, arma::colvec hVec, arma::colvec drvVec, SEXP kernFcnPtrX, SEXP kernFcnPtrT);
RcppExport SEXP _DCSmooth_KR_DSSlow(SEXP yMatSEXP, SEXP hVecSEXP, SEXP drvVecSEXP, SEXP kernFcnPtrXSEXP, SEXP kernFcnPtrTSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type yMat(yMatSEXP);
    Rcpp::traits::input_parameter< arma::colvec >::type hVec(hVecSEXP);
    Rcpp::traits::input_parameter< arma::colvec >::type drvVec(drvVecSEXP);
    Rcpp::traits::input_parameter< SEXP >::type kernFcnPtrX(kernFcnPtrXSEXP);
    Rcpp::traits::input_parameter< SEXP >::type kernFcnPtrT(kernFcnPtrTSEXP);
    rcpp_result_gen = Rcpp::wrap(KR_DSSlow(yMat, hVec, drvVec, kernFcnPtrX, kernFcnPtrT));
    return rcpp_result_gen;
END_RCPP
}
// KRSmooth_matrix2
arma::mat KRSmooth_matrix2(arma::mat yMat, double h, SEXP kernFcnPtr);
RcppExport SEXP _DCSmooth_KRSmooth_matrix2(SEXP yMatSEXP, SEXP hSEXP, SEXP kernFcnPtrSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type yMat(yMatSEXP);
    Rcpp::traits::input_parameter< double >::type h(hSEXP);
    Rcpp::traits::input_parameter< SEXP >::type kernFcnPtr(kernFcnPtrSEXP);
    rcpp_result_gen = Rcpp::wrap(KRSmooth_matrix2(yMat, h, kernFcnPtr));
    return rcpp_result_gen;
END_RCPP
}
// KR_DoubleSmooth2
arma::mat KR_DoubleSmooth2(arma::mat yMat, arma::colvec hVec, arma::colvec drvVec, SEXP kernFcnPtrX, SEXP kernFcnPtrT);
RcppExport SEXP _DCSmooth_KR_DoubleSmooth2(SEXP yMatSEXP, SEXP hVecSEXP, SEXP drvVecSEXP, SEXP kernFcnPtrXSEXP, SEXP kernFcnPtrTSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type yMat(yMatSEXP);
    Rcpp::traits::input_parameter< arma::colvec >::type hVec(hVecSEXP);
    Rcpp::traits::input_parameter< arma::colvec >::type drvVec(drvVecSEXP);
    Rcpp::traits::input_parameter< SEXP >::type kernFcnPtrX(kernFcnPtrXSEXP);
    Rcpp::traits::input_parameter< SEXP >::type kernFcnPtrT(kernFcnPtrTSEXP);
    rcpp_result_gen = Rcpp::wrap(KR_DoubleSmooth2(yMat, hVec, drvVec, kernFcnPtrX, kernFcnPtrT));
    return rcpp_result_gen;
END_RCPP
}
// KRSmooth_matrix
arma::mat KRSmooth_matrix(arma::mat yMat, double h, int drv, SEXP kernFcnPtr);
RcppExport SEXP _DCSmooth_KRSmooth_matrix(SEXP yMatSEXP, SEXP hSEXP, SEXP drvSEXP, SEXP kernFcnPtrSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type yMat(yMatSEXP);
    Rcpp::traits::input_parameter< double >::type h(hSEXP);
    Rcpp::traits::input_parameter< int >::type drv(drvSEXP);
    Rcpp::traits::input_parameter< SEXP >::type kernFcnPtr(kernFcnPtrSEXP);
    rcpp_result_gen = Rcpp::wrap(KRSmooth_matrix(yMat, h, drv, kernFcnPtr));
    return rcpp_result_gen;
END_RCPP
}
// KR_DoubleSmooth
arma::mat KR_DoubleSmooth(arma::mat yMat, arma::colvec hVec, arma::icolvec drvVec, SEXP kernFcnPtrX, SEXP kernFcnPtrT);
RcppExport SEXP _DCSmooth_KR_DoubleSmooth(SEXP yMatSEXP, SEXP hVecSEXP, SEXP drvVecSEXP, SEXP kernFcnPtrXSEXP, SEXP kernFcnPtrTSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type yMat(yMatSEXP);
    Rcpp::traits::input_parameter< arma::colvec >::type hVec(hVecSEXP);
    Rcpp::traits::input_parameter< arma::icolvec >::type drvVec(drvVecSEXP);
    Rcpp::traits::input_parameter< SEXP >::type kernFcnPtrX(kernFcnPtrXSEXP);
    Rcpp::traits::input_parameter< SEXP >::type kernFcnPtrT(kernFcnPtrTSEXP);
    rcpp_result_gen = Rcpp::wrap(KR_DoubleSmooth(yMat, hVec, drvVec, kernFcnPtrX, kernFcnPtrT));
    return rcpp_result_gen;
END_RCPP
}
// kernFkt_MW220
arma::vec kernFkt_MW220(arma::vec& uVec, double q);
RcppExport SEXP _DCSmooth_kernFkt_MW220(SEXP uVecSEXP, SEXP qSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec& >::type uVec(uVecSEXP);
    Rcpp::traits::input_parameter< double >::type q(qSEXP);
    rcpp_result_gen = Rcpp::wrap(kernFkt_MW220(uVec, q));
    return rcpp_result_gen;
END_RCPP
}
// kernFkt_MW320
arma::vec kernFkt_MW320(arma::vec& u, double q);
RcppExport SEXP _DCSmooth_kernFkt_MW320(SEXP uSEXP, SEXP qSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec& >::type u(uSEXP);
    Rcpp::traits::input_parameter< double >::type q(qSEXP);
    rcpp_result_gen = Rcpp::wrap(kernFkt_MW320(u, q));
    return rcpp_result_gen;
END_RCPP
}
// kernFkt_MW420
arma::vec kernFkt_MW420(arma::vec& uVec, double q);
RcppExport SEXP _DCSmooth_kernFkt_MW420(SEXP uVecSEXP, SEXP qSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec& >::type uVec(uVecSEXP);
    Rcpp::traits::input_parameter< double >::type q(qSEXP);
    rcpp_result_gen = Rcpp::wrap(kernFkt_MW420(uVec, q));
    return rcpp_result_gen;
END_RCPP
}
// kernFkt_MW422
arma::vec kernFkt_MW422(arma::vec& uVec, double q);
RcppExport SEXP _DCSmooth_kernFkt_MW422(SEXP uVecSEXP, SEXP qSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec& >::type uVec(uVecSEXP);
    Rcpp::traits::input_parameter< double >::type q(qSEXP);
    rcpp_result_gen = Rcpp::wrap(kernFkt_MW422(uVec, q));
    return rcpp_result_gen;
END_RCPP
}
// kernFkt_MW422_1
arma::vec kernFkt_MW422_1(arma::vec& uVec);
RcppExport SEXP _DCSmooth_kernFkt_MW422_1(SEXP uVecSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec& >::type uVec(uVecSEXP);
    rcpp_result_gen = Rcpp::wrap(kernFkt_MW422_1(uVec));
    return rcpp_result_gen;
END_RCPP
}
// kernelFcn_assign
XPtr<funcPtr> kernelFcn_assign(std::string fstr);
RcppExport SEXP _DCSmooth_kernelFcn_assign(SEXP fstrSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string >::type fstr(fstrSEXP);
    rcpp_result_gen = Rcpp::wrap(kernelFcn_assign(fstr));
    return rcpp_result_gen;
END_RCPP
}
// kernelFcn_use
arma::vec kernelFcn_use(arma::vec x, double q, SEXP xpsexp);
RcppExport SEXP _DCSmooth_kernelFcn_use(SEXP xSEXP, SEXP qSEXP, SEXP xpsexpSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type x(xSEXP);
    Rcpp::traits::input_parameter< double >::type q(qSEXP);
    Rcpp::traits::input_parameter< SEXP >::type xpsexp(xpsexpSEXP);
    rcpp_result_gen = Rcpp::wrap(kernelFcn_use(x, q, xpsexp));
    return rcpp_result_gen;
END_RCPP
}
// MWTest1
arma::vec MWTest1(arma::vec uVec, double q);
RcppExport SEXP _DCSmooth_MWTest1(SEXP uVecSEXP, SEXP qSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type uVec(uVecSEXP);
    Rcpp::traits::input_parameter< double >::type q(qSEXP);
    rcpp_result_gen = Rcpp::wrap(MWTest1(uVec, q));
    return rcpp_result_gen;
END_RCPP
}
// MWTest2
double MWTest2(double u, double q);
RcppExport SEXP _DCSmooth_MWTest2(SEXP uSEXP, SEXP qSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type u(uSEXP);
    Rcpp::traits::input_parameter< double >::type q(qSEXP);
    rcpp_result_gen = Rcpp::wrap(MWTest2(u, q));
    return rcpp_result_gen;
END_RCPP
}
// LPSmooth_matrix3
arma::mat LPSmooth_matrix3(const arma::mat yMat, const double h, const int polyOrder, const int drv);
RcppExport SEXP _DCSmooth_LPSmooth_matrix3(SEXP yMatSEXP, SEXP hSEXP, SEXP polyOrderSEXP, SEXP drvSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat >::type yMat(yMatSEXP);
    Rcpp::traits::input_parameter< const double >::type h(hSEXP);
    Rcpp::traits::input_parameter< const int >::type polyOrder(polyOrderSEXP);
    Rcpp::traits::input_parameter< const int >::type drv(drvSEXP);
    rcpp_result_gen = Rcpp::wrap(LPSmooth_matrix3(yMat, h, polyOrder, drv));
    return rcpp_result_gen;
END_RCPP
}
// LPSmooth_matrix2
arma::mat LPSmooth_matrix2(const arma::mat yMat, const double h, const int polyOrder, const int drv);
RcppExport SEXP _DCSmooth_LPSmooth_matrix2(SEXP yMatSEXP, SEXP hSEXP, SEXP polyOrderSEXP, SEXP drvSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat >::type yMat(yMatSEXP);
    Rcpp::traits::input_parameter< const double >::type h(hSEXP);
    Rcpp::traits::input_parameter< const int >::type polyOrder(polyOrderSEXP);
    Rcpp::traits::input_parameter< const int >::type drv(drvSEXP);
    rcpp_result_gen = Rcpp::wrap(LPSmooth_matrix2(yMat, h, polyOrder, drv));
    return rcpp_result_gen;
END_RCPP
}
// LP_DoubleSmooth2
arma::mat LP_DoubleSmooth2(arma::mat yMat, arma::colvec hVec, arma::icolvec polyOrderVec, arma::icolvec drvVec);
RcppExport SEXP _DCSmooth_LP_DoubleSmooth2(SEXP yMatSEXP, SEXP hVecSEXP, SEXP polyOrderVecSEXP, SEXP drvVecSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type yMat(yMatSEXP);
    Rcpp::traits::input_parameter< arma::colvec >::type hVec(hVecSEXP);
    Rcpp::traits::input_parameter< arma::icolvec >::type polyOrderVec(polyOrderVecSEXP);
    Rcpp::traits::input_parameter< arma::icolvec >::type drvVec(drvVecSEXP);
    rcpp_result_gen = Rcpp::wrap(LP_DoubleSmooth2(yMat, hVec, polyOrderVec, drvVec));
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
// LP_DoubleSmooth
arma::mat LP_DoubleSmooth(arma::mat yMat, arma::colvec hVec, arma::icolvec polyOrderVec, arma::icolvec drvVec);
RcppExport SEXP _DCSmooth_LP_DoubleSmooth(SEXP yMatSEXP, SEXP hVecSEXP, SEXP polyOrderVecSEXP, SEXP drvVecSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type yMat(yMatSEXP);
    Rcpp::traits::input_parameter< arma::colvec >::type hVec(hVecSEXP);
    Rcpp::traits::input_parameter< arma::icolvec >::type polyOrderVec(polyOrderVecSEXP);
    Rcpp::traits::input_parameter< arma::icolvec >::type drvVec(drvVecSEXP);
    rcpp_result_gen = Rcpp::wrap(LP_DoubleSmooth(yMat, hVec, polyOrderVec, drvVec));
    return rcpp_result_gen;
END_RCPP
}
// cppSample
int cppSample(arma::vec x);
RcppExport SEXP _DCSmooth_cppSample(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(cppSample(x));
    return rcpp_result_gen;
END_RCPP
}
// thinnedMat
arma::mat thinnedMat(arma::mat yMat, int seed);
RcppExport SEXP _DCSmooth_thinnedMat(SEXP yMatSEXP, SEXP seedSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type yMat(yMatSEXP);
    Rcpp::traits::input_parameter< int >::type seed(seedSEXP);
    rcpp_result_gen = Rcpp::wrap(thinnedMat(yMat, seed));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_DCSmooth_KRSlow_single", (DL_FUNC) &_DCSmooth_KRSlow_single, 3},
    {"_DCSmooth_KR_DSSlow", (DL_FUNC) &_DCSmooth_KR_DSSlow, 5},
    {"_DCSmooth_KRSmooth_matrix2", (DL_FUNC) &_DCSmooth_KRSmooth_matrix2, 3},
    {"_DCSmooth_KR_DoubleSmooth2", (DL_FUNC) &_DCSmooth_KR_DoubleSmooth2, 5},
    {"_DCSmooth_KRSmooth_matrix", (DL_FUNC) &_DCSmooth_KRSmooth_matrix, 4},
    {"_DCSmooth_KR_DoubleSmooth", (DL_FUNC) &_DCSmooth_KR_DoubleSmooth, 5},
    {"_DCSmooth_kernFkt_MW220", (DL_FUNC) &_DCSmooth_kernFkt_MW220, 2},
    {"_DCSmooth_kernFkt_MW320", (DL_FUNC) &_DCSmooth_kernFkt_MW320, 2},
    {"_DCSmooth_kernFkt_MW420", (DL_FUNC) &_DCSmooth_kernFkt_MW420, 2},
    {"_DCSmooth_kernFkt_MW422", (DL_FUNC) &_DCSmooth_kernFkt_MW422, 2},
    {"_DCSmooth_kernFkt_MW422_1", (DL_FUNC) &_DCSmooth_kernFkt_MW422_1, 1},
    {"_DCSmooth_kernelFcn_assign", (DL_FUNC) &_DCSmooth_kernelFcn_assign, 1},
    {"_DCSmooth_kernelFcn_use", (DL_FUNC) &_DCSmooth_kernelFcn_use, 3},
    {"_DCSmooth_MWTest1", (DL_FUNC) &_DCSmooth_MWTest1, 2},
    {"_DCSmooth_MWTest2", (DL_FUNC) &_DCSmooth_MWTest2, 2},
    {"_DCSmooth_LPSmooth_matrix3", (DL_FUNC) &_DCSmooth_LPSmooth_matrix3, 4},
    {"_DCSmooth_LPSmooth_matrix2", (DL_FUNC) &_DCSmooth_LPSmooth_matrix2, 4},
    {"_DCSmooth_LP_DoubleSmooth2", (DL_FUNC) &_DCSmooth_LP_DoubleSmooth2, 4},
    {"_DCSmooth_LPSmooth_matrix", (DL_FUNC) &_DCSmooth_LPSmooth_matrix, 4},
    {"_DCSmooth_LP_DoubleSmooth", (DL_FUNC) &_DCSmooth_LP_DoubleSmooth, 4},
    {"_DCSmooth_cppSample", (DL_FUNC) &_DCSmooth_cppSample, 1},
    {"_DCSmooth_thinnedMat", (DL_FUNC) &_DCSmooth_thinnedMat, 2},
    {NULL, NULL, 0}
};

RcppExport void R_init_DCSmooth(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
