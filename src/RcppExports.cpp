// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include "DCSmooth_types.h"
#include <RcppArmadillo.h>
#include <RcppEigen.h>
#include <Rcpp.h>

using namespace Rcpp;

// xMatrix
arma::mat xMatrix(arma::colvec xVector, int polyOrder);
RcppExport SEXP _DCSmooth_xMatrix(SEXP xVectorSEXP, SEXP polyOrderSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::colvec >::type xVector(xVectorSEXP);
    Rcpp::traits::input_parameter< int >::type polyOrder(polyOrderSEXP);
    rcpp_result_gen = Rcpp::wrap(xMatrix(xVector, polyOrder));
    return rcpp_result_gen;
END_RCPP
}
// np_matrix
arma::mat np_matrix(SEXP kernFcnPtr, int p, int n);
RcppExport SEXP _DCSmooth_np_matrix(SEXP kernFcnPtrSEXP, SEXP pSEXP, SEXP nSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type kernFcnPtr(kernFcnPtrSEXP);
    Rcpp::traits::input_parameter< int >::type p(pSEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    rcpp_result_gen = Rcpp::wrap(np_matrix(kernFcnPtr, p, n));
    return rcpp_result_gen;
END_RCPP
}
// m_weights
arma::vec m_weights(arma::mat npMatrix, arma::vec u, int drv);
RcppExport SEXP _DCSmooth_m_weights(SEXP npMatrixSEXP, SEXP uSEXP, SEXP drvSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type npMatrix(npMatrixSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type u(uSEXP);
    Rcpp::traits::input_parameter< int >::type drv(drvSEXP);
    rcpp_result_gen = Rcpp::wrap(m_weights(npMatrix, u, drv));
    return rcpp_result_gen;
END_RCPP
}
// acfMatrix_quarter2
arma::mat acfMatrix_quarter2(const arma::mat y_Mat);
RcppExport SEXP _DCSmooth_acfMatrix_quarter2(SEXP y_MatSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat >::type y_Mat(y_MatSEXP);
    rcpp_result_gen = Rcpp::wrap(acfMatrix_quarter2(y_Mat));
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
// KR_dcs_const0
arma::mat KR_dcs_const0(arma::mat yMat, arma::colvec hVec, arma::colvec drvVec, SEXP kernFcnPtrX, SEXP kernFcnPtrT);
RcppExport SEXP _DCSmooth_KR_dcs_const0(SEXP yMatSEXP, SEXP hVecSEXP, SEXP drvVecSEXP, SEXP kernFcnPtrXSEXP, SEXP kernFcnPtrTSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type yMat(yMatSEXP);
    Rcpp::traits::input_parameter< arma::colvec >::type hVec(hVecSEXP);
    Rcpp::traits::input_parameter< arma::colvec >::type drvVec(drvVecSEXP);
    Rcpp::traits::input_parameter< SEXP >::type kernFcnPtrX(kernFcnPtrXSEXP);
    Rcpp::traits::input_parameter< SEXP >::type kernFcnPtrT(kernFcnPtrTSEXP);
    rcpp_result_gen = Rcpp::wrap(KR_dcs_const0(yMat, hVec, drvVec, kernFcnPtrX, kernFcnPtrT));
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
// KR_dcs_const1
arma::mat KR_dcs_const1(arma::mat yMat, arma::colvec hVec, arma::icolvec drvVec, SEXP kernFcnPtrX, SEXP kernFcnPtrT);
RcppExport SEXP _DCSmooth_KR_dcs_const1(SEXP yMatSEXP, SEXP hVecSEXP, SEXP drvVecSEXP, SEXP kernFcnPtrXSEXP, SEXP kernFcnPtrTSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type yMat(yMatSEXP);
    Rcpp::traits::input_parameter< arma::colvec >::type hVec(hVecSEXP);
    Rcpp::traits::input_parameter< arma::icolvec >::type drvVec(drvVecSEXP);
    Rcpp::traits::input_parameter< SEXP >::type kernFcnPtrX(kernFcnPtrXSEXP);
    Rcpp::traits::input_parameter< SEXP >::type kernFcnPtrT(kernFcnPtrTSEXP);
    rcpp_result_gen = Rcpp::wrap(KR_dcs_const1(yMat, hVec, drvVec, kernFcnPtrX, kernFcnPtrT));
    return rcpp_result_gen;
END_RCPP
}
// kernFkt_MW210
arma::vec kernFkt_MW210(arma::vec& u, double q);
RcppExport SEXP _DCSmooth_kernFkt_MW210(SEXP uSEXP, SEXP qSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec& >::type u(uSEXP);
    Rcpp::traits::input_parameter< double >::type q(qSEXP);
    rcpp_result_gen = Rcpp::wrap(kernFkt_MW210(u, q));
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
// kernFkt_MW421
arma::vec kernFkt_MW421(arma::vec& uVec, double q);
RcppExport SEXP _DCSmooth_kernFkt_MW421(SEXP uVecSEXP, SEXP qSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec& >::type uVec(uVecSEXP);
    Rcpp::traits::input_parameter< double >::type q(qSEXP);
    rcpp_result_gen = Rcpp::wrap(kernFkt_MW421(uVec, q));
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
// kernFkt_TR420
arma::vec kernFkt_TR420(arma::vec& uVec, double q);
RcppExport SEXP _DCSmooth_kernFkt_TR420(SEXP uVecSEXP, SEXP qSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec& >::type uVec(uVecSEXP);
    Rcpp::traits::input_parameter< double >::type q(qSEXP);
    rcpp_result_gen = Rcpp::wrap(kernFkt_TR420(uVec, q));
    return rcpp_result_gen;
END_RCPP
}
// kernFkt_TR422
arma::vec kernFkt_TR422(arma::vec& uVec, double q);
RcppExport SEXP _DCSmooth_kernFkt_TR422(SEXP uVecSEXP, SEXP qSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec& >::type uVec(uVecSEXP);
    Rcpp::traits::input_parameter< double >::type q(qSEXP);
    rcpp_result_gen = Rcpp::wrap(kernFkt_TR422(uVec, q));
    return rcpp_result_gen;
END_RCPP
}
// kernel_fcn_assign
XPtr<funcPtr> kernel_fcn_assign(std::string fstr);
RcppExport SEXP _DCSmooth_kernel_fcn_assign(SEXP fstrSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string >::type fstr(fstrSEXP);
    rcpp_result_gen = Rcpp::wrap(kernel_fcn_assign(fstr));
    return rcpp_result_gen;
END_RCPP
}
// kernel_fcn_use
arma::vec kernel_fcn_use(arma::vec x, double q, SEXP xpsexp);
RcppExport SEXP _DCSmooth_kernel_fcn_use(SEXP xSEXP, SEXP qSEXP, SEXP xpsexpSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type x(xSEXP);
    Rcpp::traits::input_parameter< double >::type q(qSEXP);
    Rcpp::traits::input_parameter< SEXP >::type xpsexp(xpsexpSEXP);
    rcpp_result_gen = Rcpp::wrap(kernel_fcn_use(x, q, xpsexp));
    return rcpp_result_gen;
END_RCPP
}
// LPSmooth_matrix2
arma::mat LPSmooth_matrix2(const arma::mat yMat, const double h, const int polyOrder, const int drv, SEXP kernFcnPtr);
RcppExport SEXP _DCSmooth_LPSmooth_matrix2(SEXP yMatSEXP, SEXP hSEXP, SEXP polyOrderSEXP, SEXP drvSEXP, SEXP kernFcnPtrSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat >::type yMat(yMatSEXP);
    Rcpp::traits::input_parameter< const double >::type h(hSEXP);
    Rcpp::traits::input_parameter< const int >::type polyOrder(polyOrderSEXP);
    Rcpp::traits::input_parameter< const int >::type drv(drvSEXP);
    Rcpp::traits::input_parameter< SEXP >::type kernFcnPtr(kernFcnPtrSEXP);
    rcpp_result_gen = Rcpp::wrap(LPSmooth_matrix2(yMat, h, polyOrder, drv, kernFcnPtr));
    return rcpp_result_gen;
END_RCPP
}
// LP_dcs_const0
arma::mat LP_dcs_const0(arma::mat yMat, arma::colvec hVec, arma::icolvec polyOrderVec, arma::icolvec drvVec, SEXP kernFcnPtr_x, SEXP kernFcnPtr_t);
RcppExport SEXP _DCSmooth_LP_dcs_const0(SEXP yMatSEXP, SEXP hVecSEXP, SEXP polyOrderVecSEXP, SEXP drvVecSEXP, SEXP kernFcnPtr_xSEXP, SEXP kernFcnPtr_tSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type yMat(yMatSEXP);
    Rcpp::traits::input_parameter< arma::colvec >::type hVec(hVecSEXP);
    Rcpp::traits::input_parameter< arma::icolvec >::type polyOrderVec(polyOrderVecSEXP);
    Rcpp::traits::input_parameter< arma::icolvec >::type drvVec(drvVecSEXP);
    Rcpp::traits::input_parameter< SEXP >::type kernFcnPtr_x(kernFcnPtr_xSEXP);
    Rcpp::traits::input_parameter< SEXP >::type kernFcnPtr_t(kernFcnPtr_tSEXP);
    rcpp_result_gen = Rcpp::wrap(LP_dcs_const0(yMat, hVec, polyOrderVec, drvVec, kernFcnPtr_x, kernFcnPtr_t));
    return rcpp_result_gen;
END_RCPP
}
// LPSmooth_matrix2_BMod
arma::mat LPSmooth_matrix2_BMod(const arma::mat yMat, const double h, const int polyOrder, const int drv, SEXP kernFcnPtr);
RcppExport SEXP _DCSmooth_LPSmooth_matrix2_BMod(SEXP yMatSEXP, SEXP hSEXP, SEXP polyOrderSEXP, SEXP drvSEXP, SEXP kernFcnPtrSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat >::type yMat(yMatSEXP);
    Rcpp::traits::input_parameter< const double >::type h(hSEXP);
    Rcpp::traits::input_parameter< const int >::type polyOrder(polyOrderSEXP);
    Rcpp::traits::input_parameter< const int >::type drv(drvSEXP);
    Rcpp::traits::input_parameter< SEXP >::type kernFcnPtr(kernFcnPtrSEXP);
    rcpp_result_gen = Rcpp::wrap(LPSmooth_matrix2_BMod(yMat, h, polyOrder, drv, kernFcnPtr));
    return rcpp_result_gen;
END_RCPP
}
// LP_dcs_const0_BMod
arma::mat LP_dcs_const0_BMod(arma::mat yMat, arma::colvec hVec, arma::colvec polyOrderVec, arma::icolvec drvVec, SEXP kernFcnPtr_1, SEXP kernFcnPtr_2);
RcppExport SEXP _DCSmooth_LP_dcs_const0_BMod(SEXP yMatSEXP, SEXP hVecSEXP, SEXP polyOrderVecSEXP, SEXP drvVecSEXP, SEXP kernFcnPtr_1SEXP, SEXP kernFcnPtr_2SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type yMat(yMatSEXP);
    Rcpp::traits::input_parameter< arma::colvec >::type hVec(hVecSEXP);
    Rcpp::traits::input_parameter< arma::colvec >::type polyOrderVec(polyOrderVecSEXP);
    Rcpp::traits::input_parameter< arma::icolvec >::type drvVec(drvVecSEXP);
    Rcpp::traits::input_parameter< SEXP >::type kernFcnPtr_1(kernFcnPtr_1SEXP);
    Rcpp::traits::input_parameter< SEXP >::type kernFcnPtr_2(kernFcnPtr_2SEXP);
    rcpp_result_gen = Rcpp::wrap(LP_dcs_const0_BMod(yMat, hVec, polyOrderVec, drvVec, kernFcnPtr_1, kernFcnPtr_2));
    return rcpp_result_gen;
END_RCPP
}
// LPSmooth_matrix
arma::mat LPSmooth_matrix(const arma::mat yMat, const double h, const int polyOrder, const int drv, SEXP kernFcnPtr);
RcppExport SEXP _DCSmooth_LPSmooth_matrix(SEXP yMatSEXP, SEXP hSEXP, SEXP polyOrderSEXP, SEXP drvSEXP, SEXP kernFcnPtrSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat >::type yMat(yMatSEXP);
    Rcpp::traits::input_parameter< const double >::type h(hSEXP);
    Rcpp::traits::input_parameter< const int >::type polyOrder(polyOrderSEXP);
    Rcpp::traits::input_parameter< const int >::type drv(drvSEXP);
    Rcpp::traits::input_parameter< SEXP >::type kernFcnPtr(kernFcnPtrSEXP);
    rcpp_result_gen = Rcpp::wrap(LPSmooth_matrix(yMat, h, polyOrder, drv, kernFcnPtr));
    return rcpp_result_gen;
END_RCPP
}
// LP_dcs_const1
arma::mat LP_dcs_const1(arma::mat yMat, arma::colvec hVec, arma::icolvec polyOrderVec, arma::icolvec drvVec, SEXP kernFcnPtr_x, SEXP kernFcnPtr_t);
RcppExport SEXP _DCSmooth_LP_dcs_const1(SEXP yMatSEXP, SEXP hVecSEXP, SEXP polyOrderVecSEXP, SEXP drvVecSEXP, SEXP kernFcnPtr_xSEXP, SEXP kernFcnPtr_tSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type yMat(yMatSEXP);
    Rcpp::traits::input_parameter< arma::colvec >::type hVec(hVecSEXP);
    Rcpp::traits::input_parameter< arma::icolvec >::type polyOrderVec(polyOrderVecSEXP);
    Rcpp::traits::input_parameter< arma::icolvec >::type drvVec(drvVecSEXP);
    Rcpp::traits::input_parameter< SEXP >::type kernFcnPtr_x(kernFcnPtr_xSEXP);
    Rcpp::traits::input_parameter< SEXP >::type kernFcnPtr_t(kernFcnPtr_tSEXP);
    rcpp_result_gen = Rcpp::wrap(LP_dcs_const1(yMat, hVec, polyOrderVec, drvVec, kernFcnPtr_x, kernFcnPtr_t));
    return rcpp_result_gen;
END_RCPP
}
// LPSmooth_matrix_BMod
arma::mat LPSmooth_matrix_BMod(const arma::mat yMat, const double h, const int polyOrder, const int drv, SEXP kernFcnPtr);
RcppExport SEXP _DCSmooth_LPSmooth_matrix_BMod(SEXP yMatSEXP, SEXP hSEXP, SEXP polyOrderSEXP, SEXP drvSEXP, SEXP kernFcnPtrSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat >::type yMat(yMatSEXP);
    Rcpp::traits::input_parameter< const double >::type h(hSEXP);
    Rcpp::traits::input_parameter< const int >::type polyOrder(polyOrderSEXP);
    Rcpp::traits::input_parameter< const int >::type drv(drvSEXP);
    Rcpp::traits::input_parameter< SEXP >::type kernFcnPtr(kernFcnPtrSEXP);
    rcpp_result_gen = Rcpp::wrap(LPSmooth_matrix_BMod(yMat, h, polyOrder, drv, kernFcnPtr));
    return rcpp_result_gen;
END_RCPP
}
// LP_dcs_const1_BMod
arma::mat LP_dcs_const1_BMod(arma::mat yMat, arma::colvec hVec, arma::icolvec polyOrderVec, arma::icolvec drvVec, SEXP kernFcnPtr_x, SEXP kernFcnPtr_t);
RcppExport SEXP _DCSmooth_LP_dcs_const1_BMod(SEXP yMatSEXP, SEXP hVecSEXP, SEXP polyOrderVecSEXP, SEXP drvVecSEXP, SEXP kernFcnPtr_xSEXP, SEXP kernFcnPtr_tSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type yMat(yMatSEXP);
    Rcpp::traits::input_parameter< arma::colvec >::type hVec(hVecSEXP);
    Rcpp::traits::input_parameter< arma::icolvec >::type polyOrderVec(polyOrderVecSEXP);
    Rcpp::traits::input_parameter< arma::icolvec >::type drvVec(drvVecSEXP);
    Rcpp::traits::input_parameter< SEXP >::type kernFcnPtr_x(kernFcnPtr_xSEXP);
    Rcpp::traits::input_parameter< SEXP >::type kernFcnPtr_t(kernFcnPtr_tSEXP);
    rcpp_result_gen = Rcpp::wrap(LP_dcs_const1_BMod(yMat, hVec, polyOrderVec, drvVec, kernFcnPtr_x, kernFcnPtr_t));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_DCSmooth_xMatrix", (DL_FUNC) &_DCSmooth_xMatrix, 2},
    {"_DCSmooth_np_matrix", (DL_FUNC) &_DCSmooth_np_matrix, 3},
    {"_DCSmooth_m_weights", (DL_FUNC) &_DCSmooth_m_weights, 3},
    {"_DCSmooth_acfMatrix_quarter2", (DL_FUNC) &_DCSmooth_acfMatrix_quarter2, 1},
    {"_DCSmooth_KRSmooth_matrix2", (DL_FUNC) &_DCSmooth_KRSmooth_matrix2, 3},
    {"_DCSmooth_KR_dcs_const0", (DL_FUNC) &_DCSmooth_KR_dcs_const0, 5},
    {"_DCSmooth_KRSmooth_matrix", (DL_FUNC) &_DCSmooth_KRSmooth_matrix, 4},
    {"_DCSmooth_KR_dcs_const1", (DL_FUNC) &_DCSmooth_KR_dcs_const1, 5},
    {"_DCSmooth_kernFkt_MW210", (DL_FUNC) &_DCSmooth_kernFkt_MW210, 2},
    {"_DCSmooth_kernFkt_MW220", (DL_FUNC) &_DCSmooth_kernFkt_MW220, 2},
    {"_DCSmooth_kernFkt_MW320", (DL_FUNC) &_DCSmooth_kernFkt_MW320, 2},
    {"_DCSmooth_kernFkt_MW420", (DL_FUNC) &_DCSmooth_kernFkt_MW420, 2},
    {"_DCSmooth_kernFkt_MW421", (DL_FUNC) &_DCSmooth_kernFkt_MW421, 2},
    {"_DCSmooth_kernFkt_MW422", (DL_FUNC) &_DCSmooth_kernFkt_MW422, 2},
    {"_DCSmooth_kernFkt_TR420", (DL_FUNC) &_DCSmooth_kernFkt_TR420, 2},
    {"_DCSmooth_kernFkt_TR422", (DL_FUNC) &_DCSmooth_kernFkt_TR422, 2},
    {"_DCSmooth_kernel_fcn_assign", (DL_FUNC) &_DCSmooth_kernel_fcn_assign, 1},
    {"_DCSmooth_kernel_fcn_use", (DL_FUNC) &_DCSmooth_kernel_fcn_use, 3},
    {"_DCSmooth_LPSmooth_matrix2", (DL_FUNC) &_DCSmooth_LPSmooth_matrix2, 5},
    {"_DCSmooth_LP_dcs_const0", (DL_FUNC) &_DCSmooth_LP_dcs_const0, 6},
    {"_DCSmooth_LPSmooth_matrix2_BMod", (DL_FUNC) &_DCSmooth_LPSmooth_matrix2_BMod, 5},
    {"_DCSmooth_LP_dcs_const0_BMod", (DL_FUNC) &_DCSmooth_LP_dcs_const0_BMod, 6},
    {"_DCSmooth_LPSmooth_matrix", (DL_FUNC) &_DCSmooth_LPSmooth_matrix, 5},
    {"_DCSmooth_LP_dcs_const1", (DL_FUNC) &_DCSmooth_LP_dcs_const1, 6},
    {"_DCSmooth_LPSmooth_matrix_BMod", (DL_FUNC) &_DCSmooth_LPSmooth_matrix_BMod, 5},
    {"_DCSmooth_LP_dcs_const1_BMod", (DL_FUNC) &_DCSmooth_LP_dcs_const1_BMod, 6},
    {NULL, NULL, 0}
};

RcppExport void R_init_DCSmooth(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
