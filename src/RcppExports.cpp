// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include "DCSmooth_types.h"
#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

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
// cumsum_part_reverse
NumericVector cumsum_part_reverse(arma::rowvec vec_1, arma::colvec vec_2);
RcppExport SEXP _DCSmooth_cumsum_part_reverse(SEXP vec_1SEXP, SEXP vec_2SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::rowvec >::type vec_1(vec_1SEXP);
    Rcpp::traits::input_parameter< arma::colvec >::type vec_2(vec_2SEXP);
    rcpp_result_gen = Rcpp::wrap(cumsum_part_reverse(vec_1, vec_2));
    return rcpp_result_gen;
END_RCPP
}
// specDens_cpp
double specDens_cpp(arma::mat Y, arma::vec omega);
RcppExport SEXP _DCSmooth_specDens_cpp(SEXP YSEXP, SEXP omegaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type Y(YSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type omega(omegaSEXP);
    rcpp_result_gen = Rcpp::wrap(specDens_cpp(Y, omega));
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
// acfMatrix_quarter
arma::mat acfMatrix_quarter(const arma::mat y_Mat);
RcppExport SEXP _DCSmooth_acfMatrix_quarter(SEXP y_MatSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat >::type y_Mat(y_MatSEXP);
    rcpp_result_gen = Rcpp::wrap(acfMatrix_quarter(y_Mat));
    return rcpp_result_gen;
END_RCPP
}
// acfMatrix_cpp
arma::mat acfMatrix_cpp(const arma::mat y_Mat);
RcppExport SEXP _DCSmooth_acfMatrix_cpp(SEXP y_MatSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat >::type y_Mat(y_MatSEXP);
    rcpp_result_gen = Rcpp::wrap(acfMatrix_cpp(y_Mat));
    return rcpp_result_gen;
END_RCPP
}
// globalBndw_cpp
arma::vec globalBndw_cpp(const arma::mat acfMat, arma::vec hLag);
RcppExport SEXP _DCSmooth_globalBndw_cpp(SEXP acfMatSEXP, SEXP hLagSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat >::type acfMat(acfMatSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type hLag(hLagSEXP);
    rcpp_result_gen = Rcpp::wrap(globalBndw_cpp(acfMat, hLag));
    return rcpp_result_gen;
END_RCPP
}
// specInt00_cpp
double specInt00_cpp(arma::mat acfMat);
RcppExport SEXP _DCSmooth_specInt00_cpp(SEXP acfMatSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type acfMat(acfMatSEXP);
    rcpp_result_gen = Rcpp::wrap(specInt00_cpp(acfMat));
    return rcpp_result_gen;
END_RCPP
}
// specIntDrv_cpp
double specIntDrv_cpp(arma::mat acfMat, arma::vec drv1, arma::vec drv2, arma::vec hLag);
RcppExport SEXP _DCSmooth_specIntDrv_cpp(SEXP acfMatSEXP, SEXP drv1SEXP, SEXP drv2SEXP, SEXP hLagSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type acfMat(acfMatSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type drv1(drv1SEXP);
    Rcpp::traits::input_parameter< arma::vec >::type drv2(drv2SEXP);
    Rcpp::traits::input_parameter< arma::vec >::type hLag(hLagSEXP);
    rcpp_result_gen = Rcpp::wrap(specIntDrv_cpp(acfMat, drv1, drv2, hLag));
    return rcpp_result_gen;
END_RCPP
}
// localBndw_cpp
arma::vec localBndw_cpp(const arma::mat acfMat, arma::vec hLag, arma::vec omega);
RcppExport SEXP _DCSmooth_localBndw_cpp(SEXP acfMatSEXP, SEXP hLagSEXP, SEXP omegaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat >::type acfMat(acfMatSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type hLag(hLagSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type omega(omegaSEXP);
    rcpp_result_gen = Rcpp::wrap(localBndw_cpp(acfMat, hLag, omega));
    return rcpp_result_gen;
END_RCPP
}
// specDensEst_cpp
double specDensEst_cpp(arma::mat acfMat, arma::vec drv, arma::vec hLag, arma::vec omega);
RcppExport SEXP _DCSmooth_specDensEst_cpp(SEXP acfMatSEXP, SEXP drvSEXP, SEXP hLagSEXP, SEXP omegaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type acfMat(acfMatSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type drv(drvSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type hLag(hLagSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type omega(omegaSEXP);
    rcpp_result_gen = Rcpp::wrap(specDensEst_cpp(acfMat, drv, hLag, omega));
    return rcpp_result_gen;
END_RCPP
}
// fourier_mat
arma::cx_mat fourier_mat(int Lx, int Lt, arma::vec omega);
RcppExport SEXP _DCSmooth_fourier_mat(SEXP LxSEXP, SEXP LtSEXP, SEXP omegaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type Lx(LxSEXP);
    Rcpp::traits::input_parameter< int >::type Lt(LtSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type omega(omegaSEXP);
    rcpp_result_gen = Rcpp::wrap(fourier_mat(Lx, Lt, omega));
    return rcpp_result_gen;
END_RCPP
}
// bartlett_mat
arma::mat bartlett_mat(int Lx, int Lt, arma::vec drv1, arma::vec drv2, int power);
RcppExport SEXP _DCSmooth_bartlett_mat(SEXP LxSEXP, SEXP LtSEXP, SEXP drv1SEXP, SEXP drv2SEXP, SEXP powerSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type Lx(LxSEXP);
    Rcpp::traits::input_parameter< int >::type Lt(LtSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type drv1(drv1SEXP);
    Rcpp::traits::input_parameter< arma::vec >::type drv2(drv2SEXP);
    Rcpp::traits::input_parameter< int >::type power(powerSEXP);
    rcpp_result_gen = Rcpp::wrap(bartlett_mat(Lx, Lt, drv1, drv2, power));
    return rcpp_result_gen;
END_RCPP
}
// KRSmooth_matrix2
arma::mat KRSmooth_matrix2(arma::mat yMat, double h, int drv, SEXP kernFcnPtr);
RcppExport SEXP _DCSmooth_KRSmooth_matrix2(SEXP yMatSEXP, SEXP hSEXP, SEXP drvSEXP, SEXP kernFcnPtrSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type yMat(yMatSEXP);
    Rcpp::traits::input_parameter< double >::type h(hSEXP);
    Rcpp::traits::input_parameter< int >::type drv(drvSEXP);
    Rcpp::traits::input_parameter< SEXP >::type kernFcnPtr(kernFcnPtrSEXP);
    rcpp_result_gen = Rcpp::wrap(KRSmooth_matrix2(yMat, h, drv, kernFcnPtr));
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
// kernel_fcn_assign
Rcpp::XPtr<funcPtr> kernel_fcn_assign(std::string fstr);
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
// weight_fcn_assign
Rcpp::XPtr<weightPtr> weight_fcn_assign(std::string fstr);
RcppExport SEXP _DCSmooth_weight_fcn_assign(SEXP fstrSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string >::type fstr(fstrSEXP);
    rcpp_result_gen = Rcpp::wrap(weight_fcn_assign(fstr));
    return rcpp_result_gen;
END_RCPP
}
// kern_fcn_M200
arma::vec kern_fcn_M200(arma::vec& uVec, double q);
RcppExport SEXP _DCSmooth_kern_fcn_M200(SEXP uVecSEXP, SEXP qSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec& >::type uVec(uVecSEXP);
    Rcpp::traits::input_parameter< double >::type q(qSEXP);
    rcpp_result_gen = Rcpp::wrap(kern_fcn_M200(uVec, q));
    return rcpp_result_gen;
END_RCPP
}
// kern_fcn_M210
arma::vec kern_fcn_M210(arma::vec& uVec, double q);
RcppExport SEXP _DCSmooth_kern_fcn_M210(SEXP uVecSEXP, SEXP qSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec& >::type uVec(uVecSEXP);
    Rcpp::traits::input_parameter< double >::type q(qSEXP);
    rcpp_result_gen = Rcpp::wrap(kern_fcn_M210(uVec, q));
    return rcpp_result_gen;
END_RCPP
}
// kern_fcn_M220
arma::vec kern_fcn_M220(arma::vec& uVec, double q);
RcppExport SEXP _DCSmooth_kern_fcn_M220(SEXP uVecSEXP, SEXP qSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec& >::type uVec(uVecSEXP);
    Rcpp::traits::input_parameter< double >::type q(qSEXP);
    rcpp_result_gen = Rcpp::wrap(kern_fcn_M220(uVec, q));
    return rcpp_result_gen;
END_RCPP
}
// kern_fcn_M321
arma::vec kern_fcn_M321(arma::vec& uVec, double q);
RcppExport SEXP _DCSmooth_kern_fcn_M321(SEXP uVecSEXP, SEXP qSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec& >::type uVec(uVecSEXP);
    Rcpp::traits::input_parameter< double >::type q(qSEXP);
    rcpp_result_gen = Rcpp::wrap(kern_fcn_M321(uVec, q));
    return rcpp_result_gen;
END_RCPP
}
// kern_fcn_M420
arma::vec kern_fcn_M420(arma::vec& uVec, double q);
RcppExport SEXP _DCSmooth_kern_fcn_M420(SEXP uVecSEXP, SEXP qSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec& >::type uVec(uVecSEXP);
    Rcpp::traits::input_parameter< double >::type q(qSEXP);
    rcpp_result_gen = Rcpp::wrap(kern_fcn_M420(uVec, q));
    return rcpp_result_gen;
END_RCPP
}
// kern_fcn_M421
arma::vec kern_fcn_M421(arma::vec& uVec, double q);
RcppExport SEXP _DCSmooth_kern_fcn_M421(SEXP uVecSEXP, SEXP qSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec& >::type uVec(uVecSEXP);
    Rcpp::traits::input_parameter< double >::type q(qSEXP);
    rcpp_result_gen = Rcpp::wrap(kern_fcn_M421(uVec, q));
    return rcpp_result_gen;
END_RCPP
}
// kern_fcn_M422
arma::vec kern_fcn_M422(arma::vec& uVec, double q);
RcppExport SEXP _DCSmooth_kern_fcn_M422(SEXP uVecSEXP, SEXP qSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec& >::type uVec(uVecSEXP);
    Rcpp::traits::input_parameter< double >::type q(qSEXP);
    rcpp_result_gen = Rcpp::wrap(kern_fcn_M422(uVec, q));
    return rcpp_result_gen;
END_RCPP
}
// kern_fcn_M533
arma::vec kern_fcn_M533(arma::vec& uVec, double q);
RcppExport SEXP _DCSmooth_kern_fcn_M533(SEXP uVecSEXP, SEXP qSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec& >::type uVec(uVecSEXP);
    Rcpp::traits::input_parameter< double >::type q(qSEXP);
    rcpp_result_gen = Rcpp::wrap(kern_fcn_M533(uVec, q));
    return rcpp_result_gen;
END_RCPP
}
// kern_fcn_M644
arma::vec kern_fcn_M644(arma::vec& uVec, double q);
RcppExport SEXP _DCSmooth_kern_fcn_M644(SEXP uVecSEXP, SEXP qSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec& >::type uVec(uVecSEXP);
    Rcpp::traits::input_parameter< double >::type q(qSEXP);
    rcpp_result_gen = Rcpp::wrap(kern_fcn_M644(uVec, q));
    return rcpp_result_gen;
END_RCPP
}
// kern_fcn_MW200
arma::vec kern_fcn_MW200(arma::vec& uVec, double q);
RcppExport SEXP _DCSmooth_kern_fcn_MW200(SEXP uVecSEXP, SEXP qSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec& >::type uVec(uVecSEXP);
    Rcpp::traits::input_parameter< double >::type q(qSEXP);
    rcpp_result_gen = Rcpp::wrap(kern_fcn_MW200(uVec, q));
    return rcpp_result_gen;
END_RCPP
}
// kern_fcn_MW210
arma::vec kern_fcn_MW210(arma::vec& uVec, double q);
RcppExport SEXP _DCSmooth_kern_fcn_MW210(SEXP uVecSEXP, SEXP qSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec& >::type uVec(uVecSEXP);
    Rcpp::traits::input_parameter< double >::type q(qSEXP);
    rcpp_result_gen = Rcpp::wrap(kern_fcn_MW210(uVec, q));
    return rcpp_result_gen;
END_RCPP
}
// kern_fcn_MW220
arma::vec kern_fcn_MW220(arma::vec& uVec, double q);
RcppExport SEXP _DCSmooth_kern_fcn_MW220(SEXP uVecSEXP, SEXP qSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec& >::type uVec(uVecSEXP);
    Rcpp::traits::input_parameter< double >::type q(qSEXP);
    rcpp_result_gen = Rcpp::wrap(kern_fcn_MW220(uVec, q));
    return rcpp_result_gen;
END_RCPP
}
// kern_fcn_MW320
arma::vec kern_fcn_MW320(arma::vec& uVec, double q);
RcppExport SEXP _DCSmooth_kern_fcn_MW320(SEXP uVecSEXP, SEXP qSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec& >::type uVec(uVecSEXP);
    Rcpp::traits::input_parameter< double >::type q(qSEXP);
    rcpp_result_gen = Rcpp::wrap(kern_fcn_MW320(uVec, q));
    return rcpp_result_gen;
END_RCPP
}
// kern_fcn_MW321
arma::vec kern_fcn_MW321(arma::vec& uVec, double q);
RcppExport SEXP _DCSmooth_kern_fcn_MW321(SEXP uVecSEXP, SEXP qSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec& >::type uVec(uVecSEXP);
    Rcpp::traits::input_parameter< double >::type q(qSEXP);
    rcpp_result_gen = Rcpp::wrap(kern_fcn_MW321(uVec, q));
    return rcpp_result_gen;
END_RCPP
}
// kern_fcn_MW420
arma::vec kern_fcn_MW420(arma::vec& uVec, double q);
RcppExport SEXP _DCSmooth_kern_fcn_MW420(SEXP uVecSEXP, SEXP qSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec& >::type uVec(uVecSEXP);
    Rcpp::traits::input_parameter< double >::type q(qSEXP);
    rcpp_result_gen = Rcpp::wrap(kern_fcn_MW420(uVec, q));
    return rcpp_result_gen;
END_RCPP
}
// kern_fcn_MW421
arma::vec kern_fcn_MW421(arma::vec& uVec, double q);
RcppExport SEXP _DCSmooth_kern_fcn_MW421(SEXP uVecSEXP, SEXP qSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec& >::type uVec(uVecSEXP);
    Rcpp::traits::input_parameter< double >::type q(qSEXP);
    rcpp_result_gen = Rcpp::wrap(kern_fcn_MW421(uVec, q));
    return rcpp_result_gen;
END_RCPP
}
// kern_fcn_MW422
arma::vec kern_fcn_MW422(arma::vec& uVec, double q);
RcppExport SEXP _DCSmooth_kern_fcn_MW422(SEXP uVecSEXP, SEXP qSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec& >::type uVec(uVecSEXP);
    Rcpp::traits::input_parameter< double >::type q(qSEXP);
    rcpp_result_gen = Rcpp::wrap(kern_fcn_MW422(uVec, q));
    return rcpp_result_gen;
END_RCPP
}
// kern_fcn_MW533
arma::vec kern_fcn_MW533(arma::vec& uVec, double q);
RcppExport SEXP _DCSmooth_kern_fcn_MW533(SEXP uVecSEXP, SEXP qSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec& >::type uVec(uVecSEXP);
    Rcpp::traits::input_parameter< double >::type q(qSEXP);
    rcpp_result_gen = Rcpp::wrap(kern_fcn_MW533(uVec, q));
    return rcpp_result_gen;
END_RCPP
}
// kern_fcn_MW644
arma::vec kern_fcn_MW644(arma::vec& uVec, double q);
RcppExport SEXP _DCSmooth_kern_fcn_MW644(SEXP uVecSEXP, SEXP qSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec& >::type uVec(uVecSEXP);
    Rcpp::traits::input_parameter< double >::type q(qSEXP);
    rcpp_result_gen = Rcpp::wrap(kern_fcn_MW644(uVec, q));
    return rcpp_result_gen;
END_RCPP
}
// kern_fcn_T220
arma::vec kern_fcn_T220(arma::vec& uVec, double q);
RcppExport SEXP _DCSmooth_kern_fcn_T220(SEXP uVecSEXP, SEXP qSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec& >::type uVec(uVecSEXP);
    Rcpp::traits::input_parameter< double >::type q(qSEXP);
    rcpp_result_gen = Rcpp::wrap(kern_fcn_T220(uVec, q));
    return rcpp_result_gen;
END_RCPP
}
// kern_fcn_T321
arma::vec kern_fcn_T321(arma::vec& uVec, double q);
RcppExport SEXP _DCSmooth_kern_fcn_T321(SEXP uVecSEXP, SEXP qSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec& >::type uVec(uVecSEXP);
    Rcpp::traits::input_parameter< double >::type q(qSEXP);
    rcpp_result_gen = Rcpp::wrap(kern_fcn_T321(uVec, q));
    return rcpp_result_gen;
END_RCPP
}
// kern_fcn_T420
arma::vec kern_fcn_T420(arma::vec& uVec, double q);
RcppExport SEXP _DCSmooth_kern_fcn_T420(SEXP uVecSEXP, SEXP qSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec& >::type uVec(uVecSEXP);
    Rcpp::traits::input_parameter< double >::type q(qSEXP);
    rcpp_result_gen = Rcpp::wrap(kern_fcn_T420(uVec, q));
    return rcpp_result_gen;
END_RCPP
}
// kern_fcn_T422
arma::vec kern_fcn_T422(arma::vec& uVec, double q);
RcppExport SEXP _DCSmooth_kern_fcn_T422(SEXP uVecSEXP, SEXP qSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec& >::type uVec(uVecSEXP);
    Rcpp::traits::input_parameter< double >::type q(qSEXP);
    rcpp_result_gen = Rcpp::wrap(kern_fcn_T422(uVec, q));
    return rcpp_result_gen;
END_RCPP
}
// LPSmooth_matrix2_BMod
arma::mat LPSmooth_matrix2_BMod(const arma::mat yMat, const double h, const int polyOrder, const int drv, const int mu, SEXP weightFcnPtr);
RcppExport SEXP _DCSmooth_LPSmooth_matrix2_BMod(SEXP yMatSEXP, SEXP hSEXP, SEXP polyOrderSEXP, SEXP drvSEXP, SEXP muSEXP, SEXP weightFcnPtrSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat >::type yMat(yMatSEXP);
    Rcpp::traits::input_parameter< const double >::type h(hSEXP);
    Rcpp::traits::input_parameter< const int >::type polyOrder(polyOrderSEXP);
    Rcpp::traits::input_parameter< const int >::type drv(drvSEXP);
    Rcpp::traits::input_parameter< const int >::type mu(muSEXP);
    Rcpp::traits::input_parameter< SEXP >::type weightFcnPtr(weightFcnPtrSEXP);
    rcpp_result_gen = Rcpp::wrap(LPSmooth_matrix2_BMod(yMat, h, polyOrder, drv, mu, weightFcnPtr));
    return rcpp_result_gen;
END_RCPP
}
// LP_dcs_const0_BMod
arma::mat LP_dcs_const0_BMod(arma::mat yMat, arma::colvec hVec, arma::icolvec polyOrderVec, arma::icolvec drvVec, arma::icolvec muVec, SEXP weightFcnPtr_x, SEXP weightFcnPtr_t);
RcppExport SEXP _DCSmooth_LP_dcs_const0_BMod(SEXP yMatSEXP, SEXP hVecSEXP, SEXP polyOrderVecSEXP, SEXP drvVecSEXP, SEXP muVecSEXP, SEXP weightFcnPtr_xSEXP, SEXP weightFcnPtr_tSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type yMat(yMatSEXP);
    Rcpp::traits::input_parameter< arma::colvec >::type hVec(hVecSEXP);
    Rcpp::traits::input_parameter< arma::icolvec >::type polyOrderVec(polyOrderVecSEXP);
    Rcpp::traits::input_parameter< arma::icolvec >::type drvVec(drvVecSEXP);
    Rcpp::traits::input_parameter< arma::icolvec >::type muVec(muVecSEXP);
    Rcpp::traits::input_parameter< SEXP >::type weightFcnPtr_x(weightFcnPtr_xSEXP);
    Rcpp::traits::input_parameter< SEXP >::type weightFcnPtr_t(weightFcnPtr_tSEXP);
    rcpp_result_gen = Rcpp::wrap(LP_dcs_const0_BMod(yMat, hVec, polyOrderVec, drvVec, muVec, weightFcnPtr_x, weightFcnPtr_t));
    return rcpp_result_gen;
END_RCPP
}
// LPSmooth_matrix_BMod
arma::mat LPSmooth_matrix_BMod(const arma::mat yMat, const double h, const int polyOrder, const int drv, const int mu, SEXP weightFcnPtr);
RcppExport SEXP _DCSmooth_LPSmooth_matrix_BMod(SEXP yMatSEXP, SEXP hSEXP, SEXP polyOrderSEXP, SEXP drvSEXP, SEXP muSEXP, SEXP weightFcnPtrSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat >::type yMat(yMatSEXP);
    Rcpp::traits::input_parameter< const double >::type h(hSEXP);
    Rcpp::traits::input_parameter< const int >::type polyOrder(polyOrderSEXP);
    Rcpp::traits::input_parameter< const int >::type drv(drvSEXP);
    Rcpp::traits::input_parameter< const int >::type mu(muSEXP);
    Rcpp::traits::input_parameter< SEXP >::type weightFcnPtr(weightFcnPtrSEXP);
    rcpp_result_gen = Rcpp::wrap(LPSmooth_matrix_BMod(yMat, h, polyOrder, drv, mu, weightFcnPtr));
    return rcpp_result_gen;
END_RCPP
}
// LP_dcs_const1_BMod
arma::mat LP_dcs_const1_BMod(arma::mat yMat, arma::colvec hVec, arma::icolvec polyOrderVec, arma::icolvec drvVec, arma::icolvec muVec, SEXP weightFcnPtr_x, SEXP weightFcnPtr_t);
RcppExport SEXP _DCSmooth_LP_dcs_const1_BMod(SEXP yMatSEXP, SEXP hVecSEXP, SEXP polyOrderVecSEXP, SEXP drvVecSEXP, SEXP muVecSEXP, SEXP weightFcnPtr_xSEXP, SEXP weightFcnPtr_tSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type yMat(yMatSEXP);
    Rcpp::traits::input_parameter< arma::colvec >::type hVec(hVecSEXP);
    Rcpp::traits::input_parameter< arma::icolvec >::type polyOrderVec(polyOrderVecSEXP);
    Rcpp::traits::input_parameter< arma::icolvec >::type drvVec(drvVecSEXP);
    Rcpp::traits::input_parameter< arma::icolvec >::type muVec(muVecSEXP);
    Rcpp::traits::input_parameter< SEXP >::type weightFcnPtr_x(weightFcnPtr_xSEXP);
    Rcpp::traits::input_parameter< SEXP >::type weightFcnPtr_t(weightFcnPtr_tSEXP);
    rcpp_result_gen = Rcpp::wrap(LP_dcs_const1_BMod(yMat, hVec, polyOrderVec, drvVec, muVec, weightFcnPtr_x, weightFcnPtr_t));
    return rcpp_result_gen;
END_RCPP
}
// ar_coef
arma::vec ar_coef(const arma::vec ar, const arma::vec ma, const double d, const int k);
RcppExport SEXP _DCSmooth_ar_coef(SEXP arSEXP, SEXP maSEXP, SEXP dSEXP, SEXP kSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec >::type ar(arSEXP);
    Rcpp::traits::input_parameter< const arma::vec >::type ma(maSEXP);
    Rcpp::traits::input_parameter< const double >::type d(dSEXP);
    Rcpp::traits::input_parameter< const int >::type k(kSEXP);
    rcpp_result_gen = Rcpp::wrap(ar_coef(ar, ma, d, k));
    return rcpp_result_gen;
END_RCPP
}
// ARMA_to_AR
NumericVector ARMA_to_AR(const arma::vec phi, const arma::vec psi, const int K);
RcppExport SEXP _DCSmooth_ARMA_to_AR(SEXP phiSEXP, SEXP psiSEXP, SEXP KSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec >::type phi(phiSEXP);
    Rcpp::traits::input_parameter< const arma::vec >::type psi(psiSEXP);
    Rcpp::traits::input_parameter< const int >::type K(KSEXP);
    rcpp_result_gen = Rcpp::wrap(ARMA_to_AR(phi, psi, K));
    return rcpp_result_gen;
END_RCPP
}
// sarma_rss
double sarma_rss(const arma::vec theta, const arma::mat R_mat, const List model_order);
RcppExport SEXP _DCSmooth_sarma_rss(SEXP thetaSEXP, SEXP R_matSEXP, SEXP model_orderSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< const arma::mat >::type R_mat(R_matSEXP);
    Rcpp::traits::input_parameter< const List >::type model_order(model_orderSEXP);
    rcpp_result_gen = Rcpp::wrap(sarma_rss(theta, R_mat, model_order));
    return rcpp_result_gen;
END_RCPP
}
// sfarima_rss
double sfarima_rss(const arma::vec theta, const arma::mat R_mat, const List model_order);
RcppExport SEXP _DCSmooth_sfarima_rss(SEXP thetaSEXP, SEXP R_matSEXP, SEXP model_orderSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< const arma::mat >::type R_mat(R_matSEXP);
    Rcpp::traits::input_parameter< const List >::type model_order(model_orderSEXP);
    rcpp_result_gen = Rcpp::wrap(sfarima_rss(theta, R_mat, model_order));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_DCSmooth_xMatrix", (DL_FUNC) &_DCSmooth_xMatrix, 2},
    {"_DCSmooth_cumsum_part_reverse", (DL_FUNC) &_DCSmooth_cumsum_part_reverse, 2},
    {"_DCSmooth_specDens_cpp", (DL_FUNC) &_DCSmooth_specDens_cpp, 2},
    {"_DCSmooth_acfMatrix_quarter2", (DL_FUNC) &_DCSmooth_acfMatrix_quarter2, 1},
    {"_DCSmooth_acfMatrix_quarter", (DL_FUNC) &_DCSmooth_acfMatrix_quarter, 1},
    {"_DCSmooth_acfMatrix_cpp", (DL_FUNC) &_DCSmooth_acfMatrix_cpp, 1},
    {"_DCSmooth_globalBndw_cpp", (DL_FUNC) &_DCSmooth_globalBndw_cpp, 2},
    {"_DCSmooth_specInt00_cpp", (DL_FUNC) &_DCSmooth_specInt00_cpp, 1},
    {"_DCSmooth_specIntDrv_cpp", (DL_FUNC) &_DCSmooth_specIntDrv_cpp, 4},
    {"_DCSmooth_localBndw_cpp", (DL_FUNC) &_DCSmooth_localBndw_cpp, 3},
    {"_DCSmooth_specDensEst_cpp", (DL_FUNC) &_DCSmooth_specDensEst_cpp, 4},
    {"_DCSmooth_fourier_mat", (DL_FUNC) &_DCSmooth_fourier_mat, 3},
    {"_DCSmooth_bartlett_mat", (DL_FUNC) &_DCSmooth_bartlett_mat, 5},
    {"_DCSmooth_KRSmooth_matrix2", (DL_FUNC) &_DCSmooth_KRSmooth_matrix2, 4},
    {"_DCSmooth_KR_dcs_const0", (DL_FUNC) &_DCSmooth_KR_dcs_const0, 5},
    {"_DCSmooth_KRSmooth_matrix", (DL_FUNC) &_DCSmooth_KRSmooth_matrix, 4},
    {"_DCSmooth_KR_dcs_const1", (DL_FUNC) &_DCSmooth_KR_dcs_const1, 5},
    {"_DCSmooth_kernel_fcn_assign", (DL_FUNC) &_DCSmooth_kernel_fcn_assign, 1},
    {"_DCSmooth_kernel_fcn_use", (DL_FUNC) &_DCSmooth_kernel_fcn_use, 3},
    {"_DCSmooth_weight_fcn_assign", (DL_FUNC) &_DCSmooth_weight_fcn_assign, 1},
    {"_DCSmooth_kern_fcn_M200", (DL_FUNC) &_DCSmooth_kern_fcn_M200, 2},
    {"_DCSmooth_kern_fcn_M210", (DL_FUNC) &_DCSmooth_kern_fcn_M210, 2},
    {"_DCSmooth_kern_fcn_M220", (DL_FUNC) &_DCSmooth_kern_fcn_M220, 2},
    {"_DCSmooth_kern_fcn_M321", (DL_FUNC) &_DCSmooth_kern_fcn_M321, 2},
    {"_DCSmooth_kern_fcn_M420", (DL_FUNC) &_DCSmooth_kern_fcn_M420, 2},
    {"_DCSmooth_kern_fcn_M421", (DL_FUNC) &_DCSmooth_kern_fcn_M421, 2},
    {"_DCSmooth_kern_fcn_M422", (DL_FUNC) &_DCSmooth_kern_fcn_M422, 2},
    {"_DCSmooth_kern_fcn_M533", (DL_FUNC) &_DCSmooth_kern_fcn_M533, 2},
    {"_DCSmooth_kern_fcn_M644", (DL_FUNC) &_DCSmooth_kern_fcn_M644, 2},
    {"_DCSmooth_kern_fcn_MW200", (DL_FUNC) &_DCSmooth_kern_fcn_MW200, 2},
    {"_DCSmooth_kern_fcn_MW210", (DL_FUNC) &_DCSmooth_kern_fcn_MW210, 2},
    {"_DCSmooth_kern_fcn_MW220", (DL_FUNC) &_DCSmooth_kern_fcn_MW220, 2},
    {"_DCSmooth_kern_fcn_MW320", (DL_FUNC) &_DCSmooth_kern_fcn_MW320, 2},
    {"_DCSmooth_kern_fcn_MW321", (DL_FUNC) &_DCSmooth_kern_fcn_MW321, 2},
    {"_DCSmooth_kern_fcn_MW420", (DL_FUNC) &_DCSmooth_kern_fcn_MW420, 2},
    {"_DCSmooth_kern_fcn_MW421", (DL_FUNC) &_DCSmooth_kern_fcn_MW421, 2},
    {"_DCSmooth_kern_fcn_MW422", (DL_FUNC) &_DCSmooth_kern_fcn_MW422, 2},
    {"_DCSmooth_kern_fcn_MW533", (DL_FUNC) &_DCSmooth_kern_fcn_MW533, 2},
    {"_DCSmooth_kern_fcn_MW644", (DL_FUNC) &_DCSmooth_kern_fcn_MW644, 2},
    {"_DCSmooth_kern_fcn_T220", (DL_FUNC) &_DCSmooth_kern_fcn_T220, 2},
    {"_DCSmooth_kern_fcn_T321", (DL_FUNC) &_DCSmooth_kern_fcn_T321, 2},
    {"_DCSmooth_kern_fcn_T420", (DL_FUNC) &_DCSmooth_kern_fcn_T420, 2},
    {"_DCSmooth_kern_fcn_T422", (DL_FUNC) &_DCSmooth_kern_fcn_T422, 2},
    {"_DCSmooth_LPSmooth_matrix2_BMod", (DL_FUNC) &_DCSmooth_LPSmooth_matrix2_BMod, 6},
    {"_DCSmooth_LP_dcs_const0_BMod", (DL_FUNC) &_DCSmooth_LP_dcs_const0_BMod, 7},
    {"_DCSmooth_LPSmooth_matrix_BMod", (DL_FUNC) &_DCSmooth_LPSmooth_matrix_BMod, 6},
    {"_DCSmooth_LP_dcs_const1_BMod", (DL_FUNC) &_DCSmooth_LP_dcs_const1_BMod, 7},
    {"_DCSmooth_ar_coef", (DL_FUNC) &_DCSmooth_ar_coef, 4},
    {"_DCSmooth_ARMA_to_AR", (DL_FUNC) &_DCSmooth_ARMA_to_AR, 3},
    {"_DCSmooth_sarma_rss", (DL_FUNC) &_DCSmooth_sarma_rss, 3},
    {"_DCSmooth_sfarima_rss", (DL_FUNC) &_DCSmooth_sfarima_rss, 3},
    {NULL, NULL, 0}
};

RcppExport void R_init_DCSmooth(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
