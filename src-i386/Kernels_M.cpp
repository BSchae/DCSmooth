// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

// Müller kernels
// Form is: Type_k_mu_nu

// [[Rcpp::export]]
arma::vec kern_fcn_M200(arma::vec& uVec, double q = 1)
{
  arma::uword n_size{ uVec.size() };
  arma::vec uOut(n_size);
  
  double q2{ q * q };
  double q_denom{ pow(1 + q, 3) };
  double out{ };
  double u{ };
  
  for (int i{ 0 }; i < n_size; ++i)
  {
    u   = uVec(i);
    
    out = 2 / q_denom * (
          (2*q2 - 2*q + 2) +
          (3 - 3*q) * u);
      
      uOut(i) = out;
  }
  return uOut;
}

// [[Rcpp::export]]
arma::vec kern_fcn_M210(arma::vec& uVec, double q = 1)
{
  arma::uword n_size{ uVec.size() };
  arma::vec uOut(n_size);
  
  double q2{ q * q };
  double q_denom{ pow(1 + q, 5) };
  double out{ };
  double u{ };
  
  for (int i{ 0 }; i < n_size; ++i)
  {
    u   = uVec(i);
    
    out = 12 * (1 + u) * (q - u) / q_denom * (
          (3*q2 - 4*q + 3) +
          (- 5*q + 5) * u);
    
    uOut(i) = out;
  }
  return uOut;
}

// [[Rcpp::export]]
arma::vec kern_fcn_M220(arma::vec& uVec, double q = 1)
{
  arma::uword n_size{ uVec.size() };
  arma::vec uOut(n_size);
  
  double q2{ q * q };
  double q_denom{ pow(1 + q, 7) };
  double out{ };
  double u{ };
  
  for (int i{ 0 }; i < n_size; ++i)
  {
    u   = uVec(i);
    
    out = 60 * (1 + u) * (1 + u) * (q - u) * (q - u) / q_denom * (
          (4*q2 - 6*q + 4) +
          (- 7*q + 7) * u);
    
    uOut(i) = out;
  }
  return uOut;
}

// [[Rcpp::export]]
arma::vec kern_fcn_M321(arma::vec& uVec, double q = 1)
{
  arma::uword n_size{ uVec.size() };
  arma::vec uOut(n_size);
  
  double q2{ q * q };
  double q3{ q2 * q };
  double q_denom{ pow(1 + q, 9) };
  double out{ };
  double u{ };
  double u2{ };
  
  for (int i{ 0 }; i < n_size; ++i)
  {
    u   = uVec(i);
    u2  = u * u;
    
    out = - 840 * (1 + u) * (1 + u) * (q - u) * (q - u) / q_denom * (
          (- 5*q3 + 16*q2 - 16*q + 5) +
          (22*q2 - 40*q + 22) * u +
          ( - 21*q + 21) * u2);
    
    uOut(i) = out;
  }
  return uOut;
}

// [[Rcpp::export]]
arma::vec kern_fcn_M420(arma::vec& uVec, double q = 1)
{
  arma::uword n_size{ uVec.size() };
  arma::vec uOut(n_size);
  
  double q2{ q * q };
  double q3{ q2 * q };
  double q4{ q3 * q };
  double q5{ q4 * q };
  double q6{ q5 * q };
  double q_denom{ pow(1 + q, 11) };
  double out{ };
  double u{ };
  double u2{ };
  double u3{ };
  
  for (int i{ 0 }; i < n_size; ++i)
  {
    u   = uVec(i);
    u2  = u * u;
    u3  = u2 * u;
    
    out = 840 * (1 + u) * (1 + u) * (q - u) * (q - u) / q_denom * (
          (4*q6 - 30*q5 + 96*q4 - 136*q3 + 96*q2 - 30*q + 4) +
          (- 27*q5 + 171*q4 - 396*q3 + 396*q2 - 171*q + 27) * u + 
          (54*q4 - 300*q3 + 480*q2 - 300*q + 54) * u2 +
          (- 33*q3 + 165*q2 - 165*q + 33) * u3);
    
    uOut(i) = out;
  }
  return uOut;
}

// [[Rcpp::export]]
arma::vec kern_fcn_M421(arma::vec& uVec, double q = 1)
{
  arma::uword n_size{ uVec.size() };
  arma::vec uOut(n_size);
  
  double q2{ q * q };
  double q3{ q2 * q };
  double q4{ q3 * q };
  double q5{ q4 * q };
  double q_denom{ pow(1 + q, 11) };
  double out{ };
  double u{ };
  double u2{ };
  double u3{ };
  
  for (int i{ 0 }; i < n_size; ++i)
  {
    u   = uVec(i);
    u2  = u * u;
    u3  = u2 * u;
    
    out = -2520 * (1 + u) * (1 + u) * (q - u) * (q - u) / q_denom * (
          (- 9*q5 + 57*q4 - 132*q3 + 132*q2 - 57*q + 9) +
          (66*q4 - 292*q3 + 472*q2 - 292*q + 66) * u +
          (- 139*q3 + 455*q2 - 455*q + 139) * u2 +
          (88*q2 - 220*q + 88) * u3);
    
    uOut(i) = out;
  }
  return uOut;
}

// [[Rcpp::export]]
arma::vec kern_fcn_M422(arma::vec& uVec, double q = 1)
{
  arma::uword n_size{ uVec.size() };
  arma::vec uOut(n_size);
  
  double q2{ q * q };
  double q3{ q2 * q };
  double q4{ q3 * q };
  double q_denom{ pow(1 + q, 11) };
  double out{ };
  double u{ };
  double u2{ };
  double u3{ };
  
  for (int i{ 0 }; i < n_size; ++i)
  {
    u   = uVec(i);
    u2  = u * u;
    u3  = u2 * u;
    
    out = 5040 * (1 + u) * (1 + u) * (q - u) * (q - u) / q_denom * (
          (18*q4 - 100*q3 + 160*q2 - 100*q + 18) +
          (- 139*q3 + 455*q2 - 455*q + 139) * u +
          (304*q2 - 580*q + 304) * u2 +
          (- 198*q + 198) * u3);
    
    uOut(i) = out;
  }
  return uOut;
}