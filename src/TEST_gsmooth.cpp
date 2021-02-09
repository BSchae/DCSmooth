# include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

// FUnctions written by Dominik Schulz, 12/05/2020

// Function to help obtain a sequence from 'from' to 'to' by step 1
// in C++
// [[Rcpp::export]]
arma::vec seqCpp(int from, int to){
  int n = to - from + 1;
  arma::vec seqOut(n);
  for (int i = from; i < to + 1; ++i) {
    seqOut(i - from) = i;
  }
  return seqOut;
}

// Function to calculate the factorial of an 'int' in C++
// [[Rcpp::export]]
int factorialCpp(int k) {
  int fac = 1;
  if (k > 1) {
    for (int i = 2; i < k + 1; ++i) {
      fac *= i;
    }
  }
  return fac;
}

// C++ version of gsmoothCalc
// [[Rcpp::export]]
arma::vec gsmoothCalcCpp(arma::vec y, int v, int p, int mu, double b, int bb) {
  int n = y.size();
  arma::vec gr(n);
  int hh = trunc(n * b + 0.5);
  int htm = 2 * hh + 1;
  arma::mat ws(htm, htm);
  arma::vec wk(htm);
  arma::mat xt((p + 1), htm);
  arma::mat xw = xt;
  
  arma::vec hr = hh + bb * (hh - seqCpp(0, hh));
  arma::vec ht = 1 + seqCpp(0, hh) + hr;
  arma::mat xa = xt;
  
  return hr;
  
  for (int i = 0; i < hh + 1; ++i) {
    wk.subvec(0, ht(i) - 1) = pow(1 - pow((seqCpp(1, ht(i)) - (i + 1)) / (hr(i) + 1), 2), mu);
    
    for (int j = 0; j < p + 1; ++j) {
      xt.submat(j, 0, j, ht(i) - 1) = (pow((seqCpp(1, ht(i)) - (i + 1)) / hh, j)).t();
      xw.row(j) = xt.row(j) % wk.t();
    }
    
    xa = inv(xw * xt.t()) * xw;
    
    ws.row(i) = xa.row(v);
  }
  
  ws.rows(hh + 1, htm - 1) = pow(-1, v) * flipud(fliplr(ws.submat(0, 0, hh - 1, htm - 1)));
  //ws = factorialCpp(v) * ws * pow(n / double(hh), v);
  
  return(ws);
  
  arma::mat ym(n - htm + 1, htm);
  for (int i = hh; i < n - hh; ++i) {
    ym.row(i - hh) = (y.subvec(i - hh, i + hh)).t();
  }
  gr.subvec(0, hh - 1) = ws.rows(0, hh - 1) * y.subvec(0, htm - 1);
  gr.subvec(hh, n - hh - 1) = ym * (ws.row(hh)).t();
  gr.subvec(n - hh, n - 1) = ws.rows(hh + 1, htm - 1) * y.subvec(n - htm, n - 1);
  
  return Rcpp::NumericVector(gr.begin(), gr.end());
}

// C++ version of gsmoothCalc2
// [[Rcpp::export]]
List gsmoothCalc2Cpp(arma::vec y, int v, int p, int mu, double b, int bb) {
  int n = y.size();
  arma::vec gr(n);
  int hh = trunc(n * b + 0.5);
  int htm = 2 * hh + 1;
  arma::mat ws(htm, htm);
  arma::vec wk(htm);
  arma::mat xt((p + 1), htm);
  arma::mat xw = xt;
  
  arma::vec hr = hh + bb * (hh - seqCpp(0, hh));
  arma::vec ht = 1 + seqCpp(0, hh) + hr;
  arma::mat xa = xt;
  
  for (int i = 0; i < hh + 1; ++i) {
    
    wk.subvec(0, ht(i) - 1) = pow(1 - pow((seqCpp(1, ht(i)) - (i + 1)) / (hr(i) + 1), 2), mu);
    
    for (int j = 0; j < p + 1; ++j) {
      xt.submat(j, 0, j, ht(i) - 1) = (pow((seqCpp(1, ht(i)) - (i + 1)) / hh, j)).t();
      xw.row(j) = xt.row(j) % wk.t();
    }
    
    xa = inv(xw * xt.t()) * xw;
    
    ws.row(i) = xa.row(v);
  }
  
  ws.rows(hh + 1, htm - 1) = pow(-1, v) * flipud(fliplr(ws.submat(0, 0, hh - 1, htm - 1)));
  ws = factorialCpp(v) * ws * pow(n / hh, v);
  arma::mat ym(n - htm + 1, htm);
  for (int i = hh; i < n - hh; ++i) {
    ym.row(i - hh) = (y.subvec(i - hh, i + hh)).t();
  }
  gr.subvec(0, hh - 1) = ws.rows(0, hh - 1) * y.subvec(0, htm - 1);
  gr.subvec(hh, n - hh - 1) = ym * (ws.row(hh)).t();
  gr.subvec(n - hh, n - 1) = ws.rows(hh + 1, htm - 1) * y.subvec(n - htm, n - 1);
  
  List listOut = List::create(_["ye"] = Rcpp::NumericVector(gr.begin(), gr.end()),
                              _["ws"] = ws);
  
  return listOut;
  
}

// Updated version of the iterative process / loop in tsmooth in C++ (not working yet)
// // [[Rcpp::export]]
// double tsmoothLoopCpp(arma::vec y, double c1, double c2, int p, int n, int n1, double bStart, int mu, int bb, String bvc, String Mcf) {
//   int k = p + 1;
//   int pd = p + 2;
//   double bold1;
//   double bold;
//   double bd;
//   double bopt;
//   double bv;
//   double I2;
//   arma::vec I2pre(n - 2 * n1);
//   arma::vec yed(n);
//   arma::vec ye(n);
//   arma::vec yd(n);
//   NumericVector yd2(n);
//   double c3;
//   List cf0est;
//   double cf0LW = NA_REAL;
//   double cf0 = NA_REAL;
//   double cf0AR = NA_REAL;
//   double cf0MA = NA_REAL;
//   double cf0ARMA = NA_REAL;
//   double L0opt = NA_REAL;
//   double ARMABIC = NA_REAL;
//   double ARBIC = NA_REAL;
//   double MABIC = NA_REAL;
//   int pBIC = NA_INTEGER;
//   int qBIC = NA_INTEGER;
//   int runc = 1;
//   arma::vec steps(40);
// //  Function cfAR("cf0.AR.est");
// //  Function cfMA("cf0.MA.est");
// //  Function cfARMA("cf0.ARMA.est");
//
//   for (int i = 0; i < 40; ++i) {
//     if (runc == 1) {
//       if (i > 0) {
//         bold1 = bold;
//       }
//       if (i == 0) {
//         bold = bStart;
//       } else {
//         bold = bopt;
//       }
//
//       // to be adjusted!!!
//       bd = pow(bold, 5 / 7.0);
//
//         if (bd >= 0.49) {
//           bd = 0.49;
//         }
//         yed = gsmoothCalcCpp(y, k, pd, mu, bd, bb);
//         // !!!
//         I2pre = pow(yed.subvec(n1, n - n1 - 1), 2.0);
//         I2 = sum(I2pre) / (n - 2 * n1);
//
// // Use an enlarged bandwidth for estimating cf or not (Feng/Heiler, 2009)
//           if (bvc == "Y") {
//
// // Look up the enlarged bandwidth
//             bv = 1.4310 * bold;
//
//           } else {
//             bv = bold;
//           }
//
//           if (bv >= 0.49) {
//             bv = 0.49;
//           }
//
//           ye = gsmoothCalcCpp(y, 0, p, mu, bv, bb);
//
// // The optimal bandwidth-----------------------------------------------
//
// // Estimating the variance factor
//           yd = y - ye;
//
//           yd2 = Rcpp::NumericVector(yd.begin(), yd.end());
// // Lag-Window for estimating the variance
//           if (Mcf == "NP") {
//
//             cf0est = cf0Cpp(yd2);
//             cf0LW = cf0est["cf0.LW"];
//             cf0 = cf0LW;
//             L0opt = cf0est["L0.opt"];
//
//           } else if (Mcf == "ARMA") {
// //            cf0est = cfARMA(_["Xt"] = yd2);
//             cf0ARMA = cf0est["cf0.ARMA"];
//             cf0 = cf0ARMA;
//             ARMABIC = cf0est["ARMA.BIC"];
//             pBIC = cf0est["p.BIC"];
//             qBIC = cf0est["q.BIC"];
//
//           } else if (Mcf == "AR") {
// //            cf0est = cfAR(_["Xt"] = yd2);
//             cf0AR = cf0est["cf0.AR"];
//             cf0 = cf0AR;
//             ARBIC = cf0est["AR.BIC"];
//             pBIC = cf0est["p.BIC"];
//
//           } else if (Mcf == "MA") {
// //            cf0est = cfMA(_["Xt"] = yd2);
//             cf0MA = cf0est["cf0.MA"];
//             cf0 = cf0MA;
//             MABIC = cf0est["MA.BIC"];
//             qBIC = cf0est["q.BIC"];
//
//           }
//
//           c3 = cf0 / I2;
//
//           if (p == 1) {
//             bopt = pow(c1 * c2 * c3, 1 / 5.0) * pow(n, -1 / 5.0);
//             if (bopt < pow(n, -5 / 7.0)) {
//               bopt = pow(n, -5 / 7.0);
//             }
//           }
//           if (p == 3) {
//             bopt = pow(c1 * c2 * c3, 1 / 9.0) * pow(n, -1 / 9.0);
//             if (bopt < pow(n, -9 / 11.0)) {
//               bopt = pow(n, -9 / 11.0);
//             }
//           }
//           if (bopt > 0.49) {
//             bopt = 0.49;
//           }
//
//           steps(i) = bopt;
//
//           if (i > 1 & abs(bold - bopt) / bopt < 1 / n) {
//             runc = 0;
//           }
//           if (i > 2 & abs(bold1 - bopt) / bopt < 1 / n) {
//             bopt = (bold + bopt) / 2;
//             runc = 0;
//           }
//     }
//   }
//   return bopt;
// }

// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically
// run after the compilation.
//

/*** R
# data <- smoots::dax
# close <- data$Close
# ret <- diff(log(close))
# y <- ret - mean(ret)
# X <- log(y^2)
# X <- log(smoots::gdpUS$GDP)
# est <- gsmoothCalcCpp(X, 2, 3, 1, 0.15, 1)
# est2 <- gsmoothCalc(X, bb = 1, p = 3, v = 2)
# est
# est2

# bOld <- bench::mark(estOld <- gsmoothCalc(X, bb = 1, p = 3))
# bNew <- bench::mark(estNew <- gsmoothCalcCpp(X, 0, 3, 1, 0.15, 1))

*/
