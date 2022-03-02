// // [[Rcpp::depends(RcppArmadillo)]]
// #include <RcppArmadillo.h>
// // #include <cmath>
// #include "DCSmooth.h"
// #include "DCSmooth_types.h"
// 
// using namespace Rcpp;
// 
// arma::mat bartlett_mat(int Lx, int Lt, arma::vec drv1,
//                         arma::vec drv2, int power);
// double specInt00_cpp2(const arma::mat& acfMat);
// double specIntDrv_cpp2(const arma::mat& acfMat, arma::uvec drv1, 
//                          arma::uvec drv2, arma::uvec hLag);
// arma::cx_mat fourier_mat2(int Lx, int Lt, arma::vec omega);
// double specDensEst_cpp2(const arma::mat& acfMat, arma::vec drv, arma::uvec hLag,
//                         arma::vec omega);
// arma::uvec localBndw_cpp2(const arma::mat& acfMat, arma::uvec hLag,
//                         arma::vec omega);
// arma::vec globalBndw_cpp2(const arma::mat& acfMat, arma::uvec hLag, 
//                           arma::uword nX, arma::uword nT);
// 
// //-------------------------------Main Function--------------------------------//
// 
// //[[Rcpp::export]]
// arma::mat specDens_cpp2(const arma::mat& y0_mat, arma::vec omega)
// {
//   arma::uword nX{ y0_mat.n_rows }; // equals origin of acfMatrix
//   arma::uword nT{ y0_mat.n_cols };
//   
//   arma::mat y1_mat{ arma::reverse(y0_mat, 0) };
//   arma::mat acf_mat(2*nX - 1, 2*nT - 1, arma::fill::zeros);
//   
//   arma::vec hVec(2);
//   arma::uvec hVecInfl(2);
//   arma::vec hVecTemp(2);
//   arma::vec inflFct(2);
//   arma::mat y0;
//   arma::mat y1;
//   
//   // initial bandwidth values
//   hVec(0) = std::trunc(nX/2); hVec(1) = std::trunc(nT/2);
//   inflFct(0) = 1/pow(nX, 2.0/12.0); inflFct(1) = 1/pow(nT, 2.0/12.0);
//   arma::uword start_x{ 0 }; arma::uword start_t{ 0 };
//   arma::uword end_x{ 0 }; arma::uword end_t{ 0 };
//   
//   // global step
//   for (arma::uword g{ 0 }; g < 6; g++)
//   {
//     hVecTemp = hVec;
//     hVecInfl = arma::conv_to<arma::uvec>::from(hVec % inflFct) + 1;
//     
//     // calculate acfMatrix
//     if (hVecInfl(0) > end_x || hVecInfl(1) > end_t)
//     {
//       // calculate new max lags
//       start_x = end_x, start_t = end_t;
//       end_x = std::min(hVecInfl(0), nX - 1);
//       end_t = std::min(hVecInfl(1), nT - 1);
//       
//       for (arma::uword j{ start_x }; j < end_x + 1; j++)
//       {
//         for (arma::uword i{ start_t }; i < end_t + 1; i++)
//         {
//           // diagonal elements
//           y0 = y0_mat.submat(0, 0, nX - i - 1, nT - j - 1);
//           y1 = y0_mat.submat(i, j, nX - 1, nT - 1);
//           acf_mat.at(nX + i - 1, nT + j - 1) = arma::accu(y0 % y1)/(nX * nT);
//           acf_mat.at(nX - i - 1, nT - j - 1) =
//                                           acf_mat.at(nX + i - 1, nT + j - 1);
//           // off-diagonal elements
//           y0 = y1_mat.submat(0, 0, nX - i - 1, nT - j - 1);
//           y1 = y1_mat.submat(i, j, nX - 1, nT - 1);
//           acf_mat.at(nX - i - 1, nT + j - 1) = arma::accu(y0 % y1)/(nX * nT);
//           acf_mat.at(nX + i - 1, nT - j - 1) =
//             acf_mat.at(nX - i - 1, nT + j - 1);
//         }
//       }
//     }
//     
//     return(acf_mat);
//     
//     hVec = globalBndw_cpp2(acf_mat.submat(nX - hVecInfl(0), nT - hVecInfl(1),
//                                          nT + hVecInfl(0), nT + hVecInfl(1)),
//                                          hVecInfl, nX, nT);
//     
//     Rcout << arma::accu(acf_mat.submat(nX - hVecInfl(0), nT - hVecInfl(1),
//                          nT + hVecInfl(0), nT + hVecInfl(1))) << "---------\n";
//     
//     if (hVec(0) == hVecTemp(0) && hVec(1) == hVecTemp(1))
//     {
//       break;
//     }
//   }
//   
//   // local step
//   hVecInfl = arma::conv_to<arma::uvec>::from(arma::conv_to<arma::vec>::
//     from(hVec) % inflFct + 1);
// 
//   arma::uvec hOpt{ localBndw_cpp2(acf_mat.submat(nX - hVecInfl(0), nT - hVecInfl(1),
//                                                nT + hVecInfl(0), nT + hVecInfl(1)),
//                                                 hVecInfl, omega) };
// 
//   arma::vec drv0 = { 0, 0 };
// 
//   double specDensOut{ specDensEst_cpp2(acf_mat.submat(nX - hVecInfl(0), nT - hVecInfl(1),
//                                                      nT + hVecInfl(0), nT + hVecInfl(1)),
//                                                      drv0, hOpt, omega) };
// 
//   arma::vec returnVec(3);
//   returnVec(0) = specDensOut;
//   returnVec.subvec(1, 2) = arma::conv_to<arma::vec>::from(hOpt);
// 
//   // return returnVec;
// }
// 
// //---------------------Estimation of Global Bandwidth-------------------------//
// 
// //[[Rcpp::export]]
// arma::vec globalBndw_cpp2(const arma::mat& acfMat, arma::uvec hLag, 
//                           arma::uword nX, arma::uword nT)
// {
//   // arma::uword nX{ acfMat.n_rows/2 + 1 };
//   // arma::uword nT{ acfMat.n_cols/2 + 1 };
//   double mu_2w { 4.0/9.0 };
// 
//   double F00;
//   double F10;
//   double F01;
//   double F11;
// 
//   double Fx;
//   double Ft;
// 
//   arma::uword Lx;
//   arma::uword Lt;
// 
//   arma::uvec drv10{ 1, 0 };
//   arma::uvec drv01{ 0, 1 };
// 
//   if (hLag(0) > 1 && hLag(1) > 1)
//   {
//     F00 = specInt00_cpp2(acfMat);
//     F10 = specIntDrv_cpp2(acfMat, drv10, drv10, hLag);
//     F01 = specIntDrv_cpp2(acfMat, drv01, drv01, hLag);
//     F11 = specIntDrv_cpp2(acfMat, drv10, drv01, hLag);
// 
//     Fx = F10 * (sqrt(F10/F01) + F11/F01) / F00;
//     Ft = F01 * (sqrt(F01/F10) + F11/F10) / F00;
// 
//     Lx = std::trunc(pow(2 * Fx * nX*nT/mu_2w, 0.25));
//     Lt = std::trunc(pow(2 * Ft * nX*nT/mu_2w, 0.25));
//   }
//   else if (hLag(0) == 1 && hLag(1) > 1)
//   {
//     F00 = specInt00_cpp2(acfMat);
//     F10 = specIntDrv_cpp2(acfMat, drv10, drv10, hLag);
// 
//     Fx = F10/F00;
// 
//     Lx = std::trunc(pow(2 * Fx * nX*nT/mu_2w, 1.0/3.0));
//     Lt = 0;
//   }
//   else if (hLag(0) > 1 && hLag(1) == 1)
//   {
//     F00 = specInt00_cpp2(acfMat);
//     F01 = specIntDrv_cpp2(acfMat, drv01, drv01, hLag);
// 
//     Ft =  F01/F00;
// 
//     Lx = 0;
//     Lt = std::trunc(pow(2 * Ft * nX*nT/mu_2w, 1.0/3.0));
//   }
//   else
//   {
//     Lx = 0;
//     Lt = 0;
//   }
// 
//   arma::vec hOut = { std::min(Lx + 1, nX - 1), std::min(Lt + 1, nT - 1) };
// 
//   return hOut;
// }
// 
// //[[Rcpp::export]]
// double specInt00_cpp2(const arma::mat& acfMat) //, int nX, int nT)
// {
//   arma::uvec dX{ acfMat.n_rows };
//   arma::uvec dT{ acfMat.n_cols };
//   
//   double out{ arma::accu(acfMat % acfMat) };
//   return pow(2 * M_PI, - 2) * out; // * 0.5;
// }
// 
// //[[Rcpp::export]]
// double specIntDrv_cpp2(const arma::mat& acfMat, arma::uvec drv1, arma::uvec drv2,
//                       arma::uvec hLag)
// {
//   arma::uword origin_x{ acfMat.n_rows/2 };  // as first element has index 0,
//   arma::uword origin_t{ acfMat.n_cols/2 };  // origin is at floor(n/2)
//   int Lx = static_cast<int>(hLag(0));
//   int Lt = static_cast<int>(hLag(1));
// 
//   arma::mat w_mat(2*(Lx - 1) + 1, 2*(Lt - 1) + 1);
//   
//   arma::vec drv1v{ arma::conv_to<arma::vec>::from(drv1) };
//   arma::vec drv2v{ arma::conv_to<arma::vec>::from(drv2) };
//   
// 
//   // compute matrix of bartlett weights without "fourier factor"
//   arma::mat weights = bartlett_mat(Lx, Lt, drv1v, drv2v, 2);
// 
//   if (Lx > 1 && Lt > 1)
//   {
//     w_mat.submat(0, 0, Lx - 1, Lt - 1) =
//                   reverse(reverse(weights.submat(0, 0, Lx - 1, Lt - 1), 1), 0);
//     w_mat.submat(Lx, 0, 2*(Lx - 1), Lt - 2) =
//                   reverse(weights.submat(1, 1, Lx - 1, Lt - 1), 1);
//     w_mat.submat(0, Lt, Lx - 2, 2*(Lt - 1)) =
//                   reverse(weights.submat(1, 1, Lx - 1, Lt - 1), 0);
//     w_mat.submat(Lx - 1, Lt - 1, 2*(Lx - 1), 2*(Lt - 1)) =
//                   weights.submat(0, 0, Lx - 1, Lt - 1);
//     // are the "reverse" necessary, as w_mat is probably symmetric?
//   } else if (Lx == 1 && Lt == 1) {
//     w_mat.submat(0, 0, 0, 0) = 1;
//   }
// 
//   arma::mat acf_submat{ acfMat.submat(origin_x - Lx + 1, origin_t - Lt + 1,
//                                       origin_x + Lx - 1, origin_t + Lt - 1) };
// 
//   double F_out = arma::accu(w_mat % acf_submat % acf_submat);
// 
//   return pow(2 * M_PI, - 2) * F_out;
// }
// 
// //----------------------Estimation of Local Bandwidth-------------------------//
// 
// // [[Rcpp::export]]
// arma::uvec localBndw_cpp2(const arma::mat& acfMat, arma::uvec hLag, arma::vec omega)
// {
//   int nX{ acfMat.n_rows/2 + 1 };
//   int nT{ acfMat.n_cols/2 + 1 };
//   double mu_2w { 4.0/9.0 };
// 
//   double f00;
//   double f10;
//   double f01;
//   double f11;
// 
//   double fx;
//   double ft;
// 
//   int Lx;
//   int Lt;
// 
//   arma::vec drv00{ 0, 0 };
//   arma::vec drv10{ 1, 0 };
//   arma::vec drv01{ 0, 1 };
// 
//   if (hLag(0) > 1 && hLag(1) > 1)
//   {
//     f00 = specDensEst_cpp2(acfMat, drv00, hLag, omega);
//     f10 = specDensEst_cpp2(acfMat, drv10, hLag, omega);
//     f01 = specDensEst_cpp2(acfMat, drv01, hLag, omega);
// 
//     Lx = static_cast<int>(pow(4 * std::abs(pow(f10, 3)/f01) * 1/pow(f00, 2) *
//       nX*nT/mu_2w, 0.25));
//     Lt = static_cast<int>(pow(4 * std::abs(pow(f01, 3)/f10) * 1/pow(f00, 2) *
//       nX*nT/mu_2w, 0.25));
//   }
//   else if (hLag(0) == 1 && hLag(1) > 1)
//   {
//     f00 = specDensEst_cpp2(acfMat, drv00, hLag, omega);
//     f10 = specDensEst_cpp2(acfMat, drv10, hLag, omega);
// 
//     Lx = static_cast<int>(pow(2 * pow(f10/f00, 2) *
//       nX*nT/sqrt(mu_2w), 1.0/3.0));
//     Lt = 0;
//   }
//   else if (hLag(0) > 1 && hLag(1) == 1)
//   {
//     f00 = specDensEst_cpp2(acfMat, drv00, hLag, omega);
//     f01 = specDensEst_cpp2(acfMat, drv01, hLag, omega);
// 
//     Lx = 0;
//     Lt = static_cast<int>(pow(2 * pow(f01/f00, 2) *
//       nX*nT/sqrt(mu_2w), 1.0/3.0)); // power 1/3 is correct?
//   }
//   else
//   {
//     Lx = 0;
//     Lt = 0;
//   }
// 
//   arma::uvec hOut = { std::min(Lx + 1, nX - 1), std::min(Lt + 1, nT - 1) };
// 
//   return hOut;
// }
// 
// //[[Rcpp::export]]
// double specDensEst_cpp2(const arma::mat& acfMat, arma::vec drv, arma::uvec hLag,
//                        arma::vec omega)
// {
//   arma::uword origin_x{ acfMat.n_rows/2 };  // as first element has index 0,
//   arma::uword origin_t{ acfMat.n_cols/2 };  // origin is at floor(n/2)
//   int Lx{ static_cast<int>(hLag(0)) };
//   int Lt{ static_cast<int>(hLag(1)) };
//   arma::vec drv0 = { 0, 0 };
// 
//   arma::cx_mat w_mat(2*(Lx - 1) + 1, 2*(Lt - 1) + 1);
// 
//   // compute matrix of bartlett weights including "fourier factor"
//   arma::cx_mat weights{ bartlett_mat(Lx, Lt, drv, drv0, 1) %
//     fourier_mat2(Lx, Lt, omega) };
// 
//   if (Lx > 1 && Lt > 1)
//   {
//     w_mat.submat(0, 0, Lx - 1, Lt - 1) =
//       reverse(reverse(weights.submat(0, 0, Lx - 1, Lt - 1), 1), 0);
//     w_mat.submat(Lx, 0, 2*(Lx - 1), Lt - 2) =
//       reverse(weights.submat(1, 1, Lx - 1, Lt - 1), 1);
//     w_mat.submat(0, Lt, Lx - 2, 2*(Lt - 1)) =
//       reverse(weights.submat(1, 1, Lx - 1, Lt - 1), 0);
//     w_mat.submat(Lx - 1, Lt - 1, 2*(Lx - 1), 2*(Lt - 1)) =
//       weights.submat(0, 0, Lx - 1, Lt - 1);
//     // are the "reverse" necessary, as w_mat is probably symmetric?
//   } else if (Lx == 1 && Lt == 1) {
//     w_mat.submat(0, 0, 0, 0) = 1;
//   }
// 
//   arma::mat acf_submat{ acfMat.submat(origin_x - Lx + 1, origin_t - Lt + 1,
//                                       origin_x + Lx - 1, origin_t + Lt - 1) };
// 
//   // calculate f. Note that f(-kX, - kT) = f(kX, kT), off diagonal sub matrices
//   // are different etc.
//   // note that w is (Lx + 1)x(Lt + 1) but the outer values are all zeroes,
//   // thus summation over them is not needed.
//   double fOut{ arma::accu(arma::real(w_mat) % acf_submat) };
// 
//   return pow(2 * M_PI, -2) * fOut;
// }
// 
// //---------------------------Additional Functions-----------------------------//
// 
// // [[Rcpp::export]]
// arma::cx_mat fourier_mat2(int Lx, int Lt, arma::vec omega)
// {
//   arma::vec arg_x{ -arma::regspace(0, Lx) * omega(0) };
//   arma::vec arg_t{ -arma::regspace(0, Lt) * omega(1) };
// 
//   arma::cx_mat mat_out{ arma::cx_vec(cos(arg_x), sin(arg_x)) *
//     arma::cx_vec(cos(arg_t), sin(arg_t)).st() };
// 
//   return mat_out;
// }
// // 
// // //[[Rcpp::export]]
// // arma::mat bartlett_mat2(int Lx, int Lt, arma::uvec drv1,
// //                        arma::vec drv2, double power)
// // {
// //   arma::mat lx_mat(Lx + 1, Lt + 1);
// //   arma::mat lt_mat(Lt + 1, Lx + 1); // needs to be transposed later
// //   arma::vec lx_col{ pow(arma::regspace(0, Lx), drv1(0) + drv2(0)) };
// //   arma::vec lt_row{ pow(arma::regspace(0, Lt), drv1(1) + drv2(1)) };
// // 
// //   for (arma::uword i{ 0 }; i <= Lt; i++)
// //   {
// //     lx_mat.col(i) = lx_col;
// //   }
// // 
// //   for (arma::uword j{ 0 }; j <= Lx; j++)
// //   {
// //     lt_mat.col(j) = lt_row; // cols are rows after transposing
// //   }
// // 
// //   arma::mat weights = std::pow((1 - arma::regspace(0, Lx)/std::max(Lx, 1)) *
// //     (1 - arma::regspace(0, Lt).t()/std::max(Lt, 1)), power) %
// //     lx_mat % lt_mat.t();
// // 
// //   return weights;
// // }