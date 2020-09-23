# ################################################################################
# #                                                                              #
# #                 DCSmooth Package: stimation of QARMA for cf                  #
# #                                                                              #
# ################################################################################
# 
# #-------------------Estimation of Autocovariance Function----------------------#
# 
# acfMatrix = function(Y, lagX = 100, lagT = 100)
# {
#   nx = dim(Y)[1]; nt = dim(Y)[2]
#   if(lagX > nx - 1) { lagX = nx - 1 }
#   if(lagT > nt - 1) { lagT = nt - 1 }
#   
#   Yk = Y[1:(lagX + 1), 1:(lagT + 1)] - mean(Y[(lagX + 1), 1:(lagT + 1)])
#   R = Yk*NA
#   Yk = as.matrix(Yk); R = as.matrix(R)
#   
#   for (kx in 0:lagX)
#   {
#     for (kt in 0:lagT)
#     {
#       R[kx + 1, kt + 1] = sum(Yk[1:(nx - kx), 1:(nt - kt)] *
#                                 Yk[(kx + 1):nx, (kt + 1):nt])
#     }
#   }
#   
#   return(R/(nx*nt))
# }
# 
# #----------------------Estimation of Spectral Density--------------------------#
# 
# specDensEst = function(acfMat, drv, h, omega)
# {
#   nx = dim(as.matrix(acfMat))[1]; nt = dim(as.matrix(acfMat))[2]
#   hX = h[1]; hT = h[2]; drvX = drv[1]; drvT = drv[2]
#   bndwX = max(trunc(hX * nx), 1); bndwT = max(trunc(hT * nt), 1)
#   # absolute bandwidths
#   
#   # compute matrix of bartlett weights multiplied with "fourier factor"
#   w = (1 - 0:(bndwX - 1)/(nx*hX)) %*% t(1 - 0:(bndwT - 1)/(nt*hT)) * 
#     complex(argument = -(0:(bndwX - 1))*omega[1]) %*% t(complex(argument =
#                                           - (0:(bndwT - 1))*omega[2])) *
#     matrix((0:(bndwX - 1))^drvX, bndwX, bndwT) *
#     matrix((0:(bndwT - 1))^drvT, bndwX, bndwT, byrow = TRUE)
#   
#   # calculate f. Note that f(-kX, - kT) = f(kX, kT) etc.
#   fOut = 4 * sum(w[1:bndwX, 1:bndwT] * acfMat[1:bndwX, 1:bndwT]) -
#     2 * ifelse(bndwT > 1, sum(w[1, 2:bndwT] * acfMat[1, 2:bndwT]), 0) -
#     2 * ifelse(bndwX > 1, sum(w[2:bndwX, 1] * acfMat[2:bndwX, 1]), 0) -
#     3 * acfMat[1, 1]
#   
#   return((2*pi)^(-2) * Re(fOut))
# }
# 
# #----------------------Estimation of Integral over SD--------------------------#
# 
# # Function for spectral density
# specIntEst00 = function(acfMat)
# {
#   nx = dim(as.matrix(acfMat))[1]; nt = dim(as.matrix(acfMat))[2]
#   
#   # calculate f. Note that f(-kX, - kT) = f(kX, kT) etc.
#   iOut = 4 * sum(acfMat[1:nx, 1:nt]^2) -
#     2 * ifelse(nt > 1, sum(acfMat[1, 2:nt]^2), 0) -
#     2 * ifelse(nx > 1, sum(acfMat[2:nx, 1]^2), 0) - 3 * acfMat[1, 1]^2
#   
#   return((2*pi)^(-2) * iOut)
# }
# 
# # Function for derivativesof Spectral Density
# specIntEstDrv = function(acfMat, drv1, drv2, h)
# {
#   nx = dim(as.matrix(acfMat))[1]; nt = dim(as.matrix(acfMat))[2]
#   hX = h[1]; hT = h[2];
#   bndwX = max(trunc(hX * nx), 1); bndwT = max(trunc(hT * nt), 1) # absolute bandwidths
#   
#   # compute matrix of bartlett weights multiplied with "fourier factor"
#   w = (1 - 0:(bndwX - 1)/(nx*hX)) %*% t(1 - 0:(bndwT - 1)/(nt*hT)) *
#     matrix((0:(bndwX - 1))^(drv1[1] + drv2[1]), bndwX, bndwT) *
#     matrix((0:(bndwT - 1))^(drv1[2] + drv2[2]), bndwX, bndwT, byrow = TRUE)
#   
#   # calculate F. Note that f(-kX, - kT) = f(kX, kT) etc.
#   FOut = 4 * sum(w[1:bndwX, 1:bndwT]^2 * acfMat[1:bndwX, 1:bndwT]^2) -
#     2 * ifelse(bndwT > 1, sum(w[1, 2:bndwT]^2 * acfMat[1, 2:bndwT]^2), 0) -
#     2 * ifelse(bndwX > 1, sum(w[2:bndwX, 1]^2 * acfMat[2:bndwX, 1]^2), 0) -
#     3 * w[1, 1]^2 * acfMat[1, 1]^2
# 
#   return((2*pi)^(-2) * FOut)
# }
# 
# #-----------------------Estimation of Bandwidths-------------------------------#
# 
# localBndwEst = function(acfMat, h, omega)
# {
#   nx = dim(as.matrix(acfMat))[1]; nt = dim(as.matrix(acfMat))[2]
#   
#   f00 = specDensEst(acfMat, drv = c(0, 0), h = h, omega = omega)
#   f10 = specDensEst(acfMat, drv = c(1, 0), h = h, omega = omega)
#   f01 = specDensEst(acfMat, drv = c(0, 1), h = h, omega = omega)
#   
#   mu_2w = 4/9
#   
#   numer = mu_2w * f01 * f00^2
#   denom = 4 * nx * nt * f10^3
#   hx = (numer/denom)^0.25
#   ht = hx * (f10/f01)
#   return(c(hx, ht))
# }
# 
# globalBndwEst = function(acfMat, hInfl)
# {
#   nx = dim(as.matrix(acfMat))[1]; nt = dim(as.matrix(acfMat))[2]
#   
#   F00 = specIntEst00(acfMat)
#   F10 = specIntEstDrv(acfMat, drv1 = c(1, 0), drv2 = c(1, 0), h = hInfl)
#   F01 = specIntEstDrv(acfMat, drv1 = c(0, 1), drv2 = c(0, 1), h = hInfl)
#   F1001 = specIntEstDrv(acfMat, drv1 = c(1, 0), drv2 = c(0, 1), h = hInfl)
#   
#   mu_2w = 4/9
# 
#   F = (F10/F01)^0.5 + F1001/F01
#   numer = mu_2w * F00
#   denom = 2 * nx * nt * F10 * F
#   hx = (numer/denom)^0.25
#   ht = hx * sqrt(F10/F01)
#   
#   return(c(hx, ht))
# }
# 
# #-----------------------Estimation of Variance Factor--------------------------#
# 
# cfW = function(Y, omega)
# {
#   nx = dim(as.matrix(Y))[1]; nt = dim(as.matrix(Y))[2]
#   acfMat = acfMatrix(as.matrix(Y))
#                      
#   # initial values
#   hVec = c(0.5, 0.5)
# 
#   
#   hTest = matrix(NA, 20, 2)
#   # global step
#   for (g in 1:20)
#   {
#     hVecInfl = hVec * c(nx^(2/21), nt^(2/21))
#     hVecInfl = pmin(pmax(hVecInfl, c(2/nx, 2/nt)), c(1, 1))
#     hVec = globalBndwEst(acfMat, hVecInfl)
#     hTest[g, ] = hVec
#   }
#   
#   # local step
#   hVecInfl = hVec * c(nx^(2/21), nt^(2/21))
#   hVec = localBndwEst(acfMat, hVecInfl, omega)
# 
#   specDensOut = specDensEst(acfMat, drv = c(0, 0), h = hVec, omega = omega)
#   
#   return(list(cf = specDensOut*(2*pi)^2, h = hVec))
# }
