################################################################################
#                                                                              #
#           DCSmooth Package: Spectral Density Estimation via Buehlmann        #
#                                                                              #
################################################################################

#-------------------Estimation of Autocovariance Function----------------------#

acfMatrix = function(Y)
{
  Y = as.matrix(Y - mean(Y)) 
  nx = dim(Y)[1]; nt = dim(Y)[2]
  lagX = nx - 1; lagT = nt - 1
  # if(lagX > nx - 1) { lagX = nx - 1 }
  # if(lagT > nt - 1) { lagT = nt - 1 }
  
  Yk = matrix(Y[1:(lagX + 1), 1:(lagT + 1)], nrow = lagX + 1, ncol = lagT + 1)
  R = Yk*NA
  Yk = as.matrix(Yk); R = as.matrix(R)
  
  for (kx in 0:lagX)
  {
    for (kt in 0:lagT)
    {
      R[kx + 1, kt + 1] = sum(Yk[1:(lagX - kx), 1:(lagT - kt)] *
                                Yk[(kx + 1):lagX, (kt + 1):lagT])
    }
  }
  
  return(R/(nx*nt))
}

#----------------------Estimation of Spectral Density--------------------------#

specDensEst = function(acfMat, drv, hLag, omega)
{
  nx = dim(as.matrix(acfMat))[1]; nt = dim(as.matrix(acfMat))[2]
  Lx = hLag[1]; Lt = hLag[2]; drvX = drv[1]; drvT = drv[2]
  
  # absolute bandwidths
  
  # compute matrix of bartlett weights multiplied with "fourier factor"
  w = (1 - (0:Lx)/max(Lx, 1)) %*% t(1 - (0:Lt)/max(Lt, 1)) * 
    complex(argument = -(0:Lx) * omega[1]) %*%
    t(complex(argument = - (0:Lt) * omega[2])) *
    matrix((0:Lx)^drvX, Lx + 1, Lt + 1) *
    matrix((0:Lt)^drvT, Lx + 1, Lt + 1, byrow = TRUE)
  
  # calculate f. Note that f(-kX, - kT) = f(kX, kT) etc.
  # note that w is (Lx + 1)x(Lt + 1) but the outer values are all zeroes,
  # thus summation over them is not needed.
  fOut = 4 * sum(w[1:Lx, 1:Lt] * acfMat[1:Lx, 1:Lt]) -
    2 * ifelse(Lt > 1, sum(w[1, 2:Lt] * acfMat[1, 2:Lt]), 0) -
    2 * ifelse(Lx > 1, sum(w[2:Lx, 1] * acfMat[2:Lx, 1]), 0) -
    3 * acfMat[1, 1] * w[1, 1]
  
  return((2*pi)^(-2) * Re(fOut))
}

#----------------------Estimation of Integral over SD--------------------------#

# Function for spectral density
specIntEst00 = function(acfMat)
{
  nx = dim(as.matrix(acfMat))[1]; nt = dim(as.matrix(acfMat))[2]
  
  # calculate f. Note that f(-kX, - kT) = f(kX, kT) etc.
  iOut = 4 * sum(acfMat[1:nx, 1:nt]^2) -
    2 * ifelse(nt > 1, sum(acfMat[1, 2:nt]^2), 0) -
    2 * ifelse(nx > 1, sum(acfMat[2:nx, 1]^2), 0) - 3 * acfMat[1, 1]^2
  
  return((2*pi)^(-2) * iOut * 0.5)
}

# Function for derivativesof Spectral Density
specIntEstDrv = function(acfMat, drv1, drv2, hLag)
{
  nx = dim(as.matrix(acfMat))[1]; nt = dim(as.matrix(acfMat))[2]
  Lx = hLag[1]; Lt = hLag[2];
  
  # compute matrix of bartlett weights without "fourier factor"
  w = ((1 - (0:Lx)/max(Lx, 1)) %*% t(1 - (0:Lt)/max(Lt, 1)))^2 * 
    matrix((0:Lx)^(drv1[1] + drv2[1]), Lx + 1, Lt + 1) *
    matrix((0:Lt)^(drv1[2] + drv2[2]), Lx + 1, Lt + 1, byrow = TRUE)
  
  # calculate F. Note that f(-kX, - kT) = f(kX, kT) etc.
  FOut = 4 * sum(w[1:Lx, 1:Lt] * acfMat[1:Lx, 1:Lt]^2) -
    2 * ifelse(Lt > 1, sum(w[1, 2:Lt] * acfMat[1, 2:Lt]^2), 0) -
    2 * ifelse(Lx > 1, sum(w[2:Lx, 1] * acfMat[2:Lx, 1]^2), 0) -
    3 * w[1, 1] * acfMat[1, 1]^2
  
  return((2*pi)^(-2) * FOut)
}

#-----------------------Estimation of Bandwidths-------------------------------#

localBndwEst = function(acfMat, hLag, omega)
{
  nx = dim(as.matrix(acfMat))[1]; nt = dim(as.matrix(acfMat))[2]
  mu_2w = 4/9
  
  if(all(hLag > 1)) {
    f00 = specDensEst(acfMat, drv = c(0, 0), hLag = hLag, omega = omega)
    f10 = specDensEst(acfMat, drv = c(1, 0), hLag = hLag, omega = omega)
    f01 = specDensEst(acfMat, drv = c(0, 1), hLag = hLag, omega = omega)
  
    Lx = (4 * abs(f10^3/f01) * 1/f00^2 * nx*nt/mu_2w)^0.25
    Lt = Lx * abs(f01/f10)
  }
  else if (hLag[2] == 1)
  {
    f00 = specDensEst(acfMat, drv = c(0, 0), hLag = hLag, omega = omega)
    f10 = specDensEst(acfMat, drv = c(1, 0), hLag = hLag, omega = omega)
    
    Lx = (2 * (f10/f00)^2 * nx*nt/sqrt(mu_2w))^(1/3)
    Lt = 0
  }
  else if (hLag[1] == 1)
  {
    f00 = specDensEst(acfMat, drv = c(0, 0), hLag = hLag, omega = omega)
    f01 = specDensEst(acfMat, drv = c(0, 1), hLag = hLag, omega = omega)
    
    Lx = 0
    Lt = (2 * (f01/f00)^2 * nx*nt/sqrt(mu_2w))
  }
  
  hOut = c(trunc(Lx) + 1, trunc(Lt) + 1)
  hOut = pmin(hOut, c(nx - 1, nt - 1))
  return(hOut)
}

globalBndwEst = function(acfMat, hLag)
{
  nx = dim(as.matrix(acfMat))[1]; nt = dim(as.matrix(acfMat))[2]
  mu_2w = 4/9
  
  if(all(hLag > 1))
  {
    F00 = specIntEst00(acfMat)
    F10 = specIntEstDrv(acfMat, drv1 = c(1, 0), drv2 = c(1, 0), hLag = hLag)
    F01 = specIntEstDrv(acfMat, drv1 = c(0, 1), drv2 = c(0, 1), hLag = hLag)
    F1001 = specIntEstDrv(acfMat, drv1 = c(1, 0), drv2 = c(0, 1), hLag = hLag)
    
    Fx = F10 * (sqrt(F10/F01) + F1001/F01)/F00
    Ft = F01 * (sqrt(F01/F10) + F1001/F10)/F00
    
    Lx = (2 * Fx * nx*nt/mu_2w)^0.25
    Lt = (2 * Ft * nx*nt/mu_2w)^0.25
  } 
    else if (hLag[2] == 1)
  {
    F00 = specIntEst00(acfMat)
    F10 = specIntEstDrv(acfMat, drv1 = c(1, 0), drv2 = c(1, 0), hLag = hLag)
    
    Fx = F10/F00
    
    Lx = (2 * Fx * nx*nt/sqrt(mu_2w))^(1/3)
    Lt = 0
  }
    else if (hLag[1] == 1)
  {
    F00 = specIntEst00(acfMat)
    F01 = specIntEstDrv(acfMat, drv1 = c(0, 1), drv2 = c(0, 1), hLag = hLag)
    
    Ft = F01/F00
    
    Lx = 0
    Lt = (2 * Ft * nx*nt/mu_2w)^(1/3)
  }
  
  hLag = c(trunc(Lx) + 1, trunc(Lt) + 1)
  hLag = pmin(hLag, c(nx - 1, nt - 1))
  return(hLag)
}

#----------------------Estimation of Spectral Density--------------------------#

specDens = function(Y, omega)
{
  nx = dim(as.matrix(Y))[1]; nt = dim(as.matrix(Y))[2]
  acfMat = acfMatrix(as.matrix(Y))
  
  # initial values
  hVec = matrix(NA, 21, 2)
  hVec[1, ] = trunc(0.5*c(nx, nt))
  
  # global step
  for (g in 2:21)
  {
    hVecInfl = hVec[g - 1, ] / c(nx^(2/21), nt^(2/21))
    hVecInfl = trunc(hVecInfl) + 1
    hVec[g, ] = globalBndwEst(acfMat, hVecInfl)
    if (all(hVec[g, ] == hVec[g - 1, ])) { break() }
  }
  
  # local step
  hVecInfl = hVec[g, ] / c(nx^(2/21), nt^(2/21))
  hVecInfl = trunc(hVecInfl) + 1
  hOpt = localBndwEst(acfMat, hVecInfl, omega)
  hOpt = hOpt
  
  specDensOut = specDensEst(acfMat, drv = c(0, 0), hLag = hOpt, omega = omega)
  
  return(list(specDens = specDensOut, cf = specDensOut*(2*pi)^2, h = hOpt))
}

#--------------------------------Full SSDE-------------------------------------#

matrixSSDE = function(Y)
{
  # Preperations
  nx = dim(as.matrix(Y))[1]; nt = dim(as.matrix(Y))[2]
  acfMat = acfMatrix(as.matrix(Y))
  
  # Pre-Estimation of global bandwidth (not depending on omega)
  hVec = matrix(NA, 21, 2)
  hVec[1, ] = trunc(0.5*c(nx, nt))
  
  for (g in 2:21)
  {
    hVecInfl = hVec[g - 1, ] / c(nx^(2/21), nt^(2/21))
    hVecInfl = trunc(hVecInfl) + 1
    hVec[g, ] = globalBndwEst(acfMat, hVecInfl)
    if (all(hVec[g, ] == hVec[g - 1, ])) { break() }
  }
  
  # Set up matrices for estimation of spectrum from -pi to pi
  X = T = seq(from = -pi, to = pi, length.out = 100)
  ySpectrum = matrix(NA, nrow = 100, ncol = 100)
  
  for (i in 1:100)
  {
    for (j in 1:100)
    {
      omega = c(X[i], T[j])
      
      # estimation of local bandwidth at omega
      hVecInfl = hVec[g, ] / c(nx^(2/21), nt^(2/21))
      hVecInfl = trunc(hVecInfl) + 1
      hOpt = localBndwEst(acfMat, hVecInfl, omega)
      
      ySpectrum[i, j] = specDensEst(acfMat, drv = c(0, 0), 
                                hLag = hOpt, omega = omega)
    }
  }
  return(ySpectrum)
}