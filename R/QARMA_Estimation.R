################################################################################
#                                                                              #
#                DCSmooth Package: estimation of QARMA for cf                  #
#                                                                              #
################################################################################

### Part A: Estimation Functions

#----------------------Calculation of cf coefficient--------------------------#

QARMA.cf = function(Y, model = list(ar_x = 1, ar_t = 1, ma_x = 1, ma_t = 1))
{
  # estimation of qarma model
  qarma_model = QARMA.est(Y, model = model)
  
  # get variance of error terms
  sigma_sq = sd(qarma_model$estInnov)^2
  
  # get fractions for spectral density
  cf_out = sum(qarma_model$alpha)^2/sum(qarma_model$beta)^2 * sigma_sq
  
  return(cf_out)
}

#-------------------------------Model Selection--------------------------------#

#QARMA.order = 

#-----------------------------Estimation Function------------------------------#

QARMA.est = function(Y, model = list(ar_x = 1, ar_t = 1, ma_x = 1, ma_t = 1))
{
  # estimate mean of Y and calculate demeaned values E
  mu = mean(Y)
  E = Y - mu
  
  nX = dim(E)[1]; nT = dim(E)[2]
  u_est = E*0
  
  # set up matrices for lag estimation
  totalLag_ar = (model$ar_x + 1) * (model$ar_t + 1)
  totalLag_ma = (model$ma_x + 1) * (model$ma_t + 1)
  totalLag_x = max(model$ar_x, model$ma_x)
  totalLag_t = max(model$ar_t, model$ma_t)
  maxLag_ar = max(model$ar_x, model$ar_t)
  maxLag_ma = max(model$ma_x, model$ma_t)
  
  lag_E = matrix(NA, nrow = (nX - totalLag_x) * (nT - totalLag_t),
                 ncol = totalLag_ar)
  lag_u = matrix(NA, nrow = (nX - totalLag_x) * (nT - totalLag_t),
                 ncol = totalLag_ma)
  
  # fill AR-lag matrix
  lagNames_ar = 1:totalLag_ar*NA
  
  for (i in 1:(model$ar_x + 1))
  {
    for (j in 1:(model$ar_t + 1))
    {
      indexX = (totalLag_x - i + 2):(nX - i + 1)
      indexT = (totalLag_t - j + 2):(nT - j + 1)
      lag_E[, (i - 1)*(maxLag_ar + 1) + j] = as.vector(E[indexX, indexT])
      lagNames_ar[(i - 1)*(maxLag_ar + 1) + j] = paste0("lag", i - 1, j - 1)
    }
  }
  colnames(lag_E) = lagNames_ar
  
  # backcasting iteration
  for(s in 1:3)
  {
    # estimation of coefficients
    if (totalLag_ar == 1 && totalLag_ma > 1) {                  # QMA process
      if (s == 1) {
        coefs = rep(0, times = totalLag_ma - 1)
      } else {
        coefs = lm(lag_E[, 1] ~ lag_u[, totalLag_ma:2] + 0)$coef
      }
      
      # fill matrices alpha and beta
      alpha_est = t(matrix(c(coefs[1:(totalLag_ma - 1)], 0), 
                  nrow = (model$ma_x + 1), ncol = (model$ma_t + 1)))
      beta_est = as.matrix(0)   # no AR-parameter for lags > 0
      
    } else if (totalLag_ar > 1 && totalLag_ma == 1) {           # QAR process
      coefs = lm(lag_E[, 1] ~ lag_E[, totalLag_ar:2] + 0)$coef
      
      # fill matrices alpha and beta
      alpha_est = as.matrix(0)   # no MA-parameter for lags > 0
      beta_est = matrix(c(coefs[1:(totalLag_ar - 1)], 0), 
                        nrow = (model$ar_x + 1), ncol = (model$ar_t + 1),
                        byrow = TRUE)
      
    } else if (totalLag_ar > 1 && totalLag_ma > 1) {            # QARMA process
      if (s == 1) {
        coefs = c(lm(lag_E[, 1] ~ lag_E[, totalLag_ar:2] + 0)$coef,
                  rep(0, times = totalLag_ma - 1))
      } else {
        coefs = lm(lag_E[, 1] ~ lag_E[, totalLag_ar:2] +
                   lag_u[, totalLag_ma:2] + 0)$coef
      }
      # fill matrices alpha and beta
      beta_est = matrix(c(coefs[1:(totalLag_ar - 1)], 0), 
                        nrow = (model$ar_x + 1), ncol = (model$ar_t + 1),
                        byrow = TRUE)
      alpha_est = matrix(c(coefs[totalLag_ar:(totalLag_ar +
                        totalLag_ma - 2)], 0), nrow = (model$ma_x + 1),
                        ncol = (model$ma_t + 1))
    } else {
      coefs = 0
    }
    
    # backcasting procedure (not needed in the last iteration for s = 3)
    if (s != 3)
    {
      # backcasting of QARMA process in matrix form
      for (i in (totalLag_x + 1):nX)
      {
        for (j in (totalLag_t + 1):nT)
        {
          u_est[i, j] = E[i, j] - sum(beta_est * E[(i - model$ar_x):i,
                                                   (j - model$ar_t):j]) -
            sum(alpha_est * u_est[(i - model$ma_x):i,
                                  (j - model$ma_t):j])
        }
      }
      
      # fill MA-lag matrix
      lagNames_ma = 1:totalLag_ma*NA
      for (i in 1:(model$ma_x + 1))
      {
        for (j in 1:(model$ma_t + 1))
        {
          indexX = (totalLag_x - i + 2):(nX - i + 1)
          indexT = (totalLag_t - j + 2):(nT - j + 1)
          lag_u[,(i - 1)*(model$ma_t + 1) + j] = as.vector(u_est[indexX, indexT])
          lagNames_ma[(i - 1)*(model$ma_t + 1) + j] = paste0(i - 1, ",", j - 1)
        }
      }
      colnames(lag_u) = lagNames_ma
    }
    
    # check stationarity
    statTest = QARMA.statTest(beta_est)
    if (statTest == FALSE)
    {
      stop("QARMA model not stationary, try another order for the AR-parts.")
    }
  }
  # preparation of output
  beta_est[model$ar_x + 1, model$ar_t + 1] = -1
  alpha_est[model$ma_x + 1, model$ma_t + 1] = 1
  coefOut = list(mu = mu, beta = beta_est, alpha = alpha_est, estInnov = u_est)
  
  return(coefOut)
}

# QARMA.est2 = function(Y, model = list(ar_x = 1, ar_t = 1, ma_x = 1, ma_t = 1))
# {
#   # estimate mean of Y and calculate demeaned values E
#   mu = mean(Y)
#   E = Y - mu
#   
#   nX = dim(E)[1]; nT = dim(E)[2]
#   u_est = E*0
#   
#   ### set up matrices and parameters for lag estimation
#   totalLag_ar = (model$ar_x + 1) * (model$ar_t + 1) 
#                             # number of total lags (including 0) for the AR part
#   totalLag_ma = (model$ma_x + 1) * (model$ma_t + 1)
#                             # number of total lags (including 0) for the MA part
#   maxLag_x = max(model$ar_x, model$ma_x)  # maximum lag in x-direction
#   maxLag_t = max(model$ar_t, model$ma_t)  # maximum lag in t-direction
#   maxLag_ar = max(model$ar_x, model$ar_t) # maximum AR-lag
#   maxLag_ma = max(model$ma_x, model$ma_t) # maximum MA-lag
#   
#   lag_E = matrix(NA, nrow = (nX - maxLag_x) * (nT - maxLag_t),
#                  ncol = totalLag_ar)  # lag matrix for observations
#   lag_u = matrix(NA, nrow = (nX - maxLag_x) * (nT - maxLag_t),
#                  ncol = totalLag_ma)  # lag matrix for innovations
#   
#   ### fill AR-lag matrix
#   lagNames_ar = 1:totalLag_ar*NA # names for matrix
#   for (i in 1:(model$ar_x + 1))
#   {
#     for (j in 1:(model$ar_t + 1))
#     {
#       indexX = (maxLag_x - i + 2):(nX - i + 1)
#       indexT = (maxLag_t - j + 2):(nT - j + 1)
#       lag_E[, (i - 1)*(maxLag_ar + 1) + j] = as.vector(E[indexX, indexT])
#       lagNames_ar[(i - 1)*(maxLag_ar + 1) + j] = paste0("lag", i - 1, j - 1)
#     }
#   }
#   colnames(lag_E) = lagNames_ar
#   
#   ### auxiliary regression for calculation of innovations (from residuals)
#   aux_coefs = lm(lag_E[, 1] ~ lag_E[, totalLag_ar:2] + 0)$coef
#   aux_coefs[which(is.na(aux_coefs))] = 0
#   aux_beta = matrix(c(aux_coefs[1:(totalLag_ar - 1)], 0), 
#                 nrow = (model$ar_x + 1), ncol = (model$ar_t + 1), byrow = TRUE)
#   
#   for (s in 1:5) {
#   # backcasting of ARMA process
#   for (i in (maxLag_x + 1):nX)
#   {
#     for (j in (maxLag_t + 1):nT)
#     {
#       u_est[i, j] = E[i, j] - sum(aux_beta * E[(i - model$ar_x):i,
#                                                (j - model$ar_t):j])
#     }
#   }
#   
#   # fill MA-lag matrix
#   lagNames_ma = 1:totalLag_ma*NA
#   for (i in 1:(model$ma_x + 1))
#   {
#     for (j in 1:(model$ma_t + 1))
#     {
#       indexX = (maxLag_x - i + 2):(nX - i + 1)
#       indexT = (maxLag_t - j + 2):(nT - j + 1)
#       lag_u[ ,(i - 1)*(model$ma_t + 1) + j] = as.vector(u_est[indexX, indexT])
#       lagNames_ma[(i - 1)*(model$ma_t + 1) + j] = paste0("lag", i - 1, j - 1)
#     }
#   }
#   colnames(lag_u) = lagNames_ma
#   
#   # Parameter Estimation
#   if (totalLag_ar == 1 && totalLag_ma > 1) {
#     coefs = lm(lag_E[, 1] ~ lag_u[, totalLag_ma:2] + 0)$coef
#   } else if (totalLag_ar > 1 && totalLag_ma == 1) {
#     coefs = lm(lag_E[, 1] ~ lag_E[, totalLag_ar:2] + 0)$coef
#   } else if (totalLag_ar > 1 && totalLag_ma > 1) {
#     coefs = lm(lag_E[, 1] ~ lag_E[, totalLag_ar:2] +
#                  lag_u[, totalLag_ma:2] + 0)$coef
#   } else {
#     coefs = 0
#   }
#   
#   # calculate coefficient matrices
#   beta_est = matrix(c(coefs[1:(totalLag_ar - 1)], -1), 
#                 nrow = (model$ar_x + 1), ncol = (model$ar_t + 1), byrow = TRUE)
#   alpha_est = matrix(c(coefs[totalLag_ar:(totalLag_ar + totalLag_ma - 2)], 1), 
#                 nrow = (model$ar_x + 1), ncol = (model$ar_t + 1), byrow = TRUE)
#   
#   # check stationarity
#   statTest = QARMA.statTest(beta_est)
#   if (statTest == FALSE)
#   {
#     stop("QARMA model not stationary, try another order for the AR-parts.")
#   }
#   }
#   # preparation of output
#   coefOut = list(mu = mu, beta = beta_est, alpha = alpha_est, estInnov = u_est)
#   
#   return(coefOut)
# }

#----------------------------Simulation Function------------------------------#

QARMA.sim = function(nX, nT, model = list(alpha, beta, stdev))
{
  alpha = as.matrix(model$alpha); beta = as.matrix(model$beta)
  ar_x = dim(beta)[1] - 1; ar_t = dim(beta)[2] - 1
  ma_x = dim(alpha)[1] - 1; ma_t = dim(alpha)[2] - 1
  xInit = max(ar_x, ma_x) + 1
  tInit = max(ar_t, ma_t) + 1
  
  # check matrices alpha, beta
  alpha[ma_x + 1, ma_t + 1] = 1 # MA-coefficients
  beta[ar_x + 1, ar_t + 1] = 0  # AR-coefficients
  
  errorMat = matrix(rnorm(4*nX*nT), nrow = 2*nX, ncol = 2*nT) * model$stdev
  armaMat = errorMat
  
  for (i in xInit:(2*nX))
  {
    for (j in tInit:(2*nT))
    {
      armaMat[i, j] = sum(beta * armaMat[(i - ar_x):i, (j - ar_t):j]) +
                      sum(alpha * errorMat[(i - ma_x):i, (j - ma_t):j])
    }
  }
  
  armaOut = armaMat[(nX + 1):(2 * nX), (nT + 1):(2 * nT)]
  errorOut = errorMat[(nX + 1):(2 * nX), (nT + 1):(2 * nT)]
  
  listOut = list(qarma.ts = armaOut, innov = errorOut,
                  alpha = alpha, beta = beta, stdev = model$stdev)
  return(listOut)
}

#-------------------------Test for QARMA stationarity-------------------------#

QARMA.statTest = function(beta)
{
  beta[dim(beta)[1], dim(beta)[2]] = -1
  outValue = TRUE
  
  # compute reference sign at center (0, 0)
  signRef = sign(QARMA.auxF(0, 0, beta))
  
  # check for random numbers inside unit circle
  
  if (outValue == TRUE)
  {
    for (k in 1:1000)
    {
      # draw random point inside unit circle
      phi = runif(2) * 2*pi
      rad = sqrt(runif(2)) * 1.01
      x = rad * cos(phi)
      y = rad * sin(phi)
      z1c = complex(real = x[1], imaginary = y[1])
      z2c = complex(real = x[2], imaginary = y[2])
      
      signTest = sign(Re(QARMA.auxF(z1c, z2c, beta)))
      if (signTest != signRef)
      {
        outValue = FALSE
        break()
      }
    }
  }
  return(outValue)
}

#-----------------------Characteristic Function of AR Part--------------------#

QARMA.auxF = function(z1c, z2c, beta)
{
  ar_x = dim(beta)[1] - 1; ar_t = dim(beta)[2] - 1
  
  # set up vectors for z^i
  zx_vec = z1c^(ar_x:0)
  zt_vec = z2c^(ar_t:0)
  
  out = t(zx_vec) %*% beta %*% zt_vec
  
  return(out)
}

qarmaSSDE = function(Y, alpha, beta, sigmaSq)
{
  X = T = seq(from = -pi, to = pi, length.out = 100)
  ySpectrum = matrix(NA, nrow = 100, ncol = 100)
  
  for (i in 1:100)
  {
    for (j in 1:100)
    {
      omega = c(X[i], T[j])
      ySpectrum[i, j] = QARMA.spec(Y, alpha, beta, sigmaSq, omega)
    }
  }
  return(ySpectrum)
}

#----------------------Calculation of spectral density-------------------------#

QARMA.spec = function(Y, alpha, beta, sigmaSq, omega)
{
  lagAR = dim(beta) - 1; lagMA = dim(alpha) - 1
  z1 = complex(argument = omega[1]); z2 = complex(argument = omega[2])
  
  g00 = sd(Y)^2
  
  betaSum = Re(sum(beta * (z1^(lagAR[1]:0)) %*% t(z2^(lagAR[2]:0))) * 
                 sum(beta * ((1/z1)^(lagAR[1]:0)) %*% t((1/z2)^(lagAR[2]:0))))
  alphaSum = Re(sum(alpha * (z1^(lagMA[1]:0)) %*% t(z2^(lagMA[2]:0))) * 
                  sum(alpha * ((1/z1)^(lagMA[1]:0)) %*% t((1/z2)^(lagMA[2]:0))))
  
  specDensOut = sigmaSq/g00 * betaSum/alphaSum
  
  return(specDensOut)
}
