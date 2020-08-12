###############################################################################
#                                                                             #
#                DCSmooth Package: stimation of QARMA for cf                  #
#                                                                             #
###############################################################################


#----------------------------Estimation Function------------------------------#

QARMA.est = function(E, model = list(ar_x = 1, ar_t = 1, ma_x = 1, ma_t = 1))
{
  nX = dim(E)[1]; nT = dim(E)[2]
  
  u_est = E*0
  
  # set up matrices for lag estimation
  totalLag_ar = (model$ar_x + 1) * (model$ar_t + 1)
  totalLag_ma = (model$ma_x + 1) * (model$ma_t + 1)
  totalLag_x = max(model$ar_x, model$ma_x)
  totalLag_t = max(model$ar_t, model$ma_t)
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
      lag_E[, (i - 1)*(totalLag_t + 1) + j] = as.vector(E[indexX, indexT])
      lagNames_ar[(i - 1)*(totalLag_t + 1) + j] = paste0("lag", i - 1, j - 1)
    }
  }
  colnames(lag_E) = lagNames_ar
  
  #backcasting iteration
  for(s in 1:3)
  {
    # fill MA-lag matrix
    lagNames_ma = 1:totalLag_ma*NA
    for (i in 1:(model$ma_x+ 1))
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
    
    # preparing for backcasting
    coefs = lm(lag_E[, 1] ~ lag_E[, totalLag_ar:2] +
                            lag_u[, totalLag_ma:2] + 0)$coef
    coefs[which(is.na(coefs))] = 0
    beta_est = matrix(c(coefs[1:(totalLag_ar - 1)], 0), 
                  nrow = (model$ar_x + 1), ncol = (model$ar_t + 1), byrow = TRUE)
    alpha_est = t(matrix(c(coefs[totalLag_ar:(totalLag_ar +
                  totalLag_ma - 2)], 0), nrow = (model$ma_x + 1),
                  ncol = (model$ma_t + 1)))
    
    if (s != 3)
    {
      # backcasting of ARMA process
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
  coefOut = list(beta = beta_est, alpha = alpha_est, estInnov = u_est)
  
  return(coefOut)
}

#----------------------Calculation of cf coefficient--------------------------#

QARMA.cf = function(E, model = list(ar_x = 1, ar_t = 1, ma_x = 1, ma_t = 1))
{
  # estimation of qarma model
  qarma_model = QARMA.est(E, model = model)
  
  # get variance of error terms
  sigma_sq = sd(qarma_model$estInnov)^2
  
  # get fractions for spectral density
  cf_out = sum(qarma_model$alpha)^2/sum(qarma_model$beta)^2 * sigma_sq
  
  return(cf_out)
}

#----------------------------Simulation Function------------------------------#

QARMA.sim = function(nX, nT, model = list(alpha, beta, stdev))
{
  ar_x = dim(beta)[1] - 1; ar_t = dim(beta)[2] - 1
  ma_x = dim(alpha)[1] - 1; ma_t = dim(beta)[2] - 1
  xInit = max(ar_x, ma_x) + 1
  tInit = max(ar_t, ma_t) + 1
  
  # check matrices alpha, beta
  alpha[ma_x + 1, ma_t + 1] = 1
  beta[ar_x + 1, ar_t + 1] = 0
  
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
                  alpha = alpha, beta = beta, sd = sd)
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
