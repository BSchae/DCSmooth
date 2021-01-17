################################################################################
#                                                                              #
#                DCSmooth Package: estimation of QARMA for cf                  #
#                                                                              #
################################################################################

### Part A: Estimation Functions

#-----------------------Calculation of cf coefficient--------------------------#

QARMA.cf = function(Y, model_order = list(ar = c(1, 1), ma = c(1, 1)))
{
  # estimation of qarma model
  qarma_model = QARMA.est(Y, model_order = model_order)
  
  # get variance of error terms
  sigma_sq = sd(qarma_model$innov)^2
  
  # get fractions for spectral density
  cf_out = sum(qarma_model$ma)^2/sum(qarma_model$ar)^2 * sigma_sq
  
  return(cf_out)
}

#-------------------------------Model Selection--------------------------------#

#QARMA.order = 

#-----------------------------Estimation Function------------------------------#

QARMA.est = function(Y, model_order = list(ar = c(1, 1), ma = c(1, 1)))
{
  #Y = matrix(rnorm(10000), 100, 100)
  
  nX = dim(Y)[1]; nT = dim(Y)[2]
  
  # calculate ARMA-orders for different purposes
  max_lag_x = max(model_order$ar[1], model_order$ma[1])
  max_lag_t = max(model_order$ar[2], model_order$ma[2])
  total_lag_ar = (model_order$ar[1] + 1) * (model_order$ar[2] + 1) - 1
  total_lag_ma = (model_order$ma[1] + 1) * (model_order$ma[2] + 1) - 1
  max_lag_ar = max(model_order$ar)
  max_lag_ma = max(model_order$ma)
  
  ### auxiliary estimation model
  matrix_aux = data.frame(matrix(NA, nrow = (nX - max_lag_x) * (nT - max_lag_t),
                    ncol = total_lag_ar + 1))       # +1 is to include lag_ar_00

  for (i in 0:model_order$ar[1]) # use lag orders as loop indices
  {
    for (j in 0:model_order$ar[2])
    {
      index_x = (max_lag_x - i + 1):(nX - i)
      index_t = (max_lag_t - j + 1):(nT - j)
      matrix_aux[, (j + 1) + i * (model_order$ar[2] + 1)] =  # over colums first (inner loop)
        as.vector(Y[index_x, index_t])
      names(matrix_aux)[(j + 1) + i * (model_order$ar[2] + 1)] =
        paste0("lag_ar_", i, j)
    }
  }
  
  # linear regression with AR-part only
  aux_reg = lm(lag_ar_00 ~ . + 0, data = matrix_aux)
  residuals_arma = matrix(0, nrow = nX, ncol = nT)
  
  # set up matrix for complete arma
  matrix_arma = data.frame(matrix(NA, nrow = (nX - max_lag_x) *
                    (nT - max_lag_t), ncol = total_lag_ar + total_lag_ma + 1))
  matrix_arma[, 1:(total_lag_ar + 1)] = matrix_aux
  names(matrix_arma)[1:(total_lag_ar + 1)] = names(matrix_aux)
  
  # iteration loop (number of iterations might increase the precision)
  if (total_lag_ma > 0)
  {
    for (loop in 1:3)
    {
      ### build nX*nT matrix of residuals
      residuals_arma[(max_lag_x + 1):nX, (max_lag_t + 1):nT] = 
        matrix(aux_reg$residuals, nrow = (nX - max_lag_x),
              ncol = (nT - max_lag_t))
      
  
      ### fill matrix_arma with lagged residuals
      for (i in 0:model_order$ma[1]) # use lag orders as loop indices
      {
        for (j in 0:model_order$ma[2])
        {
          if (!(i == 0 && j == 0))  # lag_ma_00 is not needed (new residuals)
          {
            index_x = (max_lag_x - i + 1):(nX - i)
            index_t = (max_lag_t - j + 1):(nT - j)
            matrix_arma[, total_lag_ar + (j + 1) + i *
                          (model_order$ma[2] + 1)] = 
              as.vector(residuals_arma[index_x, index_t])
            names(matrix_arma)[total_lag_ar + (j + 1) + i *
                                (model_order$ma[2] + 1)] =
              paste0("lag_ma_", i, j)
          }
        }
      }
      
      # regression for iterated estimation of residuals
      aux_reg = lm(lag_ar_00 ~ . + 0, data = matrix_arma)
    }
  }
  else {
    ### build nX*nT matrix of residuals
    residuals_arma[(max_lag_x + 1):nX, (max_lag_t + 1):nT] = 
      matrix(aux_reg$residuals, nrow = (nX - max_lag_x),
             ncol = (nT - max_lag_t))
  }
  
  # byrow = TRUE and reverse needed as lm order is 01,10,11...
  ar_mat = matrix(c(aux_reg$coef[total_lag_ar:1], -1), byrow = TRUE,
                nrow = (model_order$ar[1] + 1), ncol = (model_order$ar[2] + 1))
  ar_mat[model_order$ar[1] + 1, model_order$ar[2] + 1] = -1
  ma_mat = matrix(c(aux_reg$coef[(total_lag_ar + total_lag_ma):
                  (total_lag_ar + 1)], 1), byrow = TRUE,
                  nrow = (model_order$ma[1] + 1), ncol = 
                    (model_order$ma[2] + 1))
  ma_mat[model_order$ma[1] + 1, model_order$ma[2] + 1] = 1
  innov = residuals_arma
  
  # check stationarity
  statTest = QARMA.statTest(ar_mat)
  if (statTest == FALSE)
  {
    stop("QARMA model not stationary, try another order for the AR-parts.")
  }
  
  # preparation of output
  coef_out = list(ar = ar_mat, ma = ma_mat, innov = innov)
  return(coef_out)
}

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
