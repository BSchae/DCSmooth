################################################################################
#                                                                              #
#                DCSmooth Package: estimation of QARMA for cf                  #
#                                                                              #
################################################################################

### Part A: Estimation Functions

#-----------------------Calculation of cf coefficient--------------------------#

QARMA.cf = function(Y, order_x = c(1, 1), order_t = c(1, 1))
{
  # estimation of qarma model
  qarma_model = QARMA.est(Y, order_x = order_x, order_t = order_t)
  
  # get variance of error terms
  sigma_sq = sd(qarma_model$innov)^2
  
  # get fractions for spectral density
  cf_out = sum(qarma_model$ma)^2/sum(qarma_model$ar)^2 * sigma_sq
  
  return(cf_out)
}

#-------------------------------Model Selection--------------------------------#

#QARMA.order = 

#-----------------------------Estimation Function------------------------------#

QARMA.est = function(Y, order_x = c(1, 1), order_t = c(1, 1))
{
  #Y = matrix(rnorm(10000), 100, 100)
  
  nX = dim(Y)[1]; nT = dim(Y)[2]
  u_est = Y * 0
  
  # calculate ARMA-orders for different purposes
  max_lag_x = max(order_x)
  max_lag_t = max(order_t)
  total_lag_ar = (order_x[1] + 1) * (order_t[1] + 1) - 1
  total_lag_ma = (order_x[2] + 1) * (order_t[2] + 1) - 1
  max_lag_ar = max(order_x[1], order_t[1])
  max_lag_ma = max(order_x[2], order_t[2])
  
  ### auxiliary estimation model
  matrix_aux = data.frame(matrix(NA, nrow = (nX - max_lag_x) * (nT - max_lag_t),
                    ncol = total_lag_ar + 1))       # +1 is to include lag_ar_00

  for (i in 0:order_x[1]) # use lag orders as loop indices
  {
    for (j in 0:order_t[1])
    {
      index_x = (max_lag_x - i + 1):(nX - i)
      index_t = (max_lag_t - j + 1):(nT - j)
      matrix_aux[, (j + 1) + i * (order_t[1] + 1)] =  # over colums first (inner loop)
        as.vector(Y[index_x, index_t])
      names(matrix_aux)[(j + 1) + i * (order_t[1] + 1)] =
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
  for (loop in 1:5)
  {
    ### build nX*nT matrix of residuals
    residuals_arma[(max_lag_x + 1):nX, (max_lag_t + 1):nT] = 
      matrix(aux_reg$residuals, nrow = (nX - max_lag_x),
             ncol = (nT - max_lag_t))
    
    ### fill matrix_arma with lagged residuals
    for (i in 0:order_x[2]) # use lag orders as loop indices
    {
      for (j in 0:order_t[2])
      {
        if (!(i == 0 && j == 0))  # lag_ma_00 is not needed (new residuals)
        {
          index_x = (max_lag_x - i + 1):(nX - i)
          index_t = (max_lag_t - j + 1):(nT - j)
          matrix_arma[, total_lag_ar + (j + 1) + i * (order_t[1] + 1)] =
            as.vector(residuals_arma[index_x, index_t])
          names(matrix_arma)[total_lag_ar + (j + 1) + i * (order_t[1] + 1)] =
            paste0("lag_ma_", i, j)
        }
      }
    }
    
    # regression for iterated estimation of residuals
    aux_reg = lm(lag_ar_00 ~ . + 0, data = matrix_arma)
    
  }
  
  # byrow = TRUE and reverse needed as lm order is 01,10,11...
  ar_mat = matrix(c(aux_reg$coef[total_lag_ar:1], -1), byrow = TRUE,
                  nrow = (order_x[1] + 1), ncol = (order_t[1] + 1))
  ar_mat[order_x[1] + 1, order_t[1] + 1] = -1
  ma_mat = matrix(c(aux_reg$coef[(total_lag_ar + total_lag_ma):
                  (total_lag_ar + 1)], 1), byrow = TRUE,
                  nrow = (order_x[2] + 1), ncol = (order_t[2] + 1))
  ma_mat[order_x[2] + 1, order_t[2] + 1] = 1
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

# QARMA.est = function(Y, model = list(ar_x = 1, ar_t = 1, ma_x = 1, ma_t = 1))
# {
#   # estimate mean of Y and calculate demeaned values E
#   mu = mean(Y)
#   E = Y - mu
#   
#   nX = dim(E)[1]; nT = dim(E)[2]
#   u_est = E*0
#   
#   # set up matrices for lag estimation
#   totalLag_ar = (model$ar_x + 1) * (model$ar_t + 1)
#   totalLag_ma = (model$ma_x + 1) * (model$ma_t + 1)
#   totalLag_x = max(model$ar_x, model$ma_x)
#   totalLag_t = max(model$ar_t, model$ma_t)
#   maxLag_ar = max(model$ar_x, model$ar_t)
#   maxLag_ma = max(model$ma_x, model$ma_t)
#   
#   lag_E = matrix(NA, nrow = (nX - totalLag_x) * (nT - totalLag_t),
#                  ncol = totalLag_ar)
#   lag_u = matrix(NA, nrow = (nX - totalLag_x) * (nT - totalLag_t),
#                  ncol = totalLag_ma)
#   
#   # fill AR-lag matrix
#   lagNames_ar = 1:totalLag_ar*NA
#   
#   for (i in 1:(model$ar_x + 1))
#   {
#     for (j in 1:(model$ar_t + 1))
#     {
#       indexX = (totalLag_x - i + 2):(nX - i + 1)
#       indexT = (totalLag_t - j + 2):(nT - j + 1)
#       lag_E[, (i - 1)*(maxLag_ar + 1) + j] = as.vector(E[indexX, indexT])
#       lagNames_ar[(i - 1)*(maxLag_ar + 1) + j] = paste0("lag", i - 1, j - 1)
#     }
#   }
#   colnames(lag_E) = lagNames_ar
#   
#   # backcasting iteration
#   for(s in 1:3)
#   {
#     # estimation of coefficients
#     if (totalLag_ar == 1 && totalLag_ma > 1) {                  # QMA process
#       if (s == 1) {
#         coefs = rep(0, times = totalLag_ma - 1)
#       } else {
#         coefs = lm(lag_E[, 1] ~ lag_u[, totalLag_ma:2] + 0)$coef
#       }
#       
#       # fill matrices alpha and beta
#       alpha_est = t(matrix(c(coefs[1:(totalLag_ma - 1)], 0), 
#                            nrow = (model$ma_x + 1), ncol = (model$ma_t + 1)))
#       beta_est = as.matrix(0)   # no AR-parameter for lags > 0
#       
#     } else if (totalLag_ar > 1 && totalLag_ma == 1) {           # QAR process
#       coefs = lm(lag_E[, 1] ~ lag_E[, totalLag_ar:2] + 0)$coef
#       
#       # fill matrices alpha and beta
#       alpha_est = as.matrix(0)   # no MA-parameter for lags > 0
#       beta_est = matrix(c(coefs[1:(totalLag_ar - 1)], 0), 
#                         nrow = (model$ar_x + 1), ncol = (model$ar_t + 1),
#                         byrow = TRUE)
#       
#     } else if (totalLag_ar > 1 && totalLag_ma > 1) {            # QARMA process
#       if (s == 1) {
#         coefs = c(lm(lag_E[, 1] ~ lag_E[, totalLag_ar:2] + 0)$coef,
#                   rep(0, times = totalLag_ma - 1))
#       } else {
#         coefs = lm(lag_E[, 1] ~ lag_E[, totalLag_ar:2] +
#                      lag_u[, totalLag_ma:2] + 0)$coef
#       }
#       # fill matrices alpha and beta
#       beta_est = matrix(c(coefs[1:(totalLag_ar - 1)], 0), 
#                         nrow = (model$ar_x + 1), ncol = (model$ar_t + 1),
#                         byrow = TRUE)
#       alpha_est = matrix(c(coefs[totalLag_ar:(totalLag_ar +
#                                                 totalLag_ma - 2)], 0), nrow = (model$ma_x + 1),
#                          ncol = (model$ma_t + 1))
#     } else {
#       coefs = 0
#     }
#     
#     # backcasting procedure (not needed in the last iteration for s = 3)
#     if (s != 3)
#     {
#       # backcasting of QARMA process in matrix form
#       for (i in (totalLag_x + 1):nX)
#       {
#         for (j in (totalLag_t + 1):nT)
#         {
#           u_est[i, j] = E[i, j] - sum(beta_est * E[(i - model$ar_x):i,
#                                                    (j - model$ar_t):j]) -
#             sum(alpha_est * u_est[(i - model$ma_x):i,
#                                   (j - model$ma_t):j])
#         }
#       }
#       
#       # fill MA-lag matrix
#       lagNames_ma = 1:totalLag_ma*NA
#       for (i in 1:(model$ma_x + 1))
#       {
#         for (j in 1:(model$ma_t + 1))
#         {
#           indexX = (totalLag_x - i + 2):(nX - i + 1)
#           indexT = (totalLag_t - j + 2):(nT - j + 1)
#           lag_u[,(i - 1)*(model$ma_t + 1) + j] = as.vector(u_est[indexX, indexT])
#           lagNames_ma[(i - 1)*(model$ma_t + 1) + j] = paste0(i - 1, ",", j - 1)
#         }
#       }
#       colnames(lag_u) = lagNames_ma
#     }
#     
#     # check stationarity
#     statTest = QARMA.statTest(beta_est)
#     if (statTest == FALSE)
#     {
#       stop("QARMA model not stationary, try another order for the AR-parts.")
#     }
#   }
#   # preparation of output
#   beta_est[model$ar_x + 1, model$ar_t + 1] = -1
#   alpha_est[model$ma_x + 1, model$ma_t + 1] = 1
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
