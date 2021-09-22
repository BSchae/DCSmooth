################################################################################
#                                                                              #
#                DCSmooth Package: estimation of SARMA for cf                  #
#                                                                              #
################################################################################

# This file includes all functions related to the estimation of sSARMA-processes
# for the cf coefficient of the bandwidth selection procedure

#-----------------------Calculation of cf coefficient--------------------------#

sarma2.cf = function(Y, model_order = list(ar = c(1, 1), ma = c(1, 1)))
{
  # estimation of qarma model
  sarma_model = sarma.est(Y, model_order = model_order)
  
  # get fractions for spectral density
  cf = sum(sarma_model$model$ma)^2/sum(sarma_model$model$ar)^2 *
           sarma_model$model$sigma^2
  return_list = list(cf = cf, model = sarma_model)
  
  return(return_list)
}

sarma.cf = function(Y, model_order = list(ar = c(1, 1), ma = c(1, 1)))
{
  # estimation of qarma model
  sarma_model = sarma2.est(Y, model_order = model_order)
  
  # get fractions for spectral density
  cf = sum(sarma_model$model$ma)^2/sum(sarma_model$model$ar)^2 *
    sarma_model$model$sigma^2
  return_list = list(cf = cf, model = sarma_model)
  
  return(return_list)
}

#------------------------RSS Estimation of SARMA-------------------------------#

sarma.est = function(Y, model_order = list(ar = c(1, 1), ma = c(1, 1)))
{
  model_init = suppressWarnings(sarma2.est(Y, model_order = model_order)$model)
  theta_init = c(-model_init$ar[-1, 1], -model_init$ar[1, -1],
                 model_init$ma[-1, 1], model_init$ma[1, -1])
  n_x = dim(Y)[1]; n_t = dim(Y)[2]
  # theta_init = rep(0, times = sum(unlist(model_order)))
  theta_opt  = stats::optim(theta_init, sarma_rss, R_mat = Y,
                            model_order = model_order, method = "Nelder-Mead")
  
  # put coefficients into matrices
  ar_x = c(1, -theta_opt$par[seq_len(model_order$ar[1])])
  ar_t = c(1, -theta_opt$par[model_order$ar[1] + seq_len(model_order$ar[2])])
  ma_x = c(1, theta_opt$par[sum(model_order$ar) + seq_len(model_order$ma[1])])
  ma_t = c(1, theta_opt$par[sum(model_order$ar) + model_order$ma[1] + 
                              seq_len(model_order$ma[2])])
  
  # prepare results for output
  ar_mat = ar_x %*% t(ar_t)
  ma_mat = ma_x %*% t(ma_t)
  stdev = sqrt(theta_opt$value/(n_x * n_t))
  model = list(ar = ar_mat, ma = ma_mat, sigma = stdev)
  innov = sarma.residuals(R_mat = Y, model = model)
  
  # check stationarity
  statTest = qarma.statTest(ar_mat)
  if (!statTest)
  {
    warning("SARMA model not stationary, try another order for the AR-parts.")
  }
  
  coef_out = list(Y = Y, innov = innov, model = list(ar = ar_mat, ma = ma_mat,
                  sigma = stdev), stnry = statTest)
  class(coef_out) = "qarma"
  attr(coef_out, "subclass") = "est"
  return(coef_out)
}

sarma.residuals = function(R_mat, model)
{
  n_x = dim(R_mat)[1]; n_t = dim(R_mat)[2]
  k_x = min(50, n_x); k_t = min(50, n_t)
  
  # result matrices
  E_itm = R_mat * 0   # intermediate results
  E_fnl = R_mat * 0   # final results
  
  ar_inf_x = c(1, astsa::ARMAtoAR(ar = -model$ar[-1, 1], ma = model$ma[-1, 1],
                                  lag.max = k_x))
  ar_inf_t = t(c(1, astsa::ARMAtoAR(ar = -model$ar[1, -1], ma = model$ma[1, -1],
                                    lag.max = k_t)))
  for (j in 1:n_t)
  {
    E_itm[, j] = R_mat[, j:max(1, j - k_t + 1), drop = FALSE] %*%
      ar_inf_t[1:min(j, k_t), drop = FALSE]
  }
  
  for (i in 1:n_x)
  {
    E_fnl[i, ] = ar_inf_x[1:min(i, k_x), drop = FALSE] %*%
      E_itm[i:max(1, i - k_x + 1), , drop = FALSE]
  }
  
  return(E_fnl)
}

#-----------------Fast SARMA Estimation by Hannan-Rissanen---------------------#

sarma2.est = function(Y, model_order = list(ar = c(1, 1), ma = c(1, 1)))
{
  n_x = dim(Y)[1]; n_t = dim(Y)[2]
  max_lag_x = max(model_order$ar[1], model_order$ma[1])
  max_lag_t = max(model_order$ar[2], model_order$ma[2])
  
  # estimate coefficients
  arma_x = arima(as.vector(Y),
                 order = c(model_order$ar[1], 0, model_order$ma[1]),
                 include.mean = FALSE)
  resid_x = matrix(arma_x$resid, nrow = n_x, ncol = n_t)
  arma_t = arima(as.vector(t(resid_x)),
                 order = c(model_order$ar[2], 0, model_order$ma[2]),
                 include.mean = FALSE)
  innov = matrix(arma_t$resid, nrow = n_x, ncol = n_t, byrow = TRUE)
  
  # build result matrices
  ar_mat = c(1, -arma_x$coef[seq_len(model_order$ar[1])]) %*%
           t(c(1, -arma_t$coef[seq_len(model_order$ar[2])]))
  ma_mat = c(1, arma_x$coef[model_order$ar[1] +
                             seq_len(model_order$ma[1])]) %*%
           t(c(1, arma_t$coef[model_order$ar[2] + seq_len(model_order$ma[2])]))
  
  stdev = sd(innov)
  stat_test = qarma.statTest(ar_mat)
  
  # prepare output
  coef_out = list(Y = Y, innov = innov, model = list(ar = ar_mat, ma = ma_mat,
                  sigma = stdev), stnry = stat_test)
  class(coef_out) = "qarma"
  attr(coef_out, "subclass") = "est"
  return(coef_out)
}