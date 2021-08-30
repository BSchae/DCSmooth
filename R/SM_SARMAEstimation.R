################################################################################
#                                                                              #
#                DCSmooth Package: estimation of SARMA for cf                  #
#                                                                              #
################################################################################

# This file includes all functions related to the estimation of sSARMA-processes
# for the cf coefficient of the bandwidth selection procedure

#-----------------------Calculation of cf coefficient--------------------------#

sarma.cf = function(Y, model_order = list(ar = c(1, 1), ma = c(1, 1)))
{
  # estimation of qarma model
  sarma_est = sarma.est(Y, model_order = model_order)
  
  # get fractions for spectral density
  cf = sum(sarma_est$model$ma)^2/sum(sarma_est$model$ar)^2 *
    sarma_est$model$sigma^2
  sarma_model = sarma_est$model
  return_list = list(cf = cf, model = sarma_model)
  
  return(return_list)
}

#------------------------RSS Estimation of SARMA-------------------------------#

#' @export
sarma.est = function(Y, model_order = list(ar = c(1, 1), ma = c(1, 1)))
{
  model_init = qarma.est(Y, model_order = model_order)$model
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
    warning("QARMA model not stationary, try another order for the AR-parts.")
  }
  
  coef_out = list(Y = Y, innov = innov, model = list(ar = ar_mat, ma = ma_mat,
                  sigma = stdev), stnry = statTest)
  class(coef_out) = "qarma"
  attr(coef_out, "subclass") = "est"
  return(coef_out)
}

sarma.rss = function(theta, R_mat,
                     model_order = list(ar = c(1, 1), ma = c(1, 1)))
{
  n_x = dim(R_mat)[1]; n_t = dim(R_mat)[2]
  k_x = min(50, n_x); k_t = min(50, n_t)
  
  # get coefficients from theta
  # theta = (ar_x, ar_t, ma_x, ma_t)
  ar_x = theta[seq_len(model_order$ar[1])]
  ar_t = theta[model_order$ar[1] + seq_len(model_order$ar[2])]
  ma_x = theta[sum(model_order$ar) + seq_len(model_order$ma[1])]
  ma_t = theta[sum(model_order$ar) + model_order$ma[1] + 
                 seq_len(model_order$ma[2])]
  
  # result matrices
  E_itm = R_mat * 0   # intermediate results
  E_fnl = R_mat * 0   # final results
  
  # Two-step estimation of e_ij
  ar_inf_x = c(1, astsa::ARMAtoAR(ar = ar_x, ma = ma_x, lag.max = k_x))
  ar_inf_t = t(c(1, astsa::ARMAtoAR(ar = ar_t, ma = ma_t, lag.max = k_t)))
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
  
  RSS = sum(E_fnl^2)
  
  return(RSS)
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

# sarma2.est = function(Y, model_order = list(ar = c(1, 1), ma = c(1, 1)))
# {
#   nX = dim(Y)[1]; nT = dim(Y)[2]
#   
#   # calculate ARMA-orders for different purposes
#   max_lag_x = max(model_order$ar[1], model_order$ma[1])
#   max_lag_t = max(model_order$ar[2], model_order$ma[2])
#   total_lag_ar = (model_order$ar[1] + 1) * (model_order$ar[2] + 1) - 1
#   total_lag_ma = (model_order$ma[1] + 1) * (model_order$ma[2] + 1) - 1
#   max_lag_ar = max(model_order$ar)
#   max_lag_ma = max(model_order$ma)
#   
#   ### AR-only auxiliary estimation model by YW-estimation ###
#   m1_ar = max(max_lag_x, 1) + ifelse(total_lag_ma > 0, 2, 0) 
#   m2_ar = max(max_lag_t, 1) + ifelse(total_lag_ma > 0, 2, 0)
#   
#   # ar1 = ar.yw(as.vector(Y), order.max = m1_ar, aic = FALSE, demean = FALSE)
#   # rmat1 = matrix(ar1$resid, nT, nX, byrow = TRUE)
#   # rmat1[1, 1:m1_ar] = 0
#   # ar2 = ar.yw(as.vector(rmat1), order.max = m2_ar, aic = FALSE, demean = FALSE)
#   # ar_x = c(1, -ar1$ar)
#   # ar_t = c(1, -ar2$ar)
#   
#   ar_x = c(1, -ar.yw(as.vector(Y), order = m1_ar, aic = FALSE,
#                      demean = FALSE)$ar)
#   ar_t = c(1, -ar.yw(as.vector(t(Y)), order = m2_ar, aic = FALSE,
#                      demean = FALSE)$ar)
#   ar_aux = ar_x %*% t(ar_t)
#   
#   # calculate residuals from auxiliary AR model
#   R_mat = matrix(0, nrow = nX - m1_ar, ncol = nT - m2_ar)
#   for (i in 1:(nX - m1_ar))
#   {
#     for (j in 1:(nT - m2_ar))
#     {
#       R_mat[i, j] = sum(ar_aux * Y[(i + m1_ar):i, (j + m2_ar):j])
#     }
#   }
#   
#   if (total_lag_ma > 0)
#   {
#     # set up submatrices of Y, res for use in fill_matrix procedure
#     Y_mat = Y[(m1_ar + 1):nX, (m2_ar + 1):nT]
#     
#     matrix_arma_x = .sarma.lag_matrix(Y_mat, R_mat, 
#                                       lag_ar = model_order$ar[1],
#                                       lag_ma = model_order$ma[1])
#     arma_reg_x = stats::lm(ar_0 ~ . + 0, data = matrix_arma_x)
#     R_mat_x = arma_reg
#     
#     matrix_arma_t = .sarma.lag_matrix(t(Y_mat), t(R_mat),
#                                       lag_ar = model_order$ar[2],
#                                       lag_ma = model_order$ma[2])
#                     
#     # regression for SARMA model
# 
#     stats::lm(ar_0 ~ . + 0, data = matrix_arma_t)
#     
#     # fill ar and ma matrices
#     # updated to lag-order (phi_00/psi_00 in upper left corner)
#     # ar_mat is lhs of QARMA-equation (phi_00 = 1)
#     ar_mat = matrix(c(1, -arma_reg$coef[1:total_lag_ar]), byrow = TRUE,
#                     nrow = (model_order$ar[1] + 1), 
#                     ncol = (model_order$ar[2] + 1))
#     ma_mat = matrix(c(1, arma_reg$coef[(total_lag_ar + 1):
#                                          (total_lag_ar + total_lag_ma)]), byrow = TRUE,
#                     nrow = (model_order$ma[1] + 1),
#                     ncol = (model_order$ma[2] + 1))
#     if (total_lag_ar == 0)
#     {
#       ar_mat[1, 1] = 1
#     }
#   } else {
#     arma_reg = list(ar_aux = as.vector(ar_aux)[-1])
#     
#     # fill ar matrix (and ma matrix)
#     ar_mat = matrix(c(1, arma_reg$ar_aux[1:total_lag_ar]),
#                     nrow = (model_order$ar[1] + 1), ncol = (model_order$ar[2] + 1))
#     ma_mat = matrix(1)
#     
#     arma_reg$residuals = residuals_mat
#   }
#   
#   innov = arma_reg$residuals #residuals_mat
#   
#   improve = FALSE
#   if (improve == TRUE)
#   {
#     stdev = sqrt(sum(innov^2)/((nX - max_lag_x - model_order$ar[1]) * 
#                                  (nT - max_lag_t - model_order$ar[2])))
#     
#     arma_reg$coef = qarma.est_3rdstep(Y, ar_mat, ma_mat, model_order) +
#       arma_reg$coef
#     
#     ar_mat = matrix(c(1, -arma_reg$coef[1:total_lag_ar]), byrow = TRUE,
#                     nrow = (model_order$ar[1] + 1), 
#                     ncol = (model_order$ar[2] + 1))
#     ma_mat = matrix(c(1, arma_reg$coef[(total_lag_ar + 1):
#                                          (total_lag_ar + total_lag_ma)]), byrow = TRUE,
#                     nrow = (model_order$ma[1] + 1),
#                     ncol = (model_order$ma[2] + 1))
#   }
#   
#   # check stationarity
#   statTest = qarma.statTest(ar_mat)
#   if (!statTest)
#   {
#     warning("QARMA model not stationary, try another order for the AR-parts.")
#   }
#   
#   # preparation of output
#   rownames(ar_mat) = paste0("lag ", 0:model_order$ar[1])
#   colnames(ar_mat) = paste0("lag ", 0:model_order$ar[2])
#   rownames(ma_mat) = paste0("lag ", 0:model_order$ma[1])
#   colnames(ma_mat) = paste0("lag ", 0:model_order$ma[2])
#   
#   stdev = sqrt(sum(innov^2)/((nX - max_lag_x - model_order$ar[1]) * 
#                                (nT - max_lag_t - model_order$ar[2])))
#   coef_out = list(Y = Y, innov = innov, model = list(ar = ar_mat, ma = ma_mat,
#                                                      sigma = stdev), stnry = statTest)
#   class(coef_out) = "qarma"
#   attr(coef_out, "subclass") = "est"
#   return(coef_out)
# }
# 
# 
# .sarma.lag_matrix = function(Y, R, lag_ar, lag_ma)
# {
#   # computes lag matrices in one direction over rows (x-direction)
#   nX = dim(Y)[1]; nT = dim(Y)[2]
#   lag_n = lag_ar + lag_ma + 1
#   lag_max = max(lag_ar, lag_ma)
#   lag_obs = (nX - lag_max) * nT
#   
#   matrix_out = data.frame(matrix(NA, nrow = lag_obs, ncol = lag_n))
#   matrix_out[, 1] = as.vector(Y[(1 + lag_max):nX, ])
#   names(matrix_out)[1] = "ar_0"
#   
#   for (i in seq_len(lag_ar)) # use lag orders as loop indices
#   {
#     index_ar = (lag_max - i + 1):(nX - i)
#     matrix_out[, i + 1] = as.vector(Y[index_ar, ])
#     colnames(matrix_out)[i + 1] = paste0("ar_", i)
#   }
#   for (j in seq_len(lag_ma)) # use lag orders as loop indices
#   {
#     index_ma = (lag_max - j + 1):(nX - j)
#     matrix_out[, j + lag_ar + 1] = as.vector(R[index_ma, ])
#     colnames(matrix_out)[j + lag_max + 1] = paste0("ma_", j)
#   }
# 
#   return(matrix_out)
# }