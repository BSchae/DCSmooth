################################################################################
#                                                                              #
#                     DCSmooth Package: SFARIMA Estimation                     #
#                                                                              #
################################################################################

#----------------------Formula for optimal bandwidths--------------------------#

sfarima.est = function(Y, model_order = list(ar = c(1, 1), ma = c(1, 1)))
{
  n_x = dim(Y)[1]; n_t = dim(Y)[2]
  
  # parameter estimation
  theta_init = rep(0, times = sum(unlist(model_order)) + 2)
  for (iterate in 1:3)
  {
    theta_opt  = stats::optim(theta_init, sfarima.rss, R_mat = Y, 
                              model_order = model_order, method = "BFGS")
    theta_init = theta_opt$par
  }
  
  # put coefficients into matrices
  d_vec = theta_opt$par[1:2]
  ar_x = c(1, -theta_opt$par[2 + seq_len(model_order$ar[1])])
  ar_t = c(1, -theta_opt$par[2 + model_order$ar[1] + 
                             seq_len(model_order$ar[2])])
  ma_x = c(1, theta_opt$par[2 + sum(model_order$ar) +
                            seq_len(model_order$ma[1])])
  ma_t = c(1, theta_opt$par[2 + sum(model_order$ar) + model_order$ma[1] + 
                            seq_len(model_order$ma[1])])
  
  # prepare results for output
  ar_mat = ar_x %*% t(ar_t)
  ma_mat = ma_x %*% t(ma_t)
  stdev = sqrt(theta_opt$value/(n_x * n_t))
  model = list(ar = ar_mat, ma = ma_mat, sigma = stdev)
  innov = sfarima.residuals(R_mat = Y, model = model)
  
  # check stationarity
  statTest = qarma.statTest(ar_mat)
  if (!statTest)
  {
    warning("QARMA model not stationary, try another order for the AR-parts.")
  }
  
  coef_out = list(Y = Y, innov = innov, model = list(ar = ar_mat, ma = ma_mat,
                                  d = d_vec, sigma = stdev), stnry = statTest)
  class(coef_out) = "qarma"
  attr(coef_out, "subclass") = "est"
  return(coef_out)
}

sfarima.rss = function(theta, R_mat,
                       model_order = list(ar = c(1, 1), ma = c(1, 1)))
{
  n_x = dim(R_mat)[1]; n_t = dim(R_mat)[2]
  k_x = min(50, n_x); k_t = min(100, n_t)
  
  # get coefficients from theta
  # theta = (d1, d2, ar_x, ar_t, ma_x, ma_t)
  d_vec = theta[1:2]
  ar_x = theta[2 + seq_len(model_order$ar[1])]
  ar_t = theta[2 + model_order$ar[1] + seq_len(model_order$ar[2])]
  ma_x = theta[2 + sum(model_order$ar) + seq_len(model_order$ma[1])]
  ma_t = theta[2 + sum(model_order$ar) + model_order$ma[1] + 
                 seq_len(model_order$ma[1])]
  
  # result matrices
  E_itm = R_mat * 0   # intermediate results
  E_fnl = R_mat * 0   # final results
  
  # calculate AR(inf) coefficients with long-memory
  ar_inf_x = c(1, astsa::ARMAtoAR(ar = ar_x, ma = ma_x, lag.max = k_x))
  d_x = choose(d_vec[1], 0:k_x) * ((-1)^(0:k_x))
  coef_x = cumsum_part_reverse(d_x, ar_inf_x)
  
  # a = arcoef(ar = ar_x, ma = ma_x, d = d_vec[1], k = 100)
  
  ar_inf_t = t(c(1, astsa::ARMAtoAR(ar = ar_t, ma = ma_t, lag.max = k_t)))
  d_t = choose(d_vec[2], 0:k_t) * ((-1)^(0:k_t))
  coef_t = cumsum_part_reverse(d_t, ar_inf_t)
  
  for (j in 1:n_t)
  {
    E_itm[, j] = R_mat[, j:max(1, j - k_t + 1), drop = FALSE] %*%
      coef_t[1:min(j, k_t), drop = FALSE]
  }
  
  for (i in 1:n_x)
  {
    E_fnl[i, ] = coef_x[1:min(i, k_x), drop = FALSE] %*%
      E_itm[i:max(1, i - k_x + 1), , drop = FALSE]
  }
  
  RSS = sum(E_fnl^2)
  
  return(RSS)
}

sfarima.residuals = function(R_mat, model)
{
  n_x = dim(R_mat)[1]; n_t = dim(R_mat)[2]
  k_x = min(50, n_x); k_t = min(100, n_t)
  
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