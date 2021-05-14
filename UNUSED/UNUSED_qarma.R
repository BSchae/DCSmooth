#--------------------Unused Estimation Functions for QARMA---------------------#

qarma.est2 = function(Y, model_order = list(ar = c(1, 1), ma = c(1, 1)))
{
  nX = dim(Y)[1]; nT = dim(Y)[2]
  
  # calculate ARMA-orders for different purposes
  max_lag_x = max(model_order$ar[1], model_order$ma[1])
  max_lag_t = max(model_order$ar[2], model_order$ma[2])
  total_lag_ar = (model_order$ar[1] + 1) * (model_order$ar[2] + 1) - 1
  total_lag_ma = (model_order$ma[1] + 1) * (model_order$ma[2] + 1) - 1
  max_lag_ar = max(model_order$ar)
  max_lag_ma = max(model_order$ma)
  
  ### AR-only auxiliary estimation model by YW-estimation ###
  m1_ar = max(2*max_lag_x, 2)
  m2_ar = max(2*max_lag_t, 2)
  ar_aux = qarma.yw_matrix(Y, ar_order = c(m1_ar, m2_ar))
  
  # calculate residuals from auxiliary AR model
  residuals_mat = matrix(0, nrow = nX - m1_ar, ncol = nT - m2_ar)
  for (i in 1:(nX - m1_ar))
  {
    for (j in 1:(nT - m2_ar))
    {
      residuals_mat[i, j] = -sum(ar_aux * Y[i:(i + m1_ar), j:(j + m2_ar)])
    }
  }
  
  # matrix_aux = .qarma.fill_lag_matrix(Y, lag_x = model_order$ar[1],
  #                                     lag_t = model_order$ar[2], include_00 = TRUE,
  #                                     name_prefix = "ar")
  # aux_reg = lm(ar_00 ~ . + 0, data = matrix_aux)
  # residuals_aux = matrix(aux_reg$residuals, nrow = nX - model_order$ar[1],
  #                        ncol = nT - model_order$ar[2])
  
  if (total_lag_ma > 0)
  {
    # set up submatrices of Y, res for use in fill_matrix procedure
    # Y_mat = Y[(max_lag_x + 1):nX, (max_lag_t + 1):nT]
    # residuals_mat = residuals_aux[(max_lag_x - model_order$ma[1] + 1):
    #                                 (nX - model_order$ar[1]),
    #                               (max_lag_t - model_order$ma[2] + 1):
    #                                 (nT - model_order$ar[2])]
    
    Y_mat = Y[(m1_ar + model_order$ma[1] - model_order$ar[1] + 1):nX,
              (m2_ar + model_order$ma[2] - model_order$ar[2] + 1):nT]
    
    # set up matrix for complete arma (dimension (nX - ar_x)*(nT - ar_t))
    matrix_arma = cbind(.qarma.fill_lag_matrix(Y_mat, lag_x = model_order$ar[1],
                                               lag_t = model_order$ar[2], include_00 = TRUE,
                                               name_prefix = "ar"),
                        .qarma.fill_lag_matrix(residuals_mat, lag_x = model_order$ma[1],
                                               lag_t = model_order$ma[2], include_00 = FALSE,
                                               name_prefix = "ma"))
    
    # regression for QARMA model
    arma_reg = lm(ar_00 ~ . + 0, data = matrix_arma)
    
    # fill ar and ma matrices
    # updated to lag-order (phi_00/psi_00 in upper left corner)
    # ar_mat is lhs of QARMA-equation (phi_00 = 1)
    ar_mat = matrix(c(1, -arma_reg$coef[1:total_lag_ar]), byrow = TRUE,
                    nrow = (model_order$ar[1] + 1),
                    ncol = (model_order$ar[2] + 1))
    ma_mat = matrix(c(1, arma_reg$coef[(total_lag_ar + 1):
                                         (total_lag_ar + total_lag_ma)]), byrow = TRUE,
                    nrow = (model_order$ma[1] + 1),
                    ncol = (model_order$ma[2] + 1))
    if (total_lag_ar == 0)
    {
      ar_mat[1, 1] = 1
    }
  } else {
    arma_reg = aux_reg
    
    # fill ar matrix (and ma matrix)
    ar_mat = matrix(c(1, -aux_reg$coef[1:total_lag_ar]), byrow = TRUE,
                    nrow = (model_order$ar[1] + 1), ncol = (model_order$ar[2] + 1))
    ma_mat = matrix(1)
  }
  
  innov = arma_reg$residuals #residuals_mat
  
  # check stationarity
  statTest = qarma.statTest(ar_mat)
  if (statTest == FALSE)
  {
    warning("QARMA model not stationary, try another order for the AR-parts.")
  }
  
  # preparation of output
  stdev = sqrt(sum(innov^2)/((nX - max_lag_x - model_order$ar[1]) *
                               (nT - max_lag_t - model_order$ar[2])))
  coef_out = list(ar = ar_mat, ma = ma_mat, sigma = stdev, innov = innov,
                  stationary = statTest)
  return(coef_out)
}
# 
# qarma.est2 = function(Y, model_order = list(ar = c(1, 1), ma = c(1, 1)))
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
# # 1st step: Yule-Walker estimation of AR(2*p_1, 2*p_2) model
#   m1_ar = max(1*max_lag_x, 1)
#   m2_ar = max(1*max_lag_t, 1)
#   phi_ar_matrix = qarma.yw_matrix(Y, ar_order = c(m1_ar, m2_ar))
#   
#   # calculate residuals from auxiliary AR model
#   ar_residuals_matrix = matrix(0, nrow = nX - m1_ar, ncol = nT - m2_ar)
#   for (i in 1:(nX - m1_ar))
#   {
#     for (j in 1:(nT - m2_ar))
#     {
#       ar_residuals_matrix[i, j] = -sum(phi_ar_matrix *
#                                          Y[i:(i + m1_ar), j:(j + m2_ar)])
#     }
#   }
#   
# # 2nd step: estimate ARMA model by linear regression (only for ma-orders > 0)
#   arma_lm_matrix = data.frame(matrix(NA, 
#                       nrow = (nX - m1_ar - max_lag_x) *
#                       (nT - m2_ar - max_lag_t), 
#                       ncol = total_lag_ar + total_lag_ma + 1))
#   
#   # fill arma_lm_matrix with Y-values and computed residuals
#   for (i in 0:model_order$ar[1])
#   {
#     for (j in 0:model_order$ar[2])
#     {
#       index_x = (m1_ar + max_lag_x + 1 - i):(nX - i)
#       index_t = (m2_ar + max_lag_t + 1 - j):(nT - j)
#       # over colums first (inner loop)
#       arma_lm_matrix[, (j + 1) + i * (model_order$ar[2] + 1)] =
#         as.vector(Y[index_x, index_t])
#       names(arma_lm_matrix)[(j + 1) + i * (model_order$ar[2] + 1)] =
#         paste0("lag_ar_", i, j)
#     }
#   }
#   for (i in 0:model_order$ma[1]) # use true lag orders as loop indices
#   {
#     for (j in 0:model_order$ma[2])
#     {
#       if (!(i == 0 && j == 0))  # lag_ma_00 is not needed (no coefficient)
#       {
#         index_x = (max_lag_x + 1 - i):(nX - i - m1_ar)
#         index_t = (max_lag_t + 1 - j):(nT - j - m2_ar)
#         arma_lm_matrix[, total_lag_ar + (j + 1) + i * (model_order$ma[2] + 1)] =
#           as.vector(ar_residuals_matrix[index_x, index_t])
#         names(arma_lm_matrix)[total_lag_ar + (j + 1) + i *
#                               (model_order$ma[2] + 1)] = paste0("lag_ma_", i, j)
#       }
#     }
#   }
#   
#   # linear regression & prepare output
#   arma_lm = lm(lag_ar_00 ~ . + 0, data = arma_lm_matrix)
#   
#   improve = TRUE
#   if (improve == TRUE)
#   {
#     # byrow = TRUE and reverse needed as lm order is 01,10,11...
#     ar_mat = matrix(c(arma_lm$coef[total_lag_ar:1], -1), byrow = TRUE,
#                     nrow = (model_order$ar[1] + 1), ncol = (model_order$ar[2] + 1))
#     ar_mat[model_order$ar[1] + 1, model_order$ar[2] + 1] = -1
#     ma_mat = matrix(c(arma_lm$coef[(total_lag_ar + total_lag_ma):
#                                     (total_lag_ar + 1)], 1), byrow = TRUE,
#                     nrow = (model_order$ma[1] + 1), ncol =
#                     (model_order$ma[2] + 1))
#     ma_mat[model_order$ma[1] + 1, model_order$ma[2] + 1] = 1
#     innov = matrix(arma_lm$residuals, nrow = nX - m1_ar - max_lag_x,
#                   ncol = nT - m2_ar - max_lag_t)
#     stdev = sqrt(sum(innov^2)/((nX - m1_ar - model_order$ma[1]) * 
#                                         (nX - m1_ar - model_order$ma[1])))
#     
#     arma_lm$coef = qarma.est.3rdstep(Y, ar_mat, ma_mat, model_order) + arma_lm$coef
#   }
#   
#   # byrow = TRUE and reverse needed as lm order is 01,10,11...
#   ar_mat = matrix(c(arma_lm$coef[total_lag_ar:1], -1), byrow = TRUE,
#                   nrow = (model_order$ar[1] + 1), ncol = (model_order$ar[2] + 1))
#   ar_mat[model_order$ar[1] + 1, model_order$ar[2] + 1] = -1
#   ma_mat = matrix(c(arma_lm$coef[(total_lag_ar + total_lag_ma):
#                                    (total_lag_ar + 1)], 1), byrow = TRUE,
#                   nrow = (model_order$ma[1] + 1), ncol =
#                     (model_order$ma[2] + 1))
#   ma_mat[model_order$ma[1] + 1, model_order$ma[2] + 1] = 1
#   innov = matrix(arma_lm$residuals, nrow = nX - m1_ar - max_lag_x,
#                  ncol = nT - m2_ar - max_lag_t)
#   stdev = sqrt(sum(innov^2)/((nX - m1_ar - model_order$ma[1]) * 
#                                (nX - m1_ar - model_order$ma[1])))
#   # end 3rd step
#   
#   # check stationarity
#   statTest = qarma.statTest(ar_mat)
#   if (statTest == FALSE)
#   {
#     warning("QARMA model not stationary, try another order for the AR-parts.")
#     error_out = list(ar = ar_mat, ma = ma_mat, sigma = stdev, innov = innov)
#     return(error_out)
#   }
#   
#   # preparation of output
#   coef_out = list(ar = ar_mat, ma = ma_mat, sigma = stdev, innov = innov)
#   return(coef_out)
# }
# 
qarma.est.3rdstep2 = function(Y, ar_mat, ma_mat, model_order)
{
  nX = dim(Y)[1]; nT = dim(Y)[2]
  
  # value (x,t) of ma_mat has to be zero
  ma_mat_0 = ma_mat; ma_mat_0[model_order$ma[1] + 1, model_order$ma[2] + 1] = 0
  ar_mat_0 = ar_mat; ar_mat_0[model_order$ar[1] + 1, model_order$ar[2] + 1] = 0
  
  # calculate ARMA-orders for different purposes
  max_lag_x = max(model_order$ar[1], model_order$ma[1])
  max_lag_t = max(model_order$ar[2], model_order$ma[2])
  total_lag_ar = (model_order$ar[1] + 1) * (model_order$ar[2] + 1) - 1
  total_lag_ma = (model_order$ma[1] + 1) * (model_order$ma[2] + 1) - 1
  #max_lag_ar = max(model_order$ar)
  #max_lag_ma = max(model_order$ma)
  
  z_matrix = matrix(0, nrow = nX, ncol = nT)
  v_matrix = w_matrix = z_matrix
  
  for (i in (max_lag_x + 1):nX)
  {
    for (j in (max_lag_t + 1):nT)
    {
      ind_ar_x = (i - model_order$ar[1]):i
      ind_ar_t = (j - model_order$ar[2]):j
      ind_ma_x = (i - model_order$ma[1]):i
      ind_ma_t = (j - model_order$ma[2]):j
      
      # "-" below is due to definition of ar_mat
      z_matrix[i, j] = sum(-ar_mat * Y[ind_ar_x, ind_ar_t]) -
        sum(ma_mat_0 * z_matrix[ind_ma_x, ind_ma_t])
      v_matrix[i, j] = sum(ar_mat_0 * v_matrix[ind_ar_x, ind_ar_t]) +
        z_matrix[i, j]
      w_matrix[i, j] = sum(-ma_mat_0 * w_matrix[ind_ma_x, ind_ma_t]) +
        z_matrix[i, j]
    }
  }
  
  zvw_lm_matrix = data.frame(matrix(NA,
                                    nrow = (nX - max_lag_x) * (nT - max_lag_t),
                                    ncol = total_lag_ar + total_lag_ma + 1))
  zvw_lm_matrix[, 1] =
    as.vector(z_matrix[(max_lag_x + 1):nX, (max_lag_t + 1):nT])
  names(zvw_lm_matrix)[1] = "Z"
  
  for (i in 0:model_order$ar[1])
  {
    for (j in 0:model_order$ar[2])
    {
      if (!(i == 0 && j == 0))
      {
        zvw_lm_matrix[, i*(model_order$ar[2] + 1) + j + 1] =
          as.vector(v_matrix[(max_lag_x + 1):nX - i, (max_lag_t + 1):nT - j])
        names(zvw_lm_matrix)[i*(model_order$ar[2] + 1) + j + 1] =
          paste0("V_", i, j)
      }
    }
  }
  for (i in 0:model_order$ma[1])
  {
    for (j in 0:model_order$ma[2])
    {
      if (!(i == 0 && j == 0))
      {
        zvw_lm_matrix[, total_lag_ar + i*(model_order$ma[2] + 1) + j + 1] =
          as.vector(w_matrix[(max_lag_x + 1):nX - i, (max_lag_t + 1):nT - j])
        names(zvw_lm_matrix)[total_lag_ar + i*(model_order$ma[2] + 1) + j + 1] =
          paste0("W_", i, j)
      }
    }
  }
  
  # linear regression
  zvw_lm = lm(Z ~ . + 0, data = zvw_lm_matrix)
  
  return(zvw_lm$coef)
}