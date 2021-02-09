################################################################################
#                                                                              #
#                DCSmooth Package: estimation of QARMA for cf                  #
#                                                                              #
################################################################################

# This file includes all functions related to the QARMA estimation for the cf
# coefficient of the bandwidth selection procedure


### Part A: Estimation Functions

#-----------------------Calculation of cf coefficient--------------------------#

qarma.cf = function(Y, model_order = list(ar = c(1, 1), ma = c(1, 1)))
{
  # estimation of qarma model
  qarma_model = qarma.est(Y, model_order = model_order)

  # get fractions for spectral density
  cf = sum(qarma_model$ma)^2/sum(qarma_model$ar)^2 * qarma_model$sigma^2
  return_list = list(cf = cf, qarma_model = qarma_model)
  
  return(return_list)
}

#-------------------------------Model Selection--------------------------------#

qarma.order = function(Y, order_max = list(ar = c(1, 1), ma = c(1, 1)))
{
  n = prod(dim(Y))
  ar_matrix = expand.grid(0:order_max$ar[1], 0:order_max$ar[2])
  names(ar_matrix) = NULL
  ma_matrix = expand.grid(0:order_max$ma[1], 0:order_max$ma[2])
  names(ma_matrix) = NULL
  bic_matrix = matrix(NA, nrow = dim(ar_matrix)[1], ncol = dim(ma_matrix)[1])
  
  for (i in 1:dim(ar_matrix)[1])
  {
    for (j in 1:dim(ma_matrix)[1])
    {
      print(c(i, j))
      model_order = list(ar = as.numeric(ar_matrix[i, ]),
                         ma = as.numeric(ma_matrix[j, ]))
      qarma_model = qarma.est(Y, model_order = model_order)
      log_L = -n/2 * log(2*pi*qarma_model$sigma^2) -
              sum(qarma_model$innov^2)/(2*qarma_model$sigma^2)
      bic_matrix[i, j] = -2*log_L + sum(unlist(model_order)) * log(n)
    }
  }
  
  opt_index = which(bic_matrix == min(bic_matrix, na.rm = TRUE), arr.ind = TRUE)
  model_order_opt = list(ar = as.numeric(ar_matrix[opt_index[, 1], ]),
                         ma = as.numeric(ma_matrix[opt_index[, 2], ]))
  
  return(model_order_opt)
}

#------------------------Estimation Function for QARMA-------------------------#

qarma.est = function(Y, model_order = list(ar = c(1, 1), ma = c(1, 1)))
{
  nX = dim(Y)[1]; nT = dim(Y)[2]
  
  # calculate ARMA-orders for different purposes
  max_lag_x = max(model_order$ar[1], model_order$ma[1])
  max_lag_t = max(model_order$ar[2], model_order$ma[2])
  total_lag_ar = (model_order$ar[1] + 1) * (model_order$ar[2] + 1) - 1
  total_lag_ma = (model_order$ma[1] + 1) * (model_order$ma[2] + 1) - 1
  max_lag_ar = max(model_order$ar)
  max_lag_ma = max(model_order$ma)
  
  ### AR-only auxiliary estimation model ###
  
  matrix_aux = .qarma.fill_lag_matrix(Y, lag_x = model_order$ar[1],
                                      lag_t = model_order$ar[2], include_00 = TRUE,
                                      name_prefix = "ar")
  aux_reg = lm(ar_00 ~ . + 0, data = matrix_aux)
  residuals_aux = matrix(aux_reg$residuals, nrow = nX - model_order$ar[1], 
                         ncol = nT - model_order$ar[2])
  
  if (total_lag_ma > 0)
  {
    # set up submatrices of Y, res for use in fill_matrix procedure
    Y_mat = Y[(max_lag_x + 1):nX, (max_lag_t + 1):nT]
    residuals_mat = residuals_aux[(max_lag_x - model_order$ma[1] + 1):
                                    (nX - model_order$ar[1]),
                                  (max_lag_t - model_order$ma[2] + 1):
                                    (nT - model_order$ar[2])]
    
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
    ar_mat = matrix(c(arma_reg$coef[total_lag_ar:1], -1), byrow = TRUE,
                    nrow = (model_order$ar[1] + 1), 
                    ncol = (model_order$ar[2] + 1))
    ma_mat = matrix(c(arma_reg$coef[(total_lag_ar + total_lag_ma):
                                      (total_lag_ar + 1)], 1), byrow = TRUE,
                    nrow = (model_order$ma[1] + 1),
                    ncol = (model_order$ma[2] + 1))
  } else {
    arma_reg = aux_reg
    
    # fill ar matrix (and ma matrix)
    ar_mat = matrix(c(aux_reg$coef[total_lag_ar:1], -1), byrow = TRUE,
                    nrow = (model_order$ar[1] + 1), ncol = (model_order$ar[2] + 1))
    ma_mat = matrix(1)
  }
  
  innov = arma_reg$residuals #residuals_mat
  
  # check stationarity
  statTest = qarma.statTest(ar_mat)
  if (statTest == FALSE)
  {
    stop("QARMA model not stationary, try another order for the AR-parts.")
  }
  
  # preparation of output
  stdev = sqrt(sum(innov^2)/((nX - max_lag_x - model_order$ar[1]) * 
                               (nT - max_lag_t - model_order$ar[2])))
  coef_out = list(ar = ar_mat, ma = ma_mat, sigma = stdev, innov = innov)
  return(coef_out)
}

qarma.est1 = function(Y, model_order = list(ar = c(1, 1), ma = c(1, 1)))
{
  nX = dim(Y)[1]; nT = dim(Y)[2]
  
  # calculate ARMA-orders for different purposes
  max_lag_x = max(model_order$ar[1], model_order$ma[1])
  max_lag_t = max(model_order$ar[2], model_order$ma[2])
  total_lag_ar = (model_order$ar[1] + 1) * (model_order$ar[2] + 1) - 1
  total_lag_ma = (model_order$ma[1] + 1) * (model_order$ma[2] + 1) - 1
  max_lag_ar = max(model_order$ar)
  max_lag_ma = max(model_order$ma)
  
  ### AR-only auxiliary estimation model
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
  
  ### iteration loop (number of iterations might increase the precision)
  if (total_lag_ma > 0)
  {
    for (loop in 1:3)
    {
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
      
      # build nX*nT matrix of residuals
      residuals_arma[(max_lag_x + 1):nX, (max_lag_t + 1):nT] =
        matrix(aux_reg$residuals, nrow = (nX - max_lag_x),
               ncol = (nT - max_lag_t))
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
  statTest = qarma.statTest(ar_mat)
  if (statTest == FALSE)
  {
    stop("QARMA model not stationary, try another order for the AR-parts.")
  }
  
  # preparation of output
  stdev = sd(innov)
  coef_out = list(ar = ar_mat, ma = ma_mat, sigma = stdev, innov = innov)
  return(coef_out)
}

# Estimation function using the Hannan-Rissanen Algorithm

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

# 1st step: Yule-Walker estimation of AR(2*p_1, 2*p_2) model
  m1_ar = max(2*max_lag_x, 1)
  m2_ar = max(2*max_lag_t, 1)
  phi_ar_matrix = qarma.yw_matrix(Y, ar_order = c(m1_ar, m2_ar))
  
  # calculate residuals from auxiliary AR model
  ar_residuals_matrix = matrix(0, nrow = nX - m1_ar, ncol = nT - m2_ar)
  for (i in 1:(nX - m1_ar))
  {
    for (j in 1:(nT - m2_ar))
    {
      ar_residuals_matrix[i, j] = -sum(phi_ar_matrix *
                                         Y[i:(i + m1_ar), j:(j + m2_ar)])
    }
  }
  
# 2nd step: estimate ARMA model by linear regression (only for ma-orders > 0)
  arma_lm_matrix = data.frame(matrix(NA, 
                      nrow = (nX - m1_ar - max_lag_x) *
                      (nT - m2_ar - max_lag_t), 
                      ncol = total_lag_ar + total_lag_ma + 1))
  
  # fill arma_lm_matrix with Y-values and computed residuals
  for (i in 0:model_order$ar[1])
  {
    for (j in 0:model_order$ar[2])
    {
      index_x = (m1_ar + max_lag_x + 1 - i):(nX - i)
      index_t = (m2_ar + max_lag_t + 1 - j):(nT - j)
      # over colums first (inner loop)
      arma_lm_matrix[, (j + 1) + i * (model_order$ar[2] + 1)] =
        as.vector(Y[index_x, index_t])
      names(arma_lm_matrix)[(j + 1) + i * (model_order$ar[2] + 1)] =
        paste0("lag_ar_", i, j)
    }
  }
  for (i in 0:model_order$ma[1]) # use true lag orders as loop indices
  {
    for (j in 0:model_order$ma[2])
    {
      if (!(i == 0 && j == 0))  # lag_ma_00 is not needed (no coefficient)
      {
        index_x = (max_lag_x + 1 - i):(nX - i - m1_ar)
        index_t = (max_lag_t + 1 - j):(nT - j - m2_ar)
        arma_lm_matrix[, total_lag_ar + (j + 1) + i * (model_order$ma[2] + 1)] =
          as.vector(ar_residuals_matrix[index_x, index_t])
        names(arma_lm_matrix)[total_lag_ar + (j + 1) + i *
                              (model_order$ma[2] + 1)] = paste0("lag_ma_", i, j)
      }
    }
  }
  
  # linear regression & prepare output
  arma_lm = lm(lag_ar_00 ~ . + 0, data = arma_lm_matrix)
  
  improve = TRUE
  if (improve == TRUE)
  {
    # byrow = TRUE and reverse needed as lm order is 01,10,11...
    ar_mat = matrix(c(arma_lm$coef[total_lag_ar:1], -1), byrow = TRUE,
                    nrow = (model_order$ar[1] + 1), ncol = (model_order$ar[2] + 1))
    ar_mat[model_order$ar[1] + 1, model_order$ar[2] + 1] = -1
    ma_mat = matrix(c(arma_lm$coef[(total_lag_ar + total_lag_ma):
                                    (total_lag_ar + 1)], 1), byrow = TRUE,
                    nrow = (model_order$ma[1] + 1), ncol =
                      (model_order$ma[2] + 1))
    ma_mat[model_order$ma[1] + 1, model_order$ma[2] + 1] = 1
    innov = matrix(arma_lm$residuals, nrow = nX - m1_ar - max_lag_x,
                  ncol = nT - m2_ar - max_lag_t)
    stdev = sqrt(sum(innov^2)/((nX - m1_ar - model_order$ma[1]) * 
                                        (nX - m1_ar - model_order$ma[1])))
    
    arma_lm$coef = qarma.est.3rdstep(Y, ar_mat, ma_mat, model_order) + arma_lm$coef
  }
  
  # byrow = TRUE and reverse needed as lm order is 01,10,11...
  ar_mat = matrix(c(arma_lm$coef[total_lag_ar:1], -1), byrow = TRUE,
                  nrow = (model_order$ar[1] + 1), ncol = (model_order$ar[2] + 1))
  ar_mat[model_order$ar[1] + 1, model_order$ar[2] + 1] = -1
  ma_mat = matrix(c(arma_lm$coef[(total_lag_ar + total_lag_ma):
                                   (total_lag_ar + 1)], 1), byrow = TRUE,
                  nrow = (model_order$ma[1] + 1), ncol =
                    (model_order$ma[2] + 1))
  ma_mat[model_order$ma[1] + 1, model_order$ma[2] + 1] = 1
  innov = matrix(arma_lm$residuals, nrow = nX - m1_ar - max_lag_x,
                 ncol = nT - m2_ar - max_lag_t)
  stdev = sqrt(sum(innov^2)/((nX - m1_ar - model_order$ma[1]) * 
                               (nX - m1_ar - model_order$ma[1])))
  # end 3rd step
  
  # check stationarity
  statTest = qarma.statTest(ar_mat)
  if (statTest == FALSE)
  {
    warning("QARMA model not stationary, try another order for the AR-parts.")
    error_out = list(ar = ar_mat, ma = ma_mat, sigma = stdev, innov = innov)
    return(error_out)
  }
  
  # preparation of output
  coef_out = list(ar = ar_mat, ma = ma_mat, sigma = stdev, innov = innov)
  return(coef_out)
}

### 3rd step of Hannen-Rissanen algorithm, does not seem to provide improved
### results (probably not correct implemnted)

qarma.est.3rdstep = function(Y, ar_mat, ma_mat, model_order)
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

#-----------------------------Auxiliary Functions------------------------------#

.qarma.fill_lag_matrix = function(Y, lag_x, lag_t, include_00 = TRUE, name_prefix)
{
  nX = dim(Y)[1]; nT = dim(Y)[2]
  lag_obs = (nX - lag_x) * (nT - lag_t)
  lag_n = (lag_x + 1) * (lag_t + 1)
  
  matrix_out = data.frame(matrix(NA, nrow = lag_obs, ncol = lag_n))
  
  for (i in 0:lag_x) # use lag orders as loop indices
  {
    for (j in 0:lag_t)
    {
      index_x = (lag_x - i + 1):(nX - i)
      index_t = (lag_t - j + 1):(nT - j)
      matrix_out[, (j + 1) + i * (lag_t + 1)] = as.vector(Y[index_x, index_t])
      colnames(matrix_out)[(j + 1) + i * (lag_t + 1)] = 
        paste0(name_prefix, "_", i, j)
    }
  }
  
  if (include_00 == FALSE && lag_n > 1)
  {
    matrix_out = matrix_out[2:lag_n]
  } else if (include_00 == FALSE && lag_n == 1) {
    matrix_out = matrix_out * 0
  }
  
  return(matrix_out)
}

#-------------------Yule-Walker-Estimation Function for AR---------------------#

qarma.yw_matrix = function(Y, ar_order = c(1, 1))
{
  # set up parameters, result matrices and vectors etc.
  nX = dim(Y)[1]; nT = dim(Y)[2]; n = nX * nT
  p1_ar = ar_order[1]
  p2_ar = ar_order[2]
  d_ar = (p1_ar + 1) * (p2_ar + 1) - 1

  acf_matrix_ar = matrix(NA, nrow = d_ar, ncol = d_ar)
  acf_vector_ar = vector(mode = "numeric", length = d_ar)
  # st_mat is over rows first (i, j) -> (i - 1, j), which is normal R order
  st_mat = expand.grid(p1 = c(p1_ar:0), p2 = c(p2_ar:0))

  # calculate acf_matrix and acf_vector
  for (k in 1:d_ar)
  {
    for (l in 1:d_ar)
    {
      st_vec = st_mat[k, ] - st_mat[l, ]
      acf_matrix_ar[k, l] = .matrix_acf(Y, st_vec)
    }
    acf_vector_ar[k] = .matrix_acf(Y, st_mat[k, ])
  }
  
  # solving Yule-Walker equations
  yw.estimators = solve(acf_matrix_ar) %*% acf_vector_ar
  phi_ar_matrix = matrix(c(yw.estimators, -1), nrow = (p1_ar + 1),
                         ncol = (p2_ar + 1))
  
  return(phi_ar_matrix)
}

### auxiliary functions
.matrix_acf = function(Y, st_vec)
{
  s = as.numeric(st_vec[1])
  t = as.numeric(st_vec[2])
  n = dim(Y)
  
  if ((s >= 0 && t >= 0) || (s < 0 && t < 0))
  {
    acf_out = .matrix_acf_positive(Y, abs(s), abs(t), n)
  } else if (s < 0 && t >= 0) {
    acf_out = .matrix_acf_negative(Y, s = -s, t, n)
  } else {
    acf_out = .matrix_acf_negative(Y, s, t = -t, n)
  }
  
  # unbiased estimator from Ha/Newton(1993)
  unbiased_factor = (n[1] - s)*(n[2] - t)/prod(n)
  return(acf_out/unbiased_factor)
}

.matrix_acf_positive = function(Y, s, t, n)
{
  Y0 = Y[1:(n[1] - s), 1:(n[2] - t)]
  Y1 = Y[(s + 1):n[1], (t + 1):n[2]]
  acf_out = sum(Y0 * Y1)
  return(acf_out)
}

.matrix_acf_negative = function(Y, s, t, n)
{
  Y0 = Y[1:(n[1] - s), (t + 1):n[2]]
  Y1 = Y[(s + 1):n[1], 1:(n[2] - t)]
  acf_out = sum(Y0 * Y1)
  return(acf_out)
}

#----------------------------Simulation Function------------------------------#

qarma.sim = function(nX, nT, model = list(ar, ma, sigma))
{
  ar_mat = as.matrix(model$ar); ma_mat = as.matrix(model$ma)
  ar_x = dim(ar_mat)[1] - 1; ar_t = dim(ar_mat)[2] - 1
  ma_x = dim(ma_mat)[1] - 1; ma_t = dim(ma_mat)[2] - 1
  xInit = max(ar_x, ma_x) + 1
  tInit = max(ar_t, ma_t) + 1
  
  # check matrices alpha, beta
  ma_mat[ma_x + 1, ma_t + 1] = 1 # MA-coefficients
  ar_mat[ar_x + 1, ar_t + 1] = 0  # AR-coefficients
  
  nMat = floor(1.25 * c(nX, nT))
  errorMat = matrix(rnorm(prod(nMat)), nrow = nMat[1], ncol = nMat[2]) *
    model$sigma
  armaMat = errorMat
  
  for (i in xInit:nMat[1])
  {
    for (j in tInit:nMat[2])
    {
      armaMat[i, j] = sum(ar_mat * armaMat[(i - ar_x):i, (j - ar_t):j]) +
                      sum(ma_mat * errorMat[(i - ma_x):i, (j - ma_t):j])
    }
  }
  
  armaOut = armaMat[(nMat[1] - nX + 1):nMat[1], (nMat[2] - nT + 1):nMat[2]]
  errorOut = errorMat[(nMat[1] - nX + 1):nMat[1], (nMat[2] - nT + 1):nMat[2]]
  
  listOut = list(Y = armaOut, innov = errorOut,
                  ar = ar_mat, ma = ma_mat, sigma = model$sigma)
  return(listOut)
}

#-------------------------Test for QARMA stationarity-------------------------#

qarma.statTest = function(ar)
{
  ar[dim(ar)[1], dim(ar)[2]] = -1
  outValue = TRUE
  
  # compute reference sign at center (0, 0)
  signRef = sign(qarma.auxF(0, 0, ar))
  
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
      
      signTest = sign(Re(qarma.auxF(z1c, z2c, ar)))
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

qarma.auxF = function(z1c, z2c, ar)
{
  ar_x = dim(ar)[1] - 1; ar_t = dim(ar)[2] - 1
  
  # set up vectors for z^i
  zx_vec = z1c^(ar_x:0)
  zt_vec = z2c^(ar_t:0)
  
  out = t(zx_vec) %*% ar %*% zt_vec
  
  return(out)
}

qarma.SSDE = function(Y, ar, ma, stdev)
{
  X = T = seq(from = -pi, to = pi, length.out = 100)
  ySpectrum = matrix(NA, nrow = 100, ncol = 100)
  
  for (i in 1:100)
  {
    for (j in 1:100)
    {
      omega = c(X[i], T[j])
      ySpectrum[i, j] = qarma.spec(Y, ar, ma, stdev, omega)
    }
  }
  return(ySpectrum)
}

#----------------------Calculation of spectral density-------------------------#

qarma.spec = function(Y, ar, ma, stdev, omega)
{
  lagAR = dim(ar) - 1; lagMA = dim(ma) - 1
  z1 = complex(argument = omega[1]); z2 = complex(argument = omega[2])
  
  g00 = sd(Y)^2
  
  arSum = Re(sum(ar * (z1^(lagAR[1]:0)) %*% t(z2^(lagAR[2]:0))) * 
                 sum(ar * ((1/z1)^(lagAR[1]:0)) %*% t((1/z2)^(lagAR[2]:0))))
  maSum = Re(sum(ma * (z1^(lagMA[1]:0)) %*% t(z2^(lagMA[2]:0))) * 
                  sum(ma * ((1/z1)^(lagMA[1]:0)) %*% t((1/z2)^(lagMA[2]:0))))
  
  specDensOut = stdev^2/g00 * maSum/arSum # check this again ???
  
  return(specDensOut)
}