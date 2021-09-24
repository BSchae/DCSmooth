
#---------------------------Model Order Selection------------------------------#

# Order selection by GPAC (see Illig/Truong-Van (2006))
qarma.order.gpac = function(Y, order_max = list(ar = c(1, 1), ma = c(1, 1)))
{
  # set up vectors for ar and ma
  ar_vec_test = expand.grid(0:order_max$ar[1], 0:order_max$ar[2])
  ar_vec_test = ar_vec_test[2:dim(ar_vec_test)[1], ]
  names(ar_vec_test) = NULL
  ar_vec = expand.grid(0:(order_max$ar[1] + 1), 0:(order_max$ar[2] + 1))
  ar_vec = ar_vec[2:dim(ar_vec)[1], ]
  ar_vec = ar_vec[order(ar_vec[, 1], ar_vec[, 2]), ]
  names(ar_vec) = NULL
  #ar_vec = rbind(ar_vec, ar_vec[dim(ar_vec)[1], ] + 1)
  ma_vec = expand.grid(0:(order_max$ma[1]), 0:(order_max$ma[2]))
  ma_vec = ma_vec[order(ma_vec[, 1], ma_vec[, 2]), ]
  names(ma_vec) = NULL
  #ma_vec = rbind(ma_vec, ma_vec[dim(ma_vec)[1], ] + 1)
  
  # set up matrix for results of yw estimation
  gpac_mat = matrix(NA, nrow = dim(ma_vec)[1], ncol = dim(ar_vec)[1])
  colnames(gpac_mat) = paste0("(", ar_vec[, 1], ", ", ar_vec[, 2], ")")
  rownames(gpac_mat) = paste0("(", ma_vec[, 1], ", ", ma_vec[, 2], ")")
  ar_test_names = paste0("(", ar_vec_test[, 1], ", ", ar_vec_test[, 2], ")")
  
  # calculate GPAC
  for (i in 1:dim(ar_vec)[1])
  {
    for (j in 1:dim(ma_vec)[1])
    {
      ar = ar_vec[i, ]
      ma = ma_vec[j, ]
      gpac_mat[j, i] = qarma.yw_matrix(Y, ar_order = as.numeric(ar), 
                                       ma_order = as.numeric(ma))[unlist(ar)[1] + 1,
                                                                  unlist(ar)[2] + 1]
    }
  }
  
  # find zeros of GPAC
  nRow = dim(gpac_mat)[1]; nCol = dim(gpac_mat)[2]
  gpac_count = gpac_mat*0 + 1
  
  # values chosen by hand, might be improved  
  gpac_count[which(abs(gpac_mat) < max(stats::quantile(abs(gpac_mat), 0.35), 0.1)
                   , arr.ind = TRUE)] = 0
  
  count_zeros = 1:nRow*0
  
  ar_out = vector(mode = "numeric")
  ma_out = vector(mode = "numeric")
  # score = vector(mode = "numeric")
  
  for (i in 1:nRow)
  {
    for (j in 1:(prod(order_max$ar + 1) - 1))
    {
      ar = as.vector(ar_vec_test[j, ], mode = "numeric")
      jCol = which(ar_test_names[j] == colnames(gpac_count))
      colIndex = which(as.logical(vapply(ar[1], '<=', 
                                         logical(length(ar_vec[, 1])), ar_vec[, 1]) *
                                    vapply(ar[2], '<=', 
                                           logical(length(ar_vec[, 2])), ar_vec[, 2])))
      gpac_count_sub = gpac_count[i, colIndex]
      if ((gpac_count_sub[1] == 1) & (sum(gpac_count_sub) == 1))
      {
        ar_out = rbind(ar_out, ar)
        ma_out = rbind(ma_out, ma = as.vector(ma_vec[i, ], mode = "numeric"))
        # score = rbind(score, sum(gpac_mat[i, colIndex]^2))
      }
      
    }
  }
  
  if (length(ar_out) != 0)
  {
    # ar_out = ar_out[which.min(score), ]
    # ma_out = ma_out[which.min(score), ]
    ar_out = ar_out[dim(ar_out)[1], ]
    ma_out = ma_out[dim(ma_out)[1], ]
  } else if (length(ar_out) == 0) {
    warning("No order selection by GPAC possible, iid. case is assumed")
    ar_out = c(0, 0)
    ma_out = c(0, 0)
  }
  
  if (length(ar_out) == 0) {
    warning("No order selection by GPAC possible, iid. case is assumed")
    ar_out = c(0, 0)
    ma_out = c(0, 0)
  }
  
  return(list(ar = ar_out, ma = ma_out))
}

# Order selection by BIC
qarma.order.bic = function(Y, order_max = list(ar = c(1, 1), ma = c(1, 1)))
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
      model_order = list(ar = as.numeric(ar_matrix[i, ]),
                         ma = as.numeric(ma_matrix[j, ]))
      qarma_model = qarma.est(Y, model_order = model_order)
      log_L = -n/2 * log(4*pi^2*qarma_model$model$sigma^2) -
        sum(qarma_model$innov^2)/(2*qarma_model$model$sigma^2)
      bic_matrix[i, j] = -2*log_L + sum(unlist(model_order)) * log(n)
    }
  }
  
  opt_index = which(bic_matrix == min(bic_matrix, na.rm = TRUE), arr.ind = TRUE)
  model_order_opt = list(ar = as.numeric(ar_matrix[opt_index[, 1], ]),
                         ma = as.numeric(ma_matrix[opt_index[, 2], ]))
  
  return(model_order_opt)
}

sarma.aic.bic <- function(Rmat, pmax = c(1, 1), qmax = c(1, 1), crit = "bic",
                      restr = NULL, sFUN = min) {
  
  if(crit == "bic") {
    crit.fun = stats::BIC
  }
  else if (crit == "aic") {
    crit.fun = stats::AIC
  }
  bic_x =  matrix(0, pmax[1] + 1, qmax[1] + 1)
  bic_t =  matrix(0, pmax[2] + 1, qmax[2] + 1)
  R_x = as.vector(Rmat)
  R_t = as.vector(t(Rmat))
  
  n.cores = parallel::detectCores(logical = TRUE) - 1
  #cl <- parallel::makeCluster(n.cores)
  doParallel::registerDoParallel(n.cores)
  `%dopar%` = foreach::`%dopar%`
  `%:%` = foreach::`%:%`
  
  bic_x = foreach::foreach(i = 1:(pmax[1] + 1), .combine = "rbind") %:%
    foreach::foreach(j = 1:(qmax[1] + 1), .combine = "c") %dopar% {
      bic = crit.fun(suppressWarnings(stats::arima(R_x,
                                                   order = c(i - 1, 0, j - 1), include.mean = FALSE)))
    }
  
  bic_t = foreach::foreach(i = 1:(pmax[2] + 1), .combine = "rbind") %:%
    foreach::foreach(j = 1:(qmax[2] + 1), .combine = "c") %dopar% {
      bic = crit.fun(suppressWarnings(stats::arima(R_t,
                                                   order = c(i - 1, 0, j - 1), include.mean = TRUE)))
    }
  
  doParallel::stopImplicitCluster()
  
  restr = substitute(restr)
  if(!is.null(restr)){
    ord.opt_x <- c(which(bic_x == sFUN(bic_x[eval(restr)]), arr.ind = TRUE) - 1)
    ord.opt_t <- c(which(bic_t == sFUN(bic_t[eval(restr)]), arr.ind = TRUE) - 1)
  } else {
    ord.opt_x <- c(which(bic_x == sFUN(bic_x), arr.ind = TRUE) - 1)
    ord.opt_t <- c(which(bic_t == sFUN(bic_t), arr.ind = TRUE) - 1)
  }
  # message("The optimal orders are:")
  # message("p_x = ", ord.opt_x[[1]])
  # message("p_t = ", ord.opt_t[[1]])
  # message("q_x = ", ord.opt_x[[2]])
  # message("q_t = ", ord.opt_t[[2]])
  names(ord.opt_x) = c("p", "q")
  names(ord.opt_t) = c("p", "q")
  # return(list(ord.opt_x = ord.opt_x,
  #             ord.opt_t = ord.opt_t))   
  
  return(list(ar = c(ord.opt_x[1], ord.opt_t[1]),
              ma = c(ord.opt_x[2], ord.opt_t[2])))
}