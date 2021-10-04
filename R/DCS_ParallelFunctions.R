LP.const0 = function(opts_list, Y)
{
  LP_dcs_const0_BMod(yMat = Y,
                     hVec = opts_list$h,
                     polyOrderVec = opts_list$p,
                     drvVec = opts_list$drv,
                     muVec = opts_list$mu,
                     weightFcnPtr_x = opts_list$weight_x,
                     weightFcnPtr_t = opts_list$weight_t)
}

parallel.LP.const0 = function(opts_list, Y)
{
  n_core = parallel::detectCores() - 1
  doParallel::registerDoParallel(n_core)
  `%dopar%` = foreach::`%dopar%`
  
  output = foreach::foreach(k = 1:3, .packages = c('DCSmooth')) %dopar% {
    opts_list_i = opts_list[[k]]
    opts_list_i$weight_x = weight_fcn_assign(opts_list_i$weight_x)
    opts_list_i$weight_t = weight_fcn_assign(opts_list_i$weight_t)
    result = LP.const0(opts_list_i, Y)
    return(result) }
  
  doParallel::stopImplicitCluster()
  
  return(output)
}

LP.const1 = function(opts_list, Y)
{
  LP_dcs_const1_BMod(yMat = Y,
                     hVec = opts_list$h,
                     polyOrderVec = opts_list$p,
                     drvVec = opts_list$drv,
                     muVec = opts_list$mu,
                     weightFcnPtr_x = opts_list$weight_x,
                     weightFcnPtr_t = opts_list$weight_t)
}

parallel.LP.const1 = function(opts_list, Y, n_core)
{
  n_core = parallel::detectCores() - 1
  doParallel::registerDoParallel(n_core)
  `%dopar%` = foreach::`%dopar%`
  
  output = foreach::foreach(k = 1:3, .packages = c('DCSmooth')) %dopar% {
    opts_list_i = opts_list[[k]]
    opts_list_i$weight_x = weight_fcn_assign(opts_list_i$weight_x)
    opts_list_i$weight_t = weight_fcn_assign(opts_list_i$weight_t)
    result = LP.const1(opts_list_i, Y)
    return(result) }
  
  doParallel::stopImplicitCluster()
  
  return(output)
}