###############################################################################
#                                                                             #
#         DCSmooth Package: R-Functions for KR Bandwidth Selection            #
#                                                                             #
###############################################################################

#------------------Function for the optimal bandwidth via IPI-----------------#

KR.bndw = function(Y, kernel_fcn, kernel_fcn_t, dcs_options)
{
  n_x = dim(Y)[1]; n_t = dim(Y)[2]
  n  = n_x * n_t                                  # total number of observations is needed later
  
  kernel_prop = kernel.prop.KR(kernel_fcn)         # calculate properties R and mu_2 of kernel
  
  kern_fcn_0 = kernel_fcn_assign("MW220")          # assign kernel for regression surface
  kern_fcn_2 = kernel_fcn_assign("MW422")          # assign kernel for 2nd derivative (needed in )
  
  h_opt = c(0.1, 0.1)                            # initial values for h_0, arbitrary chosen
  
  iterate = TRUE                                # iteration indicator
  iteration_count = 0
  while(iterate)                                # loop for IPI
  {
    iteration_count = iteration_count + 1
    h_opt_temp   = pmin(h_opt[1:2], c(0.45, 0.45))        # store old bandwidths for breaking condition
    h_infl  = inflation.KR(h_opt_temp, c(nX, nT), dcs_options)  # inflation of bandwidths for drv estimation
    
    if (dcs_options$const_window == TRUE)
    {
      # pre-smoothing of the surface function m(0,0) for estimation of variance
      Y_smth = KR_dcs_const1(yMat = Y, hVec = h_opt_temp, drvVec = c(0, 0),
                        kernFcnPtrX = kern_fcn_0, kernFcnPtrT = kern_fcn_0)

      # smoothing of derivatives m(2,0) and m(0,2)
      mxx = KR_dcs_const1(yMat = Y, hVec = h_infl$h_xx, drvVec = c(2, 0),
                        kernFcnPtrX = kern_fcn_2, kernFcnPtrT = kern_fcn_0)

      mtt = KR_dcs_const1(yMat = Y, hVec = h_infl$h_tt, drvVec = c(0, 2),
                    kernFcnPtrX = kern_fcn_0, kernFcnPtrT = kern_fcn_2)
    } else if (dcs_options$const_window == FALSE) {
      # pre-smoothing of the surface function m(0,0) for estimation of variance
      Y_smth = KR_dcs_const0(yMat = Y, hVec = h_opt_temp, drvVec = c(0, 0),
                             kernFcnPtrX = kern_fcn_0, kernFcnPtrT = kern_fcn_0)
      
      # smoothing of derivatives m(2,0) and m(0,2)
      mxx = KR_dcs_const0(yMat = Y, hVec = h_infl$h_xx, drvVec = c(2, 0),
                          kernFcnPtrX = kern_fcn_2, kernFcnPtrT = kern_fcn_0)
      
      mtt = KR_dcs_const0(yMat = Y, hVec = h_infl$h_tt, drvVec = c(0, 2),
                          kernFcnPtrX = kern_fcn_0, kernFcnPtrT = kern_fcn_2)
    }

    # shrink mxx, mtt from boundaries
    if (dcs_options$delta[1] != 0 || dcs_options$delta[2] != 0)
    {
      shrink_x = ceiling(dcs_options$delta[1] * n_x):
                      (n_x - floor(dcs_options$delta[1] * n_x))
      shrink_t = ceiling(dcs_options$delta[2] * n_t):
                      (n_t - floor(dcs_options$delta[2] * n_t))

      mxx = mxx[shrink_x, shrink_t]
      mtt = mtt[shrink_x, shrink_t]
      n_sub = dim(mxx)[1]*dim(mxx)[2]
    } else {
      n_sub = n
    }
      
    # calculate variance factor
    var_est = cf.estimation(Y - Y_smth, dcs_options)
    var_coef = var_est$cf_est
    stat_test = var_est$stationary
    
    # calculate optimal bandwidths for next step
    h_opt = h.opt.KR(mxx, mtt, var_coef, n, n_sub, kernel_prop)
    
    # break condition
    if( ((h_opt[1]/h_opt_temp[1] - 1 < 0.001) && (h_opt[2]/h_opt_temp[2] - 1 
        < 0.001) && (iteration_count > 3)) || (iteration_count > 15) )
    {
      iterate = FALSE
    }
  }
  return(list(h_opt = h_opt, iterations = iteration_count, var_coef = var_coef,
              qarma_model = var_est$qarma_model, stat_test = stat_test))
}
