################################################################################
#                                                                              #
#                DCSmooth Package: estimation of cf coefficients               #
#                                                                              #
################################################################################

# Different methods for estimation of the cf coefficient for bandwidth
# selection are joined in this function

# Y should be the residuals, e.g. Y - YSmth in the estimation codes
cf.estimation = function(Y, dcs_options)
{
  if (dcs_options$var_est == "iid")
  {
    cf_est = sd(Y)^2
    cf_out = list(cf_est = cf_est, var_est = "iid", stationary = TRUE)
  } else if (dcs_options$var_est == "qarma") {
    if (!is.list(dcs_options$qarma_order) && dcs_options$qarma_order == "gpac")
    {
      model_order = qarma.order_gpac(Y, order_max = dcs_options$order_max)
      qarma = qarma.cf(Y, model_order = model_order)
      cf_est = qarma$cf
      qarma_model = list(ar = qarma$qarma_model$ar, ma = qarma$qarma_model$ma)
      cf_out = list(cf_est = cf_est, qarma_model = qarma_model,
                    var_est = "qarma_gpac", stationary =
                    qarma$qarma_model$stationary)
    } else if (!is.list(dcs_options$qarma_order) &&
               dcs_options$qarma_order == "bic") {
      model_order = qarma.order_bic(Y, order_max = dcs_options$order_max)
      qarma = qarma.cf(Y, model_order = model_order)
      cf_est = qarma$cf
      qarma_model = list(ar = qarma$qarma_model$ar, ma = qarma$qarma_model$ma)
      cf_out = list(cf_est = cf_est, qarma_model = qarma_model,
                    var_est = "qarma_bic", stationary = 
                    qarma$qarma_model$stationary)
    } else {
      qarma = qarma.cf(Y, model_order = dcs_options$qarma_order)
      cf_est = qarma$cf
      qarma_model = list(ar = qarma$qarma_model$ar, ma = qarma$qarma_model$ma)
      cf_out = list(cf_est = cf_est, qarma_model = qarma_model,
                    var_est = "qarma", stationary =
                    qarma$qarma_model$stationary)
    }
  } else if (dcs_options$var_est == "nonpar") {
    cf_est = specDens((Y - YSmth), omega = c(0, 0))
    cf_out = list(cf_est = cf_est$cf, var_est = "nonpar")
  }
  
  return(cf_out)
}