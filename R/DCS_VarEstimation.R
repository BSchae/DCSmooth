################################################################################
#                                                                              #
#                DCSmooth Package: estimation of cf coefficients               #
#                                                                              #
################################################################################

# Different methods for estimation of the cf coefficient for bandwidth
# selection are joined in this function

# Y should be the residuals, e.g. Y - YSmth in the estimation codes
cf.estimation = function(Y, dcsOptions)
{
  if (dcsOptions$varEst == "iid")
  {
    cf_est = sd(Y)^2
    cf_out = list(cf_est = cf_est, varEst = "iid")
  } else if (dcsOptions$varEst == "qarma") {
    if (!is.list(dcsOptions$modelOrder) && dcsOptions$modelOrder == "gpac")
    {
      model_order = qarma.order_gpac(Y, order_max = dcsOptions$orderMax)
      qarma = qarma.cf(Y, model_order = model_order)
      cf_est = qarma$cf
      qarma_model = list(ar = qarma$qarma_model$ar, ma = qarma$qarma_model$ma)
      cf_out = list(cf_est = cf_est, qarma_model = qarma_model,
                    varEst = "qarma_gpac")
    } else if (!is.list(dcsOptions$modelOrder) &&
               dcsOptions$modelOrder == "bic") {
      model_order = qarma.order_bic(Y, order_max = dcsOptions$orderMax)
      qarma = qarma.cf(Y, model_order = model_order)
      cf_est = qarma$cf
      qarma_model = list(ar = qarma$qarma_model$ar, ma = qarma$qarma_model$ma)
      cf_out = list(cf_est = cf_est, qarma_model = qarma_model,
                    varEst = "qarma_bic")
    } else {
      qarma = qarma.cf(Y, model_order = dcsOptions$modelOrder)
      cf_est = qarma$cf
      qarma_model = list(ar = qarma$qarma_model$ar, ma = qarma$qarma_model$ma)
      cf_out = list(cf_est = cf_est, qarma_model = qarma_model,
                    varEst = "qarma")
    }
  } else if (dcsOptions$varEst == "nonpar") {
    cf_est = specDens((Y - YSmth), omega = c(0, 0))
    cf_out = list(cf_est = cf_est$cf, varEst = "nonpar")
  }
  
  return(cf_out)
}