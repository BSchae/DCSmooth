###############################################################################
#                                                                             #
#         DCSmooth Package: R-Functions for LP Bandwidth Selection            #
#                                                                             #
###############################################################################

LP_bndwSelect = function(Y, kernelFcn, dcsOptions)
{
  nX = dim(Y)[1]; nT = dim(Y)[2]
  n  = nX * nT                                  # total number of observations is needed later
  
  kernelProp = kernelPropFcn(kernelFcn)         # calculate properties R and mu_2 of kernel
  
  hOpt = c(0.1, 0.1)                            # initial values for h_0, arbitrary chosen

  iterate = TRUE                                # iteration indicator
  iterationCount = 0
  while(iterate)                                # loop for IPI
  {
    iterationCount = iterationCount + 1
    hOptTemp   = hOpt[1:2]       # store old bandwidths for breaking condition
    hInfl      = inflationFcnLP(hOptTemp, c(nX, nT), dcsOptions)  # inflation of bandwidths for drv estimation
   
    if (dcsOptions$constWindow == TRUE)
    {
      
    } else {
      # # pre-smoothing of the surface function m(0,0) for better estimation of derivatives
      # YSmth = LP_Smooth_test(yMat = Y, hVec = hOptTemp, polyOrderVec
      #         = c(dcsOptions$pOrder, dcsOptions$pOrder),
      #         drvVec = c(0, 0), kernelFcn)
      # 
      # # smoothing of derivatives m(2,0) and m(0,2)
      # mxx = LP_Smooth_test(yMat = Y, hVec = hInfl$h_xx, polyOrderVec
      #       = c(dcsOptions$pOrder + 2, dcsOptions$pOrder), drvVec = c(2, 0),
      #       kernelFcn)
      # mtt = LP_Smooth_test(yMat = Y, hVec = hInfl$h_tt, polyOrderVec
      #       = c(dcsOptions$pOrder, dcsOptions$pOrder + 2), drvVec = c(0, 2),
      #       kernelFcn)
      
      # pre-smoothing of the surface function m(0,0) for better estimation of derivatives
      YSmth = LP_DoubleSmooth2(yMat = Y, hVec = hOptTemp, polyOrderVec
                               = c(dcsOptions$pOrder, dcsOptions$pOrder),
                               drvVec = c(0, 0), kernelFcn)

      # smoothing of derivatives m(2,0) and m(0,2)
      mxx = LP_DoubleSmooth2(yMat = Y, hVec = hInfl$h_xx, polyOrderVec
                             = c(dcsOptions$pOrder + 2, dcsOptions$pOrder), drvVec = c(2, 0),
                             kernelFcn)
      mtt = LP_DoubleSmooth2(yMat = Y, hVec = hInfl$h_tt, polyOrderVec
                             = c(dcsOptions$pOrder, dcsOptions$pOrder + 2), drvVec = c(0, 2),
                             kernelFcn)
    }
      
    # shrink mxx, mtt from boundaries
    if (dcsOptions$delta[1] != 0 || dcsOptions$delta[2] != 0)
    {
      shrinkX = ceiling(dcsOptions$delta[1]*nX):
        floor((1 - dcsOptions$delta[1])*nX)
      shrinkT = ceiling(dcsOptions$delta[2]*nT):
        floor((1 - dcsOptions$delta[2])*nT)
      
      mxx = mxx[shrinkX, shrinkT]
      mtt = mtt[shrinkX, shrinkT]
      nSub = dim(mxx)[1]*dim(mxx)[2]
      
      # calculate variance factor
      if (dcsOptions$varEst == "iid")
      {
        varCoef = (sd((Y - YSmth)[shrinkX, shrinkT]))^2
        qarma_model = NA
      } else if (dcsOptions$varEst == "qarma") {
        cf_est = qarma.cf((Y - YSmth)[shrinkX, shrinkT],
                          model_order = dcsOptions$modelOrder)
        varCoef = cf_est$cf
        qarma_model = cf_est$qarma_model
      }
    } else {
      nSub = n
      
      # calculate variance factor
      if (dcsOptions$varEst == "iid")
      {
        varCoef = (sd(Y - YSmth))^2
        qarma_model = NA
      } else if (dcsOptions$varEst == "qarma") {
        cf_est = qarma.cf((Y - YSmth), model_order = dcsOptions$modelOrder)
        varCoef = cf_est$cf
        qarma_model = cf_est$qarma_model
      }
    }
    
    # # calculate variance factor
    # if (dcsOptions$varEst == "iid")
    # {
    #   varCoef = (sd(Y - YSmth))^2
    #   qarma_model = NA
    # } else if (dcsOptions$varEst == "qarma") {
    #   cf_est = qarma.cf((Y - YSmth), model_order = dcsOptions$modelOrder)
    #   varCoef = cf_est$cf
    #   qarma_model = cf_est$qarma_model
    # }
    
    # calculate optimal bandwidths for next step
    hOpt = hOptLP(mxx, mtt, varCoef, n, nSub, dcsOptions$pOrder, kernelProp)
    
    # break condition
    if( ((hOpt[1]/hOptTemp[1] - 1 < 0.001) && (hOpt[2]/hOptTemp[2] - 1 
                    < 0.001) && (iterationCount > 2)) || iterationCount > 10)
    {
      iterate = FALSE
    }
  }
  return(list(bndw = hOpt, iterations = iterationCount, varCoef = varCoef,
              qarma_model = qarma_model))
}