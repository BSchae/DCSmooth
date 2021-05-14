################################################################################
#                                                                              #
#          DCSmooth Package: R-Functions for LP Bandwidth Selection            #
#                                                                              #
################################################################################

LP_bndwSelect = function(Y, kernelFcn, dcsOptions)
{
  nX = dim(Y)[1]; nT = dim(Y)[2]
  n  = nX * nT              # total number of observations is needed later
  pOrder = dcsOptions$pOrder
  drvVec = dcsOptions$drv
  
  hOpt = c(0.1, 0.1)        # initial values for h_0, arbitrary chosen

  iterate = TRUE            # iteration indicator
  iterationCount = 0
  while(iterate)            # loop for IPI
  {
    iterationCount = iterationCount + 1
    hOptTemp   = hOpt[1:2]      # store old bandwidths for breaking condition
    hInfl      = inflationFcnLP(hOptTemp, c(nX, nT), dcsOptions)  
                                # inflation of bandwidths for drv estimation
   
    if (dcsOptions$constWindow == TRUE)
    {
      # smoothing of Y for variance factor estimation
      YSmth = LP_DoubleSmooth(yMat = Y, hVec = hOptTemp, polyOrderVec
               = pOrder, drvVec = drvVec, kernelFcn)
      
      # smoothing of derivatives m(2,0) and m(0,2)
      mxx = LP_DoubleSmooth(yMat = Y, hVec = hInfl$h_xx, polyOrderVec = 
                c(2*pOrder[1] - drvVec[1] + 1, pOrder[2]), drvVec =
                c(pOrder[1] + 1, drvVec[2]), kernelFcn)
      mtt = LP_DoubleSmooth(yMat = Y, hVec = hInfl$h_tt, polyOrderVec =
                c(pOrder[1], 2*pOrder[2] - drvVec[2] + 1), drvVec =
                c(drvVec[1], pOrder[2] + 1), kernelFcn)
    } else {
      # smoothing of Y for variance factor estimation
      YSmth = LP_DoubleSmooth2(yMat = Y, hVec = hOptTemp, polyOrderVec
                              = pOrder, drvVec = drvVec, kernelFcn)
      
      # smoothing of derivatives m(2,0) and m(0,2)
      # needs update for p even.
      mxx = LP_DoubleSmooth2(yMat = Y, hVec = hInfl$h_xx, polyOrderVec = 
                c(2*pOrder[1] - drvVec[1] + 1, pOrder[2]), drvVec =
                c(pOrder[1] + 1, drvVec[2]), kernelFcn)
      mtt = LP_DoubleSmooth2(yMat = Y, hVec = hInfl$h_tt, polyOrderVec =
                c(pOrder[1], 2*pOrder[2] - drvVec[2] + 1), drvVec =
                c(drvVec[1], pOrder[2] + 1), kernelFcn)
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
      
    } else {
      nSub = n
    }  
     
    # calculate variance factor
    varEst = cf.estimation(Y - YSmth, dcsOptions)
    varCoef = varEst$cf_est
    
    # calculate optimal bandwidths for next step
    hOpt = hOptLP(mxx, mtt, varCoef, nSub, pOrder, drvVec, kernelFcn)
    
    # break condition
    if( ((hOpt[1]/hOptTemp[1] - 1 < 0.001) && (hOpt[2]/hOptTemp[2] - 1 
                    < 0.001) && (iterationCount > 2)) || iterationCount > 10)
    {
      iterate = FALSE
    }
  }
  return(list(bndw = hOpt, iterations = iterationCount, varCoef = varCoef,
              qarma_model = varEst$qarma_model))
}