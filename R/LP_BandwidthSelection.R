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
  
  if (dcsOptions$fast == TRUE) {
    sX = floor(nX/1000) + 1
    sT = floor(nT/1000) + 1
    
    xValues = 1:(nX/sX)*sX
    tValues = 1:(nT/sT)*sT
    
    YSub = Y[xValues, tValues]
  }
  else
  {
    YSub = Y
  }
  
  nXSub = dim(YSub)[1]; nTSub = dim(YSub)[2]; nSub = nXSub * nTSub
  hOpt = c(1/nXSub, 1/nTSub)                    # initial values for h_0, arbitrary chosen

  iterate = TRUE                                # iteration indicator
  iterationCount = 0
  while(iterate)                                # loop for IPI
  {
    iterationCount = iterationCount + 1
    hOptTemp   = pmin(hOpt[1:2], c(0.45, 0.45))        # store old bandwidths for breaking condition
    hInfl      = inflationFcn(hOptTemp, c(nX, nT), dcsOptions)  # inflation of bandwidths for drv estimation
    hInfl$h_xx = pmin(hInfl$h_xx, c(0.45, 0.45))
    hInfl$h_tt = pmin(hInfl$h_tt, c(0.45, 0.45))# ensure no bandwidth > 0.5 (or 0.45)
    
    # pre-smoothing of the surface function m(0,0) for better estimation of derivatives
    YSmth = LP_DoubleSmooth(yMat = YSub, hVec = hOptTemp, polyOrderVec 
                      = c(dcsOptions$pOrder, dcsOptions$pOrder), drvVec = c(0, 0))
    
    # smoothing of derivatives m(2,0) and m(0,2)
    mxx = LP_DoubleSmooth(yMat = YSmth, hVec = hInfl$h_xx, polyOrderVec 
           = c(dcsOptions$pOrder + 2, dcsOptions$pOrder), drvVec = c(2, 0))
    mtt = LP_DoubleSmooth(yMat = YSmth, hVec = hInfl$h_tt, polyOrderVec 
           = c(dcsOptions$pOrder, dcsOptions$pOrder + 2), drvVec = c(0, 2))
    
    # calculate variance factor
    varCoef = (sd(YSub - YSmth))^2
    
    # calculate optimal bandwidths for next step
    hOpt = hOptFunction(mxx, mtt, varCoef, n, nSub, dcsOptions$pOrder, kernelProp)
    
    # break condition
    if( ((hOpt[1]/hOptTemp[1] - 1 < 0.001) && (hOpt[2]/hOptTemp[2] - 1 
                                               < 0.001)) || iterationCount > 10)
    {
      iterate = FALSE
    }
  }
  return(c(hOpt, iterationCount))
}