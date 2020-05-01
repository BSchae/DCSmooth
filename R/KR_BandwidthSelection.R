###############################################################################
#                                                                             #
#         DCSmooth Package: R-Functions for KR Bandwidth Selection            #
#                                                                             #
###############################################################################

#------------------Function for the optimal bandwidth via IPI-----------------#

KR_bndwSelect = function(Y, kernelFcn, dcsOptions)
{
  nX = dim(Y)[1]; nT = dim(Y)[2]
  n  = nX * nT                                  # total number of observations is needed later
  
  kernelProp = kernelPropFcn(kernelFcn)         # calculate properties R and mu_2 of kernel
  
  kernFcn0 = kernelFcn_assign("MW420")
  kernFcn2 = kernelFcn_assign("MW422")
  
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

    # pre-smoothing of the surface function m(0,0) for better estimation of derivatives
    YSmth = KR_DoubleSmooth2(yMat = YSub, hVec = hOptTemp, drvVec = c(0, 0),
                             kernFcnPtrX = kernFcn0, kernFcnPtrT = kernFcn0)
    
    # smoothing of derivatives m(2,0) and m(0,2)
    mxx = KR_DoubleSmooth2(yMat = YSub, hVec = hInfl$h_xx, drvVec = c(2, 0),
                           kernFcnPtrX = kernFcn2, kernFcnPtrT = kernFcn0)
    mtt = KR_DoubleSmooth2(yMat = YSub, hVec = hInfl$h_tt, drvVec = c(0, 2),
                           kernFcnPtrX = kernFcn0, kernFcnPtrT = kernFcn2)
    
    # calculate variance factor
    varCoef = (sd(YSub - YSmth))^2
    
    # calculate optimal bandwidths for next step
    hOpt = hOptKR(mxx, mtt, varCoef, n, nSub, kernelProp)
    
    # break condition
    if( ((hOpt[1]/hOptTemp[1] - 1 < 0.001) && (hOpt[2]/hOptTemp[2] - 1 
                    < 0.001) && (iterationCount > 3)) || (iterationCount > 10) )
    {
      iterate = FALSE
    }
  }
  return(c(hOpt, iterationCount))
}
