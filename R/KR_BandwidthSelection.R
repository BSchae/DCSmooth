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
  
  hOpt = c(1/nX, 1/nT)                          # initial values for h_0, arbitrary chosen
  
  iterate = TRUE                                # iteration indicator
  iterationCount = 0
  while(iterate)                                # loop for IPI
  {
    iterationCount = iterationCount + 1
    hOptTemp   = pmin(hOpt[1:2], c(0.45, 0.45))        # store old bandwidths for breaking condition
    hInfl      = inflationFcn(hOptTemp, c(nX, nT), dcsOptions)  # inflation of bandwidths for drv estimation

    # pre-smoothing of the surface function m(0,0) for better estimation of derivatives
    YSmth = KR_DoubleSmooth2(yMat = Y, hVec = hOptTemp, drvVec = c(0, 0),
                             kernFcnPtrX = kernFcn0, kernFcnPtrT = kernFcn0)
    
    # smoothing of derivatives m(2,0) and m(0,2)
    mxx = KR_DoubleSmooth2(yMat = YSmth, hVec = hInfl$h_xx, drvVec = c(2, 0),
                           kernFcnPtrX = kernFcn2, kernFcnPtrT = kernFcn0)
    mtt = KR_DoubleSmooth2(yMat = YSmth, hVec = hInfl$h_tt, drvVec = c(0, 2),
                           kernFcnPtrX = kernFcn0, kernFcnPtrT = kernFcn2)
    
    # calculate variance factor
    varCoef = (sd(Y - YSmth))^2
    
    # calculate optimal bandwidths for next step
    hOpt = hOptKR(mxx, mtt, varCoef, n, kernelProp)
    
    # break condition
    if( ((hOpt[1]/hOptTemp[1] - 1 < 0.001) && (hOpt[2]/hOptTemp[2] - 1 
                                               < 0.001)) || iterationCount > 10)
    {
      iterate = FALSE
    }
  }
  return(c(hOpt, iterationCount))
}
