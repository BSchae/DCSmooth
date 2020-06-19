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
  
  hOpt = c(0.1, 0.1)                            # initial values for h_0, arbitrary chosen
  
  iterate = TRUE                                # iteration indicator
  iterationCount = 0
  while(iterate)                                # loop for IPI
  {
    iterationCount = iterationCount + 1
    hOptTemp   = pmin(hOpt[1:2], c(0.45, 0.45))        # store old bandwidths for breaking condition
    hInfl  = inflationFcnKR(hOptTemp, c(nX, nT), dcsOptions)  # inflation of bandwidths for drv estimation
    
    # pre-smoothing of the surface function m(0,0) for better estimation of derivatives
    YSmth = KR_DoubleSmooth(yMat = Y, hVec = hOptTemp, drvVec = c(0, 0),
                      kernFcnPtrX = kernFcn0, kernFcnPtrT = kernFcn0)

        # smoothing of derivatives m(2,0) and m(0,2)
    mxx = KR_DoubleSmooth(yMat = Y, hVec = hInfl$h_xx, drvVec = c(2, 0),
                           kernFcnPtrX = kernFcn2, kernFcnPtrT = kernFcn0)

    mtt = KR_DoubleSmooth(yMat = YSmth, hVec = hInfl$h_tt, drvVec = c(0, 2),
                    kernFcnPtrX = kernFcn0, kernFcnPtrT = kernFcn2)

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
      
    # # Short Memory
    # order.arma = list(ARMAx = c(1, 1), ARMAt = c(1, 1))
    # varCoef = DCS.cf(Y - YSmth, order.arma)$cfOut

    # calculate variance factor
    varCoef = (sd(Y - YSmth))^2
    
    # calculate optimal bandwidths for next step
    hOpt = hOptKR(mxx, mtt, varCoef, n, nSub, kernelProp)
    
    # break condition
    if( ((hOpt[1]/hOptTemp[1] - 1 < 0.001) && (hOpt[2]/hOptTemp[2] - 1 
                    < 0.001) && (iterationCount > 3)) || (iterationCount > 15) )
    {
      iterate = FALSE
    }
  }
  return(list(bndw = hOpt, iterations = iterationCount))
}
