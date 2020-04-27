###############################################################################
#                                                                             #
#         DCSmooth Package: R-Functions for KR Bandwidth Selection            #
#                                                                             #
###############################################################################


#-----------------------Set Options via Function------------------------------#

.setOptions = function(    # inside function with default values in arguments
  kernPar   = c(2, 2),     # choose a kernel function with mu = 2, 
  pOrder    = 1,           # choose order of polynomials for X-/T-smoothing
  # (orders have to be the same in both directions)
  inflExp   = c(0.5, 0.5), # inflation exponent
  inflPar   = c(1.5, 0.5)  # inflation parameters c (regression), d (2nd derivative)
)
{
  return(list(kernPar = kernPar, pOrder = pOrder,
              inflExp = inflExp, inflPar = inflPar))
}

#------------------Function for the optimal bandwidth via IPI-----------------#

KR_bndwSelect = function(Y, kernelFcn, dcsOptions)
{
  nX = dim(Y)[1]; nT = dim(Y)[2]
  n  = nX * nT                                  # total number of observations is needed later
  
  kernelProp = kernelPropFcn(kernelFcn)         # calculate properties R and mu_2 of kernel
  
  kernFcn0 = kernelFcn_assign("MW420")
  kernFcn2 = kernelFcn_assign("MW422")
  
  hOpt = c(0.05, 0.05)                          # initial values for h_0, arbitrary chosen
  
  iterate = TRUE                                # iteration indicator
  iterationCount = 0
  while(iterate)                                # loop for IPI
  {
    iterationCount = iterationCount + 1
    hOptTemp   = pmin(hOpt[1:2], c(0.5, 0.5))        # store old bandwidths for breaking condition
    hInfl      = inflationFcn(hOptTemp, c(nX, nT), dcsOptions)  # inflation of bandwidths for drv estimation
    hInfl$h_xx = pmin(hInfl$h_xx, c(0.45, 0.45))
    hInfl$h_tt = pmin(hInfl$h_tt, c(0.45, 0.45))# ensure no bandwidth > 0.5 (or 0.45)
    
    # pre-smoothing of the surface function m(0,0) for better estimation of derivatives
    YSmth = KR_DoubleSmooth2(yMat = Y, hVec = hOptTemp,
                             kernFcnPtrX = kernFcn0, kernFcnPtrT = kernFcn0)
    
    # smoothing of derivatives m(2,0) and m(0,2)
    mxx = KR_DoubleSmooth2(yMat = YSmth, hVec = hInfl$h_xx,
                           kernFcnPtrX = kernFcn2, kernFcnPtrT = kernFcn0)
    mtt = KR_DoubleSmooth2(yMat = YSmth, hVec = hInfl$h_tt,
                           kernFcnPtrX = kernFcn0, kernFcnPtrT = kernFcn2)
    
    # calculate variance factor
    varCoef = (sd(Y - YSmth))^2
    
    # calculate optimal bandwidths for next step
    hOpt = hOptFunction(mxx, mtt, varCoef, n, dcsOptions$pOrder, kernelProp)
    
    # break condition
    if( ((hOpt[1]/hOptTemp[1] - 1 < 0.001) && (hOpt[2]/hOptTemp[2] - 1 
                                               < 0.001)) || iterationCount > 15)
    {
      iterate = FALSE
    }
  }
  return(c(hOpt, iterationCount))
}

#----------------------Formula for optimal bandwidths-------------------------#

hOptFunction = function(mxx, mtt, varCoef, n, p, kernelProp)
{
  integrals = intCalc(mxx, mtt, n, p)[2:4]
  i0x = intCalc(mxx, mtt, n, p)[1]
  i0t = intCalc(mtt, mxx, n, p)[1]
  
  hxOpt = (kernelProp$R^2 * varCoef)/((p + 1) * n * kernelProp$mu^2 * i0x)
  htOpt = (kernelProp$R^2 * varCoef)/((p + 1) * n * kernelProp$mu^2 * i0t)
  
  hxOpt = hxOpt^(1/(2*p + 4))
  htOpt = htOpt^(1/(2*p + 4))
  
  return(c(hxOpt, htOpt, integrals))
}

#---------------------Integrals over mxx^2, mtt^2-----------------------------#

intCalc2 = function(m11, m22, n, p)
{
  i11 = sum(m11 * m11)/n
  i22 = sum(m22 * m22)/n
  i12 = sum(m11 * m22)/n
  
  iOut = (i11/i22)^0.75 * (sqrt(i11 * i22) + i12)
  
  return(iOut)
}

intCalc = function(m11, m22, n, p)
{
  i11 = sum(m11 * m11)/n
  i22 = sum(m22 * m22)/n
  i12 = sum(m11 * m22)/n
  
  iOut = (i11/i22)^0.75 * (sqrt(i11 * i22) + i12)
  
  return(c(iOut, i11, i22, i12))
}


#------------------------Kernel property calculation--------------------------#

kernelPropFcn = function(kernelFcn, nInt = 5000)
{
  uSeq  = seq(from = -1, to = 1, length.out = (2 * nInt + 1))
  valR  = sum((kernelFcn_use(uSeq, q = 1, kernelFcn))^2) / nInt
  valMu = sum((kernelFcn_use(uSeq, q = 1, kernelFcn)) * uSeq^2) / nInt
  
  return(list(R = valR, mu = valMu))
}

#----------------bandwidth inflation for derivative estimations---------------#

inflationFcn = function(h, n, dcsOptions)
{
  Tx = 1#n[2]/sum(n)
  Tt = 1#n[1]/sum(n)
  hInflxx = c(dcsOptions$inflPar[1] * Tx * (h[1]/Tx)^dcsOptions$inflExp[1],
              dcsOptions$inflPar[2] * Tt * (h[2]/Tt)^dcsOptions$inflExp[2])
  hInfltt = c(dcsOptions$inflPar[2] * Tx * (h[1]/Tx)^dcsOptions$inflExp[2],
              dcsOptions$inflPar[1] * Tt * (h[2]/Tt)^dcsOptions$inflExp[1])
  
  return(list(h_xx = hInflxx, h_tt = hInfltt))
}
