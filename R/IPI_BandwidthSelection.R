###############################################################################
#                                                                             #
#         DCSmooth Package: R-Functions for IPI Bandwidth Selection           #
#                                                                             #
###############################################################################


#-----------------------Set Options via Function------------------------------#

.setOptions = function(    # inside function with default values in arguments
  kernPar   = c(2, 2),     # choose a kernel function with mu = 2, 
  pOrder    = 1,           # choose order of polynomials for X-/T-smoothing
                           # (orders have to be the same in both directions)
  inflExp   = 0.5,         # inflation exponent
  inflPar   = c(1.5, 0.5)  # inflation parameters c (regression), d (2nd derivative)
)
{
  return(list(kernPar = kernPar, pOrder = pOrder,
              inflExp = inflExp, inflPar = inflPar))
}

#------------------Function for the optimal bandwidth via IPI-----------------#

bndwSelect = function(Y, kernelFcn, dcsOptions)
{
  n = dim(Y)[1] * dim(Y)[2]    # total number of observations is needed later
  
  kernelProp = kernelPropFcn(kernelFcn)         # calculate properties R and mu_2 of kernel
  
  hOpt = c(0.01, 0.01)                          # initial values for h_0, arbitrary chosen

  iterate = TRUE                                # iteration indicator
  iterationCount = 0
  while(iterate)                                # loop for IPI
  {
    iterationCount = iterationCount + 1
    hOptTemp = hOpt                             # store old bandwidths for breaking condition
    hInfl = inflationFcn(hOptTemp, dcsOptions)  # inflation of bandwidths for drv estimation
    hInfl$h_xx = pmin(hInfl$h_xx, c(0.45, 0.45))
    hInfl$h_tt = pmin(hInfl$h_tt, c(0.45, 0.45))
                                                # ensure no bandwidth > 0.5 (or 0.45)
      
    # pre-smoothing of the surface function m(0,0) for better estimation of derivatives
    YSmth = FastDoubleSmooth(yMat = Y, hVec = hOptTemp, polyOrderVec 
                             = c(dcsOptions$pOrder, dcsOptions$pOrder), drvVec = c(0, 0))
    
    # smoothing of derivatives m(2,0) and m(0,2)
    mxx = FastDoubleSmooth(yMat = YSmth, hVec = hInfl$h_xx, polyOrderVec 
                      = c(dcsOptions$pOrder + 2, dcsOptions$pOrder), drvVec = c(2, 0))
    mtt = FastDoubleSmooth(yMat = YSmth, hVec = hInfl$h_tt, polyOrderVec 
      = c(dcsOptions$pOrder, dcsOptions$pOrder + 2), drvVec = c(0, 2))
    
    # calculate variance factor
    varCoef = (sd(Y - YSmth))^2
    
    # calculate optimal bandwidths for next step
    hOpt = hOptFunction(mxx, mtt, varCoef, n, dcsOptions$pOrder, kernelProp)
    
    # break condition
    if( (hOpt[1]/hOptTemp[1] - 1 < 0.001) && (hOpt[2]/hOptTemp[2] - 1 < 0.001) )
    {
      iterate = FALSE
    }
  }
  
  return(hOpt)
}

#----------------------Formula for optimal bandwidths-------------------------#
  
hOptFunction = function(mxx, mtt, varCoef, n, p, kernelProp)
{
  i0x = intCalc(mxx, mtt, n, p)
  i0t = intCalc(mtt, mxx, n, p)
  
  hxOpt = (kernelProp$R^2 * varCoef)/((p + 1) * n * kernelProp$mu^2 * i0x)
  htOpt = (kernelProp$R^2 * varCoef)/((p + 1) * n * kernelProp$mu^2 * i0x)

  hxOpt = hxOpt^(1/(2*p + 4))
  htOpt = htOpt^(1/(2*p + 4))
  
  return(c(hxOpt, htOpt))
}

#---------------------Integrals over mxx^2, mtt^2-----------------------------#

intCalc = function(m11, m22, n, p)
{
  i11 = sum(m11 * m11)/n
  i22 = sum(m22 * m22)/n
  i12 = sum(m11 * m22)/n
  
  iOut = (i11/i22)^0.75 * (sqrt(i11 * i22) + i12)

  return(iOut)
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

inflationFcn = function(h, dcsOptions)
{
  hInflxx = c(dcsOptions$inflPar[2]*h[1]^dcsOptions$inflExp,
                dcsOptions$inflPar[1]*h[2]^dcsOptions$inflExp)
  hInfltt = c(dcsOptions$inflPar[1]*h[1]^dcsOptions$inflExp,
           dcsOptions$inflPar[2]*h[2]^dcsOptions$inflExp)
  
  return(list(h_xx = hInflxx, h_tt = hInfltt))
}

#----------------function for fast calc of true bndw--------------------------#


trueBndw = function(I, kernProp, n, p, varCoef)
{
  I0 = (I[1]/I[2])^0.75 * (sqrt(I[1]*I[2]) + I[3])
  hx = ((kernProp$R^2 * varCoef) / (n * (p + 1) * kernProp$mu^2 * I0))^(1/(4 + 2*p))
  ht = (I[1]/I[2])^0.25 * hx
  
  return(c(hx, ht))
}
