###############################################################################
#                                                                             #
#       DCSmooth Package: Auxiliary Funcions for Bandwidth Selection          #
#                                                                             #
###############################################################################

#-----------------------Set Options via Function------------------------------#

.setOptions = function(    # inside function with default values in arguments
  kernPar   = c(2, 2),     # choose a kernel function with mu = 2, 
  pOrder    = 1,           # choose order of polynomials for X-/T-smoothing
                           # if pOrder == 0, kernel regression will be used.
  # (orders have to be the same in both directions)
  inflExp   = c(0.5, 0.5), # inflation exponent
  inflPar   = c(1.5, 0.5)  # inflation parameters c (regression), d (2nd derivative)
  
)
{
  return(list(kernPar = kernPar, pOrder = pOrder,
              inflExp = inflExp, inflPar = inflPar))
}

#----------------------Formula for optimal bandwidths-------------------------#

hOptFunction = function(mxx, mtt, varCoef, n, p, kernelProp)
{
  i0x = intCalc(mxx, mtt, n)
  i0t = intCalc(mtt, mxx, n)
  
  hxOpt = (kernelProp$R^2 * varCoef)/((p + 1) * n * kernelProp$mu^2 * i0x)
  htOpt = (kernelProp$R^2 * varCoef)/((p + 1) * n * kernelProp$mu^2 * i0t)
  
  hxOpt = hxOpt^(1/(2*p + 4))
  htOpt = htOpt^(1/(2*p + 4))
  
  return(c(hxOpt, htOpt))
}

hOptKR = function(mxx, mtt, varCoef, n, kernelProp)
{
  i0x = intCalc(mxx, mtt, n)
  i0t = intCalc(mtt, mxx, n)
  
  hxOpt = (kernelProp$R^2 * varCoef)/(n * kernelProp$mu^2 * i0x)
  htOpt = (kernelProp$R^2 * varCoef)/(n * kernelProp$mu^2 * i0t)
  
  hxOpt = hxOpt^(1/6)
  htOpt = htOpt^(1/6)
  
  return(c(hxOpt, htOpt))
}

#---------------------Integrals over mxx^2, mtt^2-----------------------------#

intCalc = function(m11, m22, n)
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

inflationFcn = function(h, n, dcsOptions)
{
  Tx = 1#n[2]/sum(n)
  Tt = 1#n[1]/sum(n)
  hInflxx = c(dcsOptions$inflPar[1] * Tx * (h[1]/Tx)^dcsOptions$inflExp[1],
              dcsOptions$inflPar[2] * Tt * (h[2]/Tt)^dcsOptions$inflExp[2])
  hInfltt = c(dcsOptions$inflPar[2] * Tx * (h[1]/Tx)^dcsOptions$inflExp[2],
              dcsOptions$inflPar[1] * Tt * (h[2]/Tt)^dcsOptions$inflExp[1])
  
  hInflxx = pmin(hInflxx, c(0.45, 0.45))
  hInfltt = pmin(hInfltt, c(0.45, 0.45))  # ensure no bandwidth > 0.5 (or 0.45)
  
  return(list(h_xx = hInflxx, h_tt = hInfltt))
}