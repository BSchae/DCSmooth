################################################################################
#                                                                              #
#       DCSmooth Package: Auxiliary Funcions for Bandwidth Selection           #
#                                                                              #
################################################################################

#----------------------Formula for optimal bandwidths--------------------------#

# Local Polynomial Regression
hOptLP = function(mxx, mtt, varCoef, n, nSub, p, kernelProp)
{
  i0x = intCalcLP(mxx, mtt, p, nSub)[1]
  i0t = intCalcLP(mtt, mxx, p, nSub)[1]
  
  # change nSub to n back later?
  hxOpt = (kernelProp$R^2 * varCoef)/((p + 1) * nSub * kernelProp$mu^2 * i0x)
  htOpt = (kernelProp$R^2 * varCoef)/((p + 1) * nSub * kernelProp$mu^2 * i0t)
  
  hxOpt = hxOpt^(1/(2*p + 4))
  htOpt = htOpt^(1/(2*p + 4))
  
  return(c(hxOpt, htOpt))
}

# Kernel Regression
hOptKR = function(mxx, mtt, varCoef, n, nSub, kernelProp)
{
  i0x = intCalcKR(mxx, mtt, nSub)[1]
  i0t = intCalcKR(mtt, mxx, nSub)[1]
  
  hxOpt = (kernelProp$R^2 * varCoef)/(n * kernelProp$mu^2 * i0x)
  htOpt = (kernelProp$R^2 * varCoef)/(n * kernelProp$mu^2 * i0t)
  
  hxOpt = hxOpt^(1/6)
  htOpt = htOpt^(1/6)

  return(c(hxOpt, htOpt))
}

#----------------------Integrals over mxx^2, mtt^2-----------------------------#

intCalcKR = function(m11, m22, nSub)
{
  i11 = sum(m11 * m11)/nSub
  i22 = sum(m22 * m22)/nSub
  i12 = sum(m11 * m22)/nSub
  #print(c(i11, i22, i12))
  
  iOut = (i11/i22)^0.75 * (sqrt(i11 * i22) + i12)
  
  return(c(iOut, i11/i22))
}

intCalcLP = function(m11, m22, p, nSub)
{
  i11 = sum(m11 * m11)/nSub
  i22 = sum(m22 * m22)/nSub
  i12 = sum(m11 * m22)/nSub
  #print(c(i11, i22, i12))
  
  iOut = (i11/i22)^(1/(2*p + 2)) * (2*i11 + sqrt(i11 / i22) + i12)
  
  return(c(iOut, i11/i22))
}

#-------------------------Kernel property calculation--------------------------#

#' @export
kernelPropFcn = function(kernelFcn, nInt = 5000)
{
  uSeq  = seq(from = -1, to = 1, length.out = (2 * nInt + 1))
  valR  = sum((kernelFcn_use(uSeq, q = 1, kernelFcn))^2) / nInt
  valMu = sum((kernelFcn_use(uSeq, q = 1, kernelFcn)) * uSeq^2) / nInt
  
  return(list(R = valR, mu = valMu))
}

#-----------------bandwidth inflation for derivative estimations---------------#

inflationFcnKR = function(h, n, dcsOptions)
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

inflationFcnLP = function(h, n, dcsOptions)
{
  Tx = 1#n[2]/sum(n)
  Tt = 1#n[1]/sum(n)
  hInflxx = c(dcsOptions$inflPar[1] * Tx * (h[1]/Tx)^dcsOptions$inflExp[1],
              dcsOptions$inflPar[2] * Tt * (h[2]/Tt)^dcsOptions$inflExp[2])
  hInfltt = c(dcsOptions$inflPar[2] * Tx * (h[1]/Tx)^dcsOptions$inflExp[2],
              dcsOptions$inflPar[1] * Tt * (h[2]/Tt)^dcsOptions$inflExp[1])

  return(list(h_xx = hInflxx, h_tt = hInfltt))
}