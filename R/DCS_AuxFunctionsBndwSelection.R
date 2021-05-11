################################################################################
#                                                                              #
#       DCSmooth Package: Auxiliary Funcions for Bandwidth Selection           #
#                                                                              #
################################################################################

#----------------------Formula for optimal bandwidths--------------------------#

# Local Polynomial Regression
hOptLP = function(mxx, mtt, varCoef, nSub, pVec, drvVec, kernelFcn)
{
  # calculation of integrals
  i11 = sum(mxx^2)/nSub; i22 = sum(mtt^2)/nSub; i12 = sum(mxx * mtt)/nSub
  
  # kernel constants (kernel Functions may also depend on p, drv)
  kernelProp1 = kernelPropLP(kernelFcn, pVec[1], drvVec[1])
  kernelProp2 = kernelPropLP(kernelFcn, pVec[2], drvVec[2])
  
  # relation factor gamma_21 (h_1 = gamma_21 * h_2)
  delta = (pVec - drvVec)[1] # should be the same for both entries
  gamma_21 = (kernelProp1$mu/kernelProp2$mu)^(1/(delta + 1)) *
              ( i12/i11 * (drvVec[1] - drvVec[2])/(2*drvVec[2] + 1) +
              sqrt(i12^2/i11^2 * (drvVec[1] - drvVec[2])^2/(2*drvVec[2] + 1)^2 +
              i22/i11 * (2*drvVec[1] + 1)/(2*drvVec[2] + 1)) )^(1/(delta + 1))
  gamma_12 = 1/gamma_21
  
  # optimal bandwidths
  I1 = kernelProp2$mu^2 * i22 + kernelProp1$mu * kernelProp2$mu * 
       i12 * gamma_21^(delta + 1)
  hxOpt = (2*drvVec[1] + 1)/(2*(delta + 1)) * (kernelProp1$R * kernelProp2$R*
          varCoef)/(nSub * gamma_21^(2*drvVec[1] + 1) * I1) 
  I2 = kernelProp1$mu^2 * i11 + kernelProp1$mu * kernelProp2$mu * 
       i12 * gamma_12^(delta + 1)
  htOpt = (2*drvVec[2] + 1)/(2*(delta + 1)) * (kernelProp1$R * kernelProp2$R*
          varCoef)/(nSub * gamma_12^(2*drvVec[2] + 1) * I2) 
  
  hxOpt = hxOpt^(1/(2*(delta + sum(drvVec) + 2)))
  htOpt = htOpt^(1/(2*(delta + sum(drvVec) + 2)))
  
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

intCalcLP = function(m11, m22, p, nSub)
{
  i11 = sum(m11 * m11)/nSub
  i22 = sum(m22 * m22)/nSub
  i12 = sum(m11 * m22)/nSub
  #print(c(i11, i22, i12))
  
  iOut = (i11/i22)^(1/(2*p + 2)) * (2*i11 + sqrt(i11 / i22) * i12)
  
  return(c(iOut, i11/i22))
}

intCalcKR = function(m11, m22, nSub)
{
  i11 = sum(m11 * m11)/nSub
  i22 = sum(m22 * m22)/nSub
  i12 = sum(m11 * m22)/nSub
  #print(c(i11, i22, i12))
  
  iOut = (i11/i22)^0.75 * (sqrt(i11 * i22) + i12)
  
  return(c(iOut, i11/i22))
}

#-------------------------Kernel property calculation--------------------------#

#' @export
kernelPropLP = function(kernelFcn, p, drv, nInt = 5000)
{
  uSeq  = seq(from = -1, to = 1, length.out = (2 * nInt + 1))
  npMatrix = npMatrix(kernelFcn, p, nInt)
  addWeights = mWeights(npMatrix, uSeq, drv)
  
  valR  = sum((addWeights^2 * kernelFcn_use(uSeq, q = 1, kernelFcn))^2) / nInt
  valMu = sum((addWeights * kernelFcn_use(uSeq, q = 1, kernelFcn)) *
                uSeq^(p + 1)) / (nInt * factorial(p + 1))
  
  return(list(R = valR, mu = valMu))
}

kernelPropKR = function(kernelFcn, nInt = 5000)
{
  uSeq  = seq(from = -1, to = 1, length.out = (2 * nInt + 1))
  valR  = sum((kernelFcn_use(uSeq, q = 1, kernelFcn))^2) / nInt
  valMu = sum((kernelFcn_use(uSeq, q = 1, kernelFcn)) * uSeq^2) / nInt
  
  return(list(R = valR, mu = valMu))
}

#-----------------bandwidth inflation for derivative estimations---------------#

inflationFcnLP = function(h, n, dcsOptions)
{
  hInflxx = c(dcsOptions$inflPar[1] * h[1]^dcsOptions$inflExp[1],
              dcsOptions$inflPar[2] * h[2]^dcsOptions$inflExp[2])
  hInfltt = c(dcsOptions$inflPar[2] * h[1]^dcsOptions$inflExp[2],
              dcsOptions$inflPar[1] * h[2]^dcsOptions$inflExp[1])
  
  return(list(h_xx = hInflxx, h_tt = hInfltt))
}

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