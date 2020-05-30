###############################################################################
#                                                                             #
#       DCSmooth Package: Auxiliary Funcions for Bandwidth Selection          #
#                                                                             #
###############################################################################

#-----------------------Set Options via Function------------------------------#

.setOptions = function(    # inside function with default values in arguments
  kernPar   = c(2, 2),     # choose a kernel function with mu = 2, 
  pOrder    = 0,           # choose order of polynomials for X-/T-smoothing
                           # if pOrder == 0, kernel regression will be used.
  # (orders have to be the same in both directions)
  inflExp   = c(0.5, 0.5), # inflation exponent
  inflPar   = c(2, 1), # inflation parameters c (regression), d (2nd derivative)
  fast      = FALSE
)
{
  return(list(kernPar = kernPar, pOrder = pOrder,
              inflExp = inflExp, inflPar = inflPar, fast = fast))
}

#----------------------Formula for optimal bandwidths-------------------------#

hOptFunction = function(mxx, mtt, varCoef, n, nSub, p, kernelProp)
{
  i0x = intCalc(mxx, mtt, nSub)
  i0t = intCalc(mtt, mxx, nSub)
  
  hxOpt = (kernelProp$R^2 * varCoef)/((p + 1) * n * kernelProp$mu^2 * i0x)
  htOpt = (kernelProp$R^2 * varCoef)/((p + 1) * n * kernelProp$mu^2 * i0t)
  
  hxOpt = hxOpt^(1/(2*p + 4))
  htOpt = htOpt^(1/(2*p + 4))
  
  return(c(hxOpt, htOpt))
}

hOptKR = function(mxx, mtt, varCoef, n, nSub, kernelProp)
{
  i0x = intCalc(mxx, mtt, nSub)
  i0t = intCalc(mtt, mxx, nSub)
  
  hxOpt = (kernelProp$R^2 * varCoef)/(n * kernelProp$mu^2 * i0x)
  htOpt = (kernelProp$R^2 * varCoef)/(n * kernelProp$mu^2 * i0t)
  
  hxOpt = hxOpt^(1/6)
  htOpt = htOpt^(1/6)
  
  return(c(hxOpt, htOpt))
}

#---------------------Integrals over mxx^2, mtt^2-----------------------------#

intCalc = function(m11, m22, nSub)
{
  i11 = sum(m11 * m11)/nSub
  i22 = sum(m22 * m22)/nSub
  i12 = sum(m11 * m22)/nSub
  
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

#-------------------Function for plotting smoothed surface--------------------#

.plotDCS = function(DCSobj, X = 1, T = 1, color = c("#000000", "#FFFFFF"),
                    xlab = "X", ylab = "T", zlab = "Y")
{
  if (length(X) == 1) { X = DCSobj$X }
  if (length(T) == 1) { T = DCSobj$T }
  M = DCSobj$M
  
  colValue = plotDCScolFcn(DCSobj$M, myColor = color)

  par3d(windowRect = c(20, 30, 600, 600), viewport = c(50, 50, 500, 500))
  view3d(userMatrix = rotationMatrix(-1.4, 14, -4.3, -5), fov = 20)
  persp3d(X, T, M, col = colValue, xlab = xlab, ylab = ylab, zlab = zlab)
}

plotDCScolFcn = function(Y, myColor, nColor = 100)
{
  YCol = Y - min(Y)
  YCol = YCol/max(YCol) * (nColor - 1) + 1
  colorFcn = colorRampPalette(colors = myColor)
  colorGrad = colorFcn(nColor)
  col = matrix(colorGrad[YCol])
  return(col)
}

#-------------------------Estimation of cf (temporary)------------------------#

DCS.cf = function(Y_0, X, T, order.arma) {
  n_x = length(X)
  n_t = length(T)
  
  Y_vec_x = c(Y_0)
  Y_vec_t = c(t(Y_0))
  
  # Estimation of ARMA in X direction
  BIC_x = matrix(NA, nrow = (order.arma[1] + 1), ncol = (order.arma[3] + 1))
  for (p in 0:order.arma[1]) {
    for (q in 0:order.arma[3]) {
      arma.est = arima(Y_vec_x, order = c(p, 0, q))
      BIC_x[p + 1, q + 1] = -2 * arma.est$loglik + log(n_x*n_t) * (p + q)
    }
  }
  order.opt_x = which(BIC_x == min(BIC_x), arr.ind = TRUE) - 1
  ARMA_x = arima(Y_vec_x, order = c(order.opt_x[1], 0, order.opt_x[2]))
  
  # Store sum of coefficients
  sc.AR_x = ifelse(order.opt_x[1] != 0, sum(ARMA_x$coef[1:order.opt_x[1]]), 0)
  sc.MA_x = ifelse(order.opt_x[2] != 0, sum(ARMA_x$coef[(order.opt_x[1] + 1):
                                                          (order.opt_x[1] + order.opt_x[2])]), 0)
  
  # Estimation of ARMA in T direction
  BIC_t = matrix(NA, nrow = (order.arma[2] + 1), ncol = (order.arma[4] + 1))
  for (p in 0:order.arma[2]) {
    for (q in 0:order.arma[4]) {
      arma.est = arima(Y_vec_t, order = c(p, 0, q))
      BIC_t[p + 1, q + 1] = -2 * arma.est$loglik + log(n_x*n_t) * (p + q)
    }
  }
  order.opt_t = which(BIC_t == min(BIC_t), arr.ind = TRUE) - 1
  ARMA_t = arima(Y_vec_t, order = c(order.opt_t[1], 0, order.opt_t[2]))
  
  sc.AR_t = ifelse(order.opt_t[1] != 0, sum(ARMA_t$coef[1:order.opt_t[1]]), 0)
  sc.MA_t = ifelse(order.opt_t[2] != 0, sum(ARMA_t$coef[(order.opt_t[1] + 1):
                                                          (order.opt_t[1] + order.opt_t[2])]), 0)
  
  fraction_MA_AR = (((1 + sc.MA_x)*(1 + sc.MA_t)) / ((1 - sc.AR_x)*(1 - sc.AR_t)))
  cf_out = fraction_MA_AR^2 * sqrt(ARMA_x$sigma2 * ARMA_t$sigma2)
  return(list(cf_out = cf_out, arma_x = order.opt_x, arma_t = order.opt_t))
}
