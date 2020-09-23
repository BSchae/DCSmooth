################################################################################
#                                                                              #
#                DCSmooth Package: General Auxiliary Funcions                  #
#                                                                              #
################################################################################

#------------------------Set Options via Function------------------------------#

.setOptions = function(    # inside function with default values in arguments
  kernPar   = c(2, 2),     # choose a kernel function with mu = 2, 
  pOrder    = 1,           # choose order of polynomials for X-/T-smoothing
  # if pOrder == 0, kernel regression will be used.
  # (orders have to be the same in both directions)
  inflExp   = c(0.5, 0.5), # inflation exponent
  inflPar   = c(2, 1),     # inflation parameters c (regression),
  # d (2nd derivative)
  delta     = c(0.0, 0.0), # parameter for shrinking the derivatives
  constWindow = FALSE,
  varEst    = "iid"
)
{
  return(list(kernPar = kernPar, pOrder = pOrder,
              inflExp = inflExp, inflPar = inflPar, delta = delta, 
              constWindow = constWindow, varEst = varEst))
}

#--------------------Function for color ramping for plots----------------------#

plotDCScolFcn = function(Y, myColor, nColor = 100)
{
  YCol = Y - min(Y)
  YCol = YCol/max(YCol) * (nColor - 1) + 1
  colorFcn = grDevices::colorRampPalette(colors = myColor)
  colorGrad = colorFcn(nColor)
  col = matrix(colorGrad[YCol])
  return(col)
}