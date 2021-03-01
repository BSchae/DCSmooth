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
  inflPar   = c(3, 1),     # inflation parameters c (regression),
  # d (2nd derivative)
  delta     = c(0.0, 0.0), # parameter for shrinking the derivatives
  constWindow = FALSE,
  varEst    = "iid",
  modelOrder = list(ar = c(1, 1), ma = c(1, 1))
)
{
  return(list(kernPar = kernPar, pOrder = pOrder,
              inflExp = inflExp, inflPar = inflPar, delta = delta, 
              constWindow = constWindow, varEst = varEst,
              modelOrder = modelOrder))
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

#---------------------------Function for 3d plots------------------------------#

# Function should be a S3 method, not working currently

.persp3ddcs = function(DCSobj, X = 1, T = 1, color = c("#000000", "#FFFFFF"),
                    xlab = "X", ylab = "T", zlab = "Y")
{
  if (length(X) == 1) { X = DCSobj$X }
  if (length(T) == 1) { T = DCSobj$T }
  M = DCSobj$M

  colValue = plotDCScolFcn(DCSobj$M, myColor = color)

  rgl::par3d(windowRect = c(20, 30, 600, 600), viewport = c(50, 50, 500, 500))
  rgl::view3d(userMatrix = rgl::rotationMatrix(-1.4, 14, -4.3, -5), fov = 20)
  rgl::persp3d(X, T, M, col = colValue, xlab = xlab, ylab = ylab, zlab = zlab)
}