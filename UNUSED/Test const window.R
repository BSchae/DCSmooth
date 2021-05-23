library(rgl)
Y = y.norm3 + 0.5*rnorm(101^2)
dcs_options = set.options()
kernel_0 = kernel_fcn_assign("MW220")
kernel_2 = kernel_fcn_assign("MW422")

a = LP_dcs_const1(Y, c(0.4, 0.2), c(3, 1), c(2, 0), kernel_0, kernel_0)
b = LP_dcs_const0(Y, c(0.4, 0.2), c(3, 1), c(2, 0), kernel_0)
c = LP_dcs_const1_BMod(Y, c(0.4, 0.2), c(3, 1), c(2, 0), kernel_2, kernel_0)
d = LP_dcs_const0_BMod(Y, c(0.4, 0.2), c(3, 1), c(2, 0), kernel_2, kernel_0)

mean(((a^2 - m20val^2)/101^2)^2)
mean(((b^2 - m20val^2)/101^2)^2)
mean(((c^2 - m20val^2)/101^2)^2)
mean(((d^2 - m20val^2)/101^2)^2)

persp3d(m20val, col = 2)






#X = 1:1000/1000
X = 0:100/100
Y0 = y.norm1[50, ]
Y2 = m20val[50, ]
# Y0 = 1/sqrt(2*pi*0.2^2)*exp(-(X - 0.5)^2/(2*0.2^2))
# Y2 = 1/sqrt(2*pi*0.2^2)*exp(-(X - 0.5)^2/(2*0.2^2)) * (625*X^2 - 625*X + 131.25)
# Y0 = 0.05*X^3 + 0.2*X^2 - 0.1*X
# Y1 = 0.15*X^2 + 0.4*X - 0.1
# Y2 = 0.3*X + 0.4
set.seed(123)
Y = Y0 + 0.5*rnorm(length(Y0)) +5
h = 0.3
# i = 1

#plot(X, Y, xlim = c(0.1, 0.3))
plot(X, Y)
lines(X, Y0, col = 2, lwd = 2)
test1 = LPSmooth_cw(X, Y, 2, kernFkt_MW422, 0.4)
test2 = LPSmooth(X, Y, 2, kernFkt_MW422, 0.4)
# test3 = LPSmooth_R(X, Y, 3, 2, 3, 0.4)$Ye.LP
lines(X, Ytest, col = 4, lwd = 2)
lines(X, test2, col = 6, lwd = 2)
lines(X, test3, col = 8, lwd = 2)

Y = as.matrix(t(Y))
a = LPSmooth_matrix(Y, 0.3, 3, 2, kernel_0)
b = LPSmooth_matrix2(Y, 0.3, 3, 2, kernel_0)
c = LPSmooth_matrix_BMod(Y, 0.3, 3, 2, kernel_2)
d = LPSmooth_matrix2_BMod(Y, 0.3, 3, 2, kernel_2)
plot(X, Y2, type = "l", lty = 2, ylim = range(a, b, c, d))
lines(X, a, col = 2)
lines(X, b, col = 4)
lines(X, c, col = 5)
lines(X, d, col = 6)
mean((a - Y2)^2)
mean((b - Y2)^2)
mean((c - Y2)^2)
mean((d - Y2)^2)

LPSmooth_cw = function(X, Y, drv, kernFcn, h)
{
  p = drv + 1
  n = length(X)
  bndw = trunc(n*h)
  windowWidth = 2*bndw + 1
  yEst = 1:n*0
  
  # boundaries
  for (i in 1:bndw)
  {
    if (i > trunc(n/2) + 1)
    {
      break()
    }
    Xi = X[1:windowWidth] - X[i]
    h_window = 2*h - X[i]
    q = 2*h/h_window - 1
    
    u = (X[i] - X[1:windowWidth])/h_window
    WK = as.vector(kernFcn(u, q))
    
    xMat = t(t(matrix(rep(Xi, times = p + 1), ncol = p + 1))^(0:p))
    weights = (solve(t(WK*xMat) %*% xMat) %*% t(WK*xMat))[drv + 1, ]
    lines(weights, col = 0, lty = 2)
    
    yEst[i] = (weights %*% Y[1:min(length(Xi), n)]) * factorial(drv)
    yEst[n - i + 1] = (rev(weights) %*% Y[(n - min(length(Xi), n) + 1):n]) * factorial(drv) * (-1)^drv
  }
  
  # interior values
  if (h <= 0.5)
  {
    x = -bndw:bndw / n
    u = - x / h
    WK = as.vector(kernFcn(u, 1))
    
    xMat = t(t(matrix(rep(x, times = p + 1), ncol = p + 1))^(0:p))
    weights = (solve(t(WK*xMat) %*% xMat) %*% t(WK*xMat))[drv + 1, ]
    
    for (i in (bndw + 1):(n - bndw))
    {
      yEst[i] = (weights %*% Y[(i - bndw):(i + bndw)]) * factorial(drv)
    }
  }  
  return(yEst)
}

LPSmooth = function(X, Y, drv, kernFcn, h)
{
  p = drv + 1
  n = length(X)
  bndw = trunc(n*h)
  yEst = 1:n*0
  
  # boundaries
  for (i in 1:bndw)
  {
    x = -(i - 1):bndw / n
    u = - x / h
    q = (i - 1) / (bndw - 1)
    WK = as.vector(kernFcn(u, q))
    
    xMat = t(t(matrix(rep(x, times = p + 1), ncol = p + 1))^(0:p))
    weights = (solve(t(WK*xMat) %*% xMat) %*% t(WK*xMat))[drv + 1, ]
    
    yEst[i] = (weights %*% Y[1:length(x)]) * factorial(drv)
    yEst[n - i + 1] = (rev(weights) %*% Y[(n - length(x) + 1):n]) *
                                            factorial(drv) * (-1)^drv
  }
  
  # interior values
  x = -bndw:bndw / n
  u = - x / h
  WK = as.vector(kernFcn(u, 1))
  
  xMat = t(t(matrix(rep(x, times = p + 1), ncol = p + 1))^(0:p))
  weights = (solve(t(WK*xMat) %*% xMat) %*% t(WK*xMat))[drv + 1, ]
  
  for (i in (bndw + 1):(n - bndw))
  {
    yEst[i] = (weights %*% Y[(i - bndw):(i + bndw)]) * factorial(drv)
  }

  return(yEst)
}

n = 4000
h = 0.3
x = 0:(n - 1)/(n - 1)*2*h

p = 1
drv = 0


color = c("#444C5C", "#78A5A3", "#E1B16A", "#CE5A57")
color_fcn = grDevices::colorRampPalette(colors = color)
colors_krn = color_fcn(trunc(n/2))
kernel_function = kernFkt_MW220

c = b = a = 1:trunc(n/2)

for (i in 1868:2068)
{
  xi = x - x[i]
  h_window = 2*h - x[i]
  q = 2*h/h_window - 1
  
  u = (x[i] - x)/h_window
  
  
  a[i] = q; b[i] = u[1]; c[i] = h_window
  
  WK = as.vector(kernel_function(u, q))
  xMat = t(t(matrix(rep(xi, times = p + 1), ncol = p + 1))^(0:p))
  weights = (solve(t(WK*xMat) %*% xMat) %*% t(WK*xMat))[drv + 1, ]
  if(i == 1)
  {
    #plot(u, WK, type = "l", ylim = c(-1.5, 5), xlim = c(-1, 1))
    plot(x, weights, type = "l", ylim = c(-0.1, 0.2), xlim = c(0, 0.6))
  }
  if (i %% 10 == 0)
  {
    #lines(u, WK, col = colors_krn[i])
    lines(x, weights, col = colors_krn[(i - 1867)/10])
    #points(x, weights, col = colors_krn[i], pch = 16)
    text(x = q, y = 0, labels = i)
  }
}

i = 1

