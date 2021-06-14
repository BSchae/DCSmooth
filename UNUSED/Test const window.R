library(rgl)
Y = y.norm1 + rnorm(101^2)
dcs_options = set.options()
kernel_0 = kernel_fcn_assign("MW220")
kernel_2 = kernel_fcn_assign("MW422")

Ysmth = LP_dcs_const0(Y, c(0.2, 0.2), c(1, 1), c(0, 0), kernel_0, kernel_0)

Ytest = Ysmth; #Ytest[30:72, 20:82] = Y[30:72, 20:82]
a = LP_dcs_const0(Y, c(0.55, 0.35), c(3, 1), c(2, 0), kernel_0, kernel_0)
b = LP_dcs_const0(Ytest, c(0.35, 0.55), c(3, 1), c(2, 0), kernel_0, kernel_0)

sum(a^2)/101^2
sum(b^2)/101^2
sum(m20val^2)/101^2

persp3d(a, color = "coral1")
persp3d(b, color = "steelblue")
persp3d(m20val, color = "aquamarine")

plot(m20val[50, ], type = "l", ylim = range(a[50, ], b[50, ], m20val[50, ]))
lines(a[50, ], col = 2)
lines(b[50, ], col = 4)
plot(m20val[, 50], type = "l", ylim = range(a[, 50], b[, 50], m20val[, 50]))
lines(a[, 50], col = 2)
lines(b[, 50], col = 4)

sum(a^2)/101^2
sum(b^2)/101^2

mean(((a^2 - m20val^2)/101^2)^2)
mean(((b^2 - m20val^2)/101^2)^2)
mean(((c^2 - m20val^2)/101^2)^2)
mean(((d^2 - m20val^2)/101^2)^2)

persp3d(m20val, col = 2)


#------------------------------------------------------------------------------#

nSim = 100
h = 0.3
p = 3
drv = 0
kernFcn = kernFkt_MW220

x = 0:(nSim - 1)/(nSim - 1)*h*2
y = -(x - 0.2)^2 + 0.01*rnorm(nSim)
#y = 0.2*x + 0.01*rnorm(nSim)
m1 = m2 = m3 = m4 = qvec = 1:trunc(nSim/2)*NA
windowWidth = nSim


for (i in 1:(trunc(nSim/2)))
{
  xi = x[1:windowWidth] - x[i]
  h_window = 2*h - x[i]
  q = 2*h/h_window - 1
  
  u = (x[i] - x[1:windowWidth])/h_window
  WK = as.vector(kernFcn(u, 1))
  
  xMat = t(t(matrix(rep(xi, times = p + 1), ncol = p + 1))^(0:p))
  weights = (solve(t(WK*xMat) %*% xMat) %*% t(WK*xMat))[drv + 1, ]
  
  m1[i] = (weights %*% y[1:nSim]) * factorial(drv)
}

for (i in 1:(trunc(nSim/2)))
{
  xi = x[1:(i + trunc(nSim/2) - 1)] - x[i]
  u = - xi / h
  WK = as.vector(kernFcn(u, 1))
  
  xMat = t(t(matrix(rep(xi, times = p + 1), ncol = p + 1))^(0:p))
  weights = (solve(t(WK*xMat) %*% xMat) %*% t(WK*xMat))[drv + 1, ]
  
  m2[i] = (weights %*% y[1:length(weights)]) * factorial(drv)
}
for (i in 1:(trunc(nSim/2)))
{
  xi = x[1:windowWidth] - x[i]
  h_window = 2*h - x[i]
  q = 2*h/h_window - 1
  
  u = (x[i] - x[1:windowWidth])/h_window
  WK = as.vector(kernFcn(u, q))
  
  xMat = t(t(matrix(rep(xi, times = p + 1), ncol = p + 1))^(0:p))
  weights = (solve(t(WK*xMat) %*% xMat) %*% t(WK*xMat))[drv + 1, ]
  
  m3[i] = (weights %*% y[1:nSim]) * factorial(drv)
}

for (i in 1:(trunc(nSim/2)))
{
  xi = x[1:(i + trunc(nSim/2) - 1)] - x[i]
  q = (i - 1) / (trunc(nSim/2) - 1)
  WK = as.vector(kernFcn(u, q))
  qvec[i] = q
  
  xMat = t(t(matrix(rep(xi, times = p + 1), ncol = p + 1))^(0:p))
  weights = (solve(t(WK*xMat) %*% xMat) %*% t(WK*xMat))[drv + 1, ]
  
  m4[i] = (weights %*% y[1:length(weights)]) * factorial(drv)
}


library(tidyverse)
m_df = data.frame(x = x[1:trunc(nSim/2)], y = y[1:trunc(nSim/2)], m1, m2, m3, m4)
m_df = pivot_longer(m_df, cols = 3:6, names_to = "type", values_to = "m")
m_df$type = factor(m_df$type, labels = c("const_1", "const_0", "const_1_bmod",
                                         "const_0_bmod"))
ggplot(m_df) + geom_point(aes(x = x, y = y)) +
  geom_line(aes(x = x, y = m, color = type))


plot(x[1:trunc(nSim/2)], m1, type = "l", ylim = range(m1, m2, m3, m4))
lines(x[1:trunc(nSim/2)], m2, col = 2)
lines(x[1:trunc(nSim/2)], m3, col = 4)
lines(x[1:trunc(nSim/2)], m4, col = 6)
plot(diff(m1))


#------------------------------------------------------------------------------#


#X = 1:1000/1000
X = 0:100/100
# Y0 = y.norm3[, 50]
# Y2 = m20val[, 50]
Y0 = 1/sqrt(2*pi*0.2^2)*exp(-(X - 0.5)^2/(2*0.2^2))
Y2 = 1/sqrt(2*pi*0.2^2)*exp(-(X - 0.5)^2/(2*0.2^2)) * (625*X^2 - 625*X + 131.25)
#Y0 = 0.5*X^3 + 1*X^2 - 2*X
#Y1 = 1.5*X^2 + 2*X - 2
#Y2 = 3*X + 2
set.seed(123)
Y = Y0 + 0.05*rnorm(length(Y0))
h = 0.3
# i = 1

#plot(X, Y, xlim = c(0.1, 0.3))
plot(X, Y)
lines(X, Y0, col = 2)
lines(X, Y1, col = 4)
lines(X, Y2, col = 6)
test0 = LPSmooth(X, Y, 0, kernFkt_MW220, 0.2)
test1a = LPSmooth(X, Y, 1, kernFkt_MW220, 0.2)
test1b = LPSmooth(X, test0, 1, kernFkt_MW220, 0.07)
test2a = LPSmooth(X, Y, 2, kernFkt_MW220, 0.1)
test2b = LPSmooth(X, test1b, 1, kernFkt_MW220, 0.07)
test2c = LPSmooth(X, test0, 2, kernFkt_MW220, 0.07)
testY = test0; testY[11:92] = Y[11:92]
test2d = LPSmooth(X, testY, 2, kernFkt_MW220, 0.1)
# test3 = LPSmooth_R(X, Y, 3, 2, 3, 0.4)$Ye.LP
lines(X, test0, col = 2, lty = 2)
lines(X, test1a, col = 4, lty = 2)
lines(X, test1b, col = 4, lty = 4)
plot(X, Y2, type = "l", ylim = range(Y2, test2b, test2c))
lines(X, test2a, col = 2, lty = 2)
lines(X, test2b, col = 4, lty = 3)
lines(X, test2c, col = 5, lty = 3)
lines(X, test2d, col = 3, lty = 3)

kernel_0 = kernel_fcn_assign("MW220")
kernel_2 = kernel_fcn_assign("MW422")

Y = as.matrix(t(Y[1:101]))
a = LPSmooth_matrix_BMod(Y, 0.3, 3, 2, kernel_0)
b = LPSmooth_matrix2_BMod(Y, 0.3, 3, 2, kernel_0)
c = LPSmooth_cw(X, as.vector(Y), 2, k2, 0.3)
#d = LPSmooth(X, Y, 2, kernFkt_MW220, 0.35)
# c = LPSmooth_matrix_BMod(Y, 0.3, 3, 2, kernel_2)
# d = LPSmooth_matrix2_BMod(Y, 0.3, 3, 2, kernel_2)
plot(X, Y2, type = "l", lty = 2, ylim = range(a))
lines(X, a, col = 2)
lines(X, b, col = 4)
lines(X, c, col = 5)
#lines(X, d, col = 6)
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
    WK = -as.vector(kernFcn(u, q))
    
    xMat = t(t(matrix(rep(Xi, times = p + 1), ncol = p + 1))^(0:p))
    weights = (solve(t(WK*xMat) %*% xMat) %*% t(WK*xMat))[drv + 1, ]
    #lines(weights, col = 0, lty = 2)
    
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
  
  # # boundaries
  # for (i in (n - bndw + 1):n)
  # {
  #   Xi = X[(n - windowWidth + 1):n] - X[i]
  #   h_window = 2*h - (1 - X[i])
  #   q = 2*h/h_window - 1
  #   
  #   u = (X[i] - X[(n - windowWidth + 1):n])/h_window
  #   WK = as.vector(kernFcn(-u, 1))
  #   
  #   xMat = t(t(matrix(rep(Xi, times = p + 1), ncol = p + 1))^(0:p))
  #   weights = (solve(t(WK*xMat) %*% xMat) %*% t(WK*xMat))[drv + 1, ]
  #   lines(weights, col = 2, lty = 2)
  #   
  #   yEst[i] = (weights %*% Y[(n - windowWidth + 1):n]) * factorial(drv)
  #   yEst[n - i + 1] = (rev(weights) %*% Y[(n - min(length(Xi), n) + 1):n]) * factorial(drv) * (-1)^drv
  # }
  
  
  
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

#------------------------------------------------------------------------------#

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

