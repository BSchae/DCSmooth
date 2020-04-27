# library(mvtnorm)
# library(rgl)
# library(MASS)
# 
# ### General settings
# n_sim = 100
# seed = 42
# 
# #-------------------------------------------------------------------------------#
# #                           Sine-Function                                       #
# #-------------------------------------------------------------------------------#
# 
# ### Set up Function
# n_x = 500
# n_t = 500
# sigma_sq = 0.5
# 
# X = seq(from = 0, to = 1, length.out = n_x)
# T = seq(from = 0, to = 1, length.out = n_t)
# 
# Y0 = matrix(NA, nrow = n_x, ncol = n_t)
# 
# for (i in 1:n_x) {
#   for (j in 1:n_t) {
#     Y0[i, j] =   sin(3*pi*X[i]) * sin(2*pi*T[j]) * exp(-2*X[i]*T[j])
#   }
# }
# 
# ### Smoothing
# K_sine  = matrix(NA, nrow = n_sim, ncol = 3)
# set.seed(seed)
# 
# t0 = Sys.time()
# for(i in 1:n_sim) {
#   Y = Y0 + matrix(sqrt(sigma_sq)*rnorm(n_x*n_t), nrow = n_x, ncol = n_t)
#   Smooth = DCSmooth(Y, dcsOptions = setOptions(inflPar = c(1, 1), inflExp = 0.5))
#   K_sine[i, 1:2] = Smooth$bndw
#   K_sine[i, 3]   = mean((Smooth$M - Y0)^2)
#   print(i)
# }
# t1 = Sys.time()
# t1 - t0
# 
# plotDCS(Smooth)
# 
# #-------------------------------------------------------------------------------#
# #         Normal Distribution with single Peak (Herrmann et al. 1995)           #
# #-------------------------------------------------------------------------------#
# 
# ### Set up Function
# n_x = 400
# n_t = 250
# sigma_sq = 0.1
# 
# m1 = c(0.5, 0.5)
# s1 = 0.2*diag(2)
# 
# X = seq(from = 0, to = 1, length.out = n_x)
# T = seq(from = 0, to = 1, length.out = n_t)
# XT_combine = expand.grid(X, T)
# 
# Y0 = dmvnorm(XT_combine, m1, s1)
# Y0 = matrix(Y0, nrow = n_x, ncol = n_t)
# 
# ### Smoothing
# K_normal  = matrix(NA, nrow = n_sim, ncol = 4)
# set.seed(seed)
# 
# for(i in 1:n_sim) {
#   Y = Y0 + matrix(sqrt(sigma_sq)*rnorm(n_x*n_t), nrow = n_x, ncol = n_t)
#   Smooth = DCSmooth(Y, dcsOptions = setOptions(inflPar = c(1.5, 0.5)))
#   K_normal[i, 1:3] = Smooth$bndw
#   K_normal[i, 4] = mean((Smooth$M - Y0)^2)
#   print(i)
# }
# 
# par(mfrow = c(2, 2))
# truehist(K_normal[, 1])
# abline(v = 0.1065, col = "red")
# truehist(K_normal[, 2])
# abline(v = 0.1065, col = "red")
# truehist(K_normal[, 3])
# truehist(K_normal[, 4])
# 
# #-------------------------------------------------------------------------------#
# #                           Saddle Function                                     #
# #-------------------------------------------------------------------------------#
# 
# ### Set up Function
# Y.fkt = function(x, t) { 0.5*(t^3 - t) - (x^3 - x) }
# 
# sigma_sq = 0.01
# n_x = 1500
# n_t = 510
# 
# X = 0:(n_x - 1)/(n_x - 1)
# T = 0:(n_t - 1)/(n_t - 1)
# 
# Y0 = matrix(NA, nrow = n_x, ncol = n_t)
# 
# for (i in 1:n_x) {
#   for (j in 1:n_t) {
#     Y0[i, j] = Y.fkt(X[i], T[j])
#   }
# }
# 
# ### Smoothing
# K_saddle = matrix(NA, nrow = n_sim, ncol = 4)
# set.seed(seed)
# 
# for (i in 1:n_sim) {
#   Y = Y0 + sqrt(sigma_sq) * matrix(rnorm(n_x*n_t), nrow = n_x, ncol = n_t)
#   DCS_out = DCSmooth(Y)
#   K_saddle[i, 1:3] = DCS_out$bndw
#   K_saddle[i, 4] = mean((DCS_out$est - Y0)^2)
#   print(i)
# }
# 
# #-----------------------------------------------------------------------------#
# #                       General sine function                                 #
# #-----------------------------------------------------------------------------#
# 
# ### functions needed
# 
# sineMat = function(a, b, c, nX, nT)
# {
#   X = seq(from = 0, to = 1, length.out = nX)
#   T = seq(from = 0, to = 1, length.out = nT)
#   YOut = matrix(NA, nrow = nX, ncol = nT)
#   for (i in 1:nX)
#   {
#     for (j in 1:nT)
#     {
#       YOut[i, j] = sin(a*pi*X[i]) * sin(b*pi*T[j]) * exp(-c*X[i]*T[j])
#     }
#   }
# 
#   return(YOut)
# }
# 
# iFunc = function(a, b, c)
# {
#   n = 1000
#   X = seq(from = 0, to = 1, length.out = n)
#   T = seq(from = 0, to = 1, length.out = n)
#   mxx = matrix(NA, nrow = n, ncol = n)
#   mtt = matrix(NA, nrow = n, ncol = n)
#   for (i in 1:n)
#   {
#     for (j in 1:n)
#     {
#       x = X[i]; t = T[j]
#       mxx[i, j] = exp(-c*x*t) * sin(b*pi*t) * ((c^2*t^2 - a^2*pi^2) * sin(a*pi*x)
#                                                - 2*a*pi*c*t * cos(a*pi*x))
#       mtt[i, j] = exp(-c*x*t) * sin(a*pi*x) * ((c^2*x^2 - b^2*pi^2) * sin(b*pi*t)
#                                                - 2*b*pi*c*x * cos(b*pi*t))
#     }
#   }
# 
#   Ixx = sum(mxx^2)/n^2
#   Itt = sum(mtt^2)/n^2
#   Ixt = sum(mxx*mtt)/n^2
# 
#   IOut = (Ixx/Itt)^0.75 * (sqrt(Ixx*Itt) + Ixt)
# 
#   return(c(IOut, Ixx/Itt, Ixx, Itt, Ixt))
# }
# 
# hFunc = function(a, b, c, nX, nT, sd)
# {
#   IntFunc = iFunc(a, b, c)
#   upper = (0.71429^2*sd^2)
#   lower = c(2*nX*nT*0.14286^2 * IntFunc[1])
#   hx = (upper/lower)^(1/6)
#   ht = IntFunc[2]^0.25 * hx
# 
#   return(c(hx, ht))
# }
# 
# ### simulation study
# 
# K_sineMult = matrix(NA, nrow = n_sim, ncol = 17)
# colnames(K_sineMult) = c("hX", "hT", "hx", "ht", "iter",
#                          "a", "b", "c", "nX", "nT", "sd", "Ixx", "Itt", "Ixt",
#                          "IXX", "ITT", "IXT")
# set.seed(seed)
# 
# for (k in 1:n_sim)
# {
#   print(k)
#   # draw random values for parameters
#   a  = sample(1:100/20, size = 1)
#   b  = sample(1:100/20, size = 1)
#   c  = runif(n = 1, min = 0, max = 2)
#   nX = sample(100:600, size = 1)
#   nT = sample(100:600, size = 1)
#   sd = runif(1)
# 
#   K_sineMult[k, 6:11] = c(a, b, c, nX, nT, sd)
# 
#   K_sineMult[k, 1:2] = hFunc(a, b, c, nX, nT, sd)
# 
#   Y = sineMat(a, b, c, nX, nT) + matrix(rnorm(nX*nT)*sd, nrow = nX, ncol = nT)
#   Smooth = DCSmooth(Y)
#   K_sineMult[k, 3:5] = Smooth$bndw[c(1, 2, 6)]
#   K_sineMult[k, 12:14] = Smooth$bndw[3:5]
#   K_sineMult[k, 15:17] = iFunc(a, b, c)[3:5]
# }
# 
# ### evaluation
# K = data.frame(K_sineMult)
# 
# delta_hx = K$hx - K$hX
# delta_ht = K$ht - K$hT
# 
# par(mfrow = c(2, 2))
# truehist(delta_hx)
# truehist(delta_ht)
# 
# plot(K$hX, K$hx)
# abline(a = 0, b = 1, col = 2)
# abline(lm(K$hx~K$hX), col = 4)
# plot(K$hT, K$ht)
# abline(a = 0, b = 1, col = 2)
# abline(lm(K$ht~K$hT), col = 4)
# 
# #-----------------------------------------------------------------------------#
# #                           Polynomial Funcion                                #
# #-----------------------------------------------------------------------------#
# 
# ### functions
# 
# polyFcn = function(x, t, pX, pT, coefMat)
# {
#   xVec = x^(0:pX)
#   tVec = t^(0:pT)
#   value = xVec %*% coefMat %*% tVec
#   return(value)
# }
# 
# polyDrv20 = function(x, t, pX, pT, coefMat)
# {
#   xVec = c(0, 0, (2:pX)*(1:(pX - 1))*x^(0:(pX - 2)))
#   tVec = t^(0:pT)
#   value = xVec %*% coefMat %*% tVec
#   return(value)
# }
# 
# polyDrv02 = function(x, t, pX, pT, coefMat)
# {
#   xVec = x^(0:pX)
#   tVec = c(0, 0, (2:pT)*(1:(pT - 1))*t^(0:(pT - 2)))
#   value = xVec %*% coefMat %*% tVec
#   return(value)
# }
# 
# intPolyFunc = function(pX, pT, coefMat)
# {
#   nInt = 500
#   X = T = seq(from = 0, to = 1, length.out = nInt)
#   mxx = mtt = matrix(NA, nInt, nInt)
#   for (i in 1:nInt)
#   {
#     for (j in 1:nInt)
#     {
#       mxx[i, j] = polyDrv20(X[i], T[j], pX, pT, coefMat)
#       mtt[i, j] = polyDrv02(X[i], T[j], pX, pT, coefMat)
#     }
#   }
# 
#   Ixx = sum(mxx^2)/nInt^2
#   Itt = sum(mtt^2)/nInt^2
#   Ixt = sum(mxx*mtt)/nInt^2
# 
#   IOut = (Ixx/Itt)^0.75 * (sqrt(Ixx*Itt) + Ixt)
# 
#   return(c(IOut, Ixx/Itt, Ixx, Itt, Ixt))
# }
# 
# hPolyFunc = function(pX, pT, coefMat, nX, nT, sd)
# {
#   IntFunc = intPolyFunc(pX, pT, coefMat)
#   upper = (0.71429^2*sd^2)
#   lower = c(2*nX*nT*0.14286^2 * IntFunc[1])
#   hx = (upper/lower)^(1/6)
#   ht = IntFunc[2]^0.25 * hx
# 
#   return(c(hx, ht))
# }
# 
# ### Simulation Study
# 
# K_PolyMult = matrix(NA, n_sim, 16)
# colnames(K_PolyMult) = c("hX", "hT", "hx", "ht", "iter",
#                          "nX", "nT", "pX", "pT", "sd", "Ixx", "Itt", "Ixt",
#                          "IXX", "ITT", "IXT")
# set.seed(seed)
# myOpt = setOptions(inflPar = c(1, 0.5), inflExp = c(0.3, 0.5))
# 
# for (s in 1:n_sim)
# {
#   print(s)
# 
#   # draw parameters
#   pX = sample(2:8, 1); pT = sample(2:8, 1)
#   nX = sample(100:600, 1); nT = sample(100:600, 1)
#   coefN = (pX + 1)*(pT + 1)
#   coefMat = matrix(rnorm(coefN), pX + 1, pT + 1)
#   sd = 2*runif(1)
#   X = seq(from = 0, to = 1, length.out = nX)
#   T = seq(from = 0, to = 1, length.out = nT)
#   Y = matrix(NA, nrow = nX, ncol = nT)
# 
#   # simulate function
#   for (i in 1:nX)
#   {
#     for (j in 1:nT)
#     {
#       Y[i, j] = polyFcn(X[i], T[j], pX, pT, coefMat) + rnorm(1)*sd
#     }
#   }
# 
#   # true bandwidths
#   K_PolyMult[s, 1:2] = hPolyFunc(pX, pT, coefMat, nX, nT, sd)
# 
#   Smooth = DCSmooth(Y)
#   K_PolyMult[s, 3:5] = Smooth$bndw[c(1, 2, 6)]
#   K_PolyMult[s, 6:7] = c(nX, nT)
#   K_PolyMult[s, 8:9] = c(pX, pT)
#   K_PolyMult[s, 10]  = sd
#   K_PolyMult[s, 11:13] = Smooth$bndw[3:5]
#   K_PolyMult[s, 14:16] = intPolyFunc(pX, pT, coefMat)[3:5]
# }
# 
# ## evaluation
# 
# K = data.frame(K_PolyMult[, ])
# 
# #write.csv(K, file = "KPoly_allInt.csv")
# #K = read.csv("KPoly_all.csv")
# 
# delta_hx = K$hx - K$hX
# delta_ht = K$ht - K$hT
# 
# Iout = (K$Ixx/K$Itt)^0.75 * (sqrt(K$Ixx*K$Itt) + K$Ixt)
# IOUT = (K$IXX/K$ITT)^0.75 * (sqrt(K$IXX*K$ITT) + K$IXT)
# 
# # Picture 1
# par(mfrow = c(2, 2))
# truehist(delta_hx)
# truehist(delta_ht)
# plot(K$hX, K$hx)
# abline(a = 0, b = 1, col = 2)
# a = lm(K$hx~K$hX)
# abline(a, col = 4)
# plot(K$hT, K$ht)
# abline(a = 0, b = 1, col = 2)
# a = lm(K$ht~K$hT)
# abline(a, col = 4)
# 
# # picture 2
# par(mfrow = c(2, 2))
# plot(K$IXX, K$Ixx)
# abline(a = 0, b = 1, col = 2)
# abline(lm(K$Ixx~K$IXX), col = 4)
# plot(K$ITT, K$Itt)
# abline(a = 0, b = 1, col = 2)
# abline(lm(K$Itt~K$ITT), col = 4)
# plot(K$IXT, K$Ixt)
# abline(a = 0, b = 1, col = 2)
# abline(lm(K$Ixt~K$IXT), col = 4)
# plot(IOUT, Iout)
# abline(a = 0, b = 1, col = 2)
# abline(lm(Iout~IOUT), col = 4)
# 
# # picture 2
# par(mfrow = c(2, 2))
# plot(K$IXX, K$Ixx, ylim = c(0, 500), xlim = c(0, 300))
# abline(a = 0, b = 1, col = 2)
# abline(lm(K$Ixx~K$IXX), col = 4)
# plot(K$ITT, K$Itt, ylim = c(0, 500), xlim = c(0, 300))
# abline(a = 0, b = 1, col = 2)
# abline(lm(K$Itt~K$ITT), col = 4)
# plot(K$IXT, K$Ixt, ylim = c(-200, 500), xlim = c(-200, 500))
# abline(a = 0, b = 1, col = 2)
# abline(lm(K$Ixt~K$IXT), col = 4)
# plot(IOUT, Iout, ylim = c(0, 200), xlim = c(0, 200))
# abline(a = 0, b = 1, col = 2)
# abline(lm(Iout~IOUT), col = 4)
# 
# plot(K$Ixx, K$Itt)
# 
# summary(lm(K$Ixx~K$IXX))
# names(K)
