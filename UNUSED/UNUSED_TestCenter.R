library(rgl)
library(MASS)
library(microbenchmark)

### General settings
n_sim = 500
seed = 42

#-----------------------------------------------------------------------------#
#                             Sine function                                   #
#-----------------------------------------------------------------------------#

### functions needed

sineMat = function(a, b, c, nX, nT)
{
  X = seq(from = 0, to = 1, length.out = nX)
  T = seq(from = 0, to = 1, length.out = nT)
  YOut = matrix(NA, nrow = nX, ncol = nT)
  for (i in 1:nX)
  {
    for (j in 1:nT)
    {
      YOut[i, j] = sin(a*pi*X[i]) * sin(b*pi*T[j]) * exp(-c*X[i]*T[j])
    }
  }

  return(YOut)
}

iFunc = function(a, b, c)
{
  n = 1000
  X = seq(from = 0, to = 1, length.out = n)
  T = seq(from = 0, to = 1, length.out = n)
  mxx = matrix(NA, nrow = n, ncol = n)
  mtt = matrix(NA, nrow = n, ncol = n)
  for (i in 1:n)
  {
    for (j in 1:n)
    {
      x = X[i]; t = T[j]
      mxx[i, j] = exp(-c*x*t) * sin(b*pi*t) * ((c^2*t^2 - a^2*pi^2) * sin(a*pi*x)
                                               - 2*a*pi*c*t * cos(a*pi*x))
      mtt[i, j] = exp(-c*x*t) * sin(a*pi*x) * ((c^2*x^2 - b^2*pi^2) * sin(b*pi*t)
                                               - 2*b*pi*c*x * cos(b*pi*t))
    }
  }

  Ixx = sum(mxx^2)/n^2
  Itt = sum(mtt^2)/n^2
  Ixt = sum(mxx*mtt)/n^2

  IOut = (Ixx/Itt)^0.75 * (sqrt(Ixx*Itt) + Ixt)

  return(c(IOut, Ixx/Itt, Ixx, Itt, Ixt))
}

hFunc = function(a, b, c, nX, nT, sd)
{
  IntFunc = iFunc(a, b, c)
  upper = (0.71429^2*sd^2)
  lower = c(nX*nT*0.14286^2 * IntFunc[1])
  hx = (upper/lower)^(1/6)
  ht = IntFunc[2]^0.25 * hx

  return(c(hx, ht))
}

### simulation study

K_sineMult = matrix(NA, nrow = n_sim, ncol = 12)
colnames(K_sineMult) = c("hX", "hT", "hx", "ht", "a", "b", "c", "nX", "nT",
                         "sd", "iter", "time")
set.seed(seed)
myOpt = setOptions(pOrder = 0, fast = TRUE)

for (k in 1:n_sim)
{
  # draw random values for parameters
  a  = sample(1:100/20, size = 1)
  b  = sample(1:100/20, size = 1)
  c  = runif(n = 1, min = 0, max = 2)
  nX = sample(500:2000, size = 1)
  nT = sample(500:2000, size = 1)
  sd = runif(1)

  K_sineMult[k, 5:10] = c(a, b, c, nX, nT, sd)

  hOpt = hFunc(a, b, c, nX, nT, sd)

  K_sineMult[k, 1:2] = hOpt

  Y = sineMat(a, b, c, nX, nT) + matrix(rnorm(nX*nT)*sd, nrow = nX, ncol = nT)

  t0 = Sys.time()
  Smooth = DCSmooth(Y, dcsOptions = myOpt)
  K_sineMult[k, 12] = difftime(Sys.time(), t0, units = "secs")

  K_sineMult[k, 3:4] = Smooth$bndw[1:2]
  K_sineMult[k, 11] = Smooth$iterations

  print(k)
}

K = data.frame(K_sineMult)
plot(K[, 2], K[, 4])

abline(a = 0, b= 1)
abline(lm(K[, 4]~K[, 2]))

#-----------------------------------------------------------------------------#
#                           Polynomial Funcion                                #
#-----------------------------------------------------------------------------#

### functions

polyFcn = function(x, t, pX, pT, coefMat)
{
  xVec = x^(0:pX)
  tVec = t^(0:pT)
  value = xVec %*% coefMat %*% tVec
  return(value)
}

polyDrv20 = function(x, t, pX, pT, coefMat)
{
  xVec = c(0, 0, (2:pX)*(1:(pX - 1))*x^(0:(pX - 2)))
  tVec = t^(0:pT)
  value = xVec %*% coefMat %*% tVec
  return(value)
}

polyDrv02 = function(x, t, pX, pT, coefMat)
{
  xVec = x^(0:pX)
  tVec = c(0, 0, (2:pT)*(1:(pT - 1))*t^(0:(pT - 2)))
  value = xVec %*% coefMat %*% tVec
  return(value)
}

intPolyFunc = function(pX, pT, coefMat)
{
  nInt = 500
  X = T = seq(from = 0, to = 1, length.out = nInt)
  mxx = mtt = matrix(NA, nInt, nInt)
  for (i in 1:nInt)
  {
    for (j in 1:nInt)
    {
      mxx[i, j] = polyDrv20(X[i], T[j], pX, pT, coefMat)
      mtt[i, j] = polyDrv02(X[i], T[j], pX, pT, coefMat)
    }
  }

  Ixx = sum(mxx^2)/nInt^2
  Itt = sum(mtt^2)/nInt^2
  Ixt = sum(mxx*mtt)/nInt^2

  IOut = (Ixx/Itt)^0.75 * (sqrt(Ixx*Itt) + Ixt)

  return(c(IOut, Ixx/Itt, Ixx, Itt, Ixt))
}

hPolyFunc = function(pX, pT, coefMat, nX, nT, sd)
{
  IntFunc = intPolyFunc(pX, pT, coefMat)
  upper = (0.71429^2*sd^2)
  lower = c(2*nX*nT*0.14286^2 * IntFunc[1])
  hx = (upper/lower)^(1/6)
  ht = IntFunc[2]^0.25 * hx

  return(c(hx, ht))
}

### Simulation Study

K_PolyMult = matrix(NA, n_sim, 16)
colnames(K_PolyMult) = c("hX", "hT", "hx", "ht", "iter",
                         "nX", "nT", "pX", "pT", "sd", "Ixx", "Itt", "Ixt",
                         "IXX", "ITT", "IXT")
set.seed(seed)
myOpt = setOptions(inflPar = c(1, 0.5), inflExp = c(0.3, 0.5))

for (s in 1:n_sim)
{
  print(s)

  # draw parameters
  pX = sample(2:8, 1); pT = sample(2:8, 1)
  nX = sample(100:600, 1); nT = sample(100:600, 1)
  coefN = (pX + 1)*(pT + 1)
  coefMat = matrix(rnorm(coefN), pX + 1, pT + 1)
  sd = 2*runif(1)
  X = seq(from = 0, to = 1, length.out = nX)
  T = seq(from = 0, to = 1, length.out = nT)
  Y = matrix(NA, nrow = nX, ncol = nT)

  # simulate function
  for (i in 1:nX)
  {
    for (j in 1:nT)
    {
      Y[i, j] = polyFcn(X[i], T[j], pX, pT, coefMat) + rnorm(1)*sd
    }
  }

  # true bandwidths
  K_PolyMult[s, 1:2] = hPolyFunc(pX, pT, coefMat, nX, nT, sd)

  Smooth = DCSmooth(Y)
  K_PolyMult[s, 3:5] = Smooth$bndw[c(1, 2, 6)]
  K_PolyMult[s, 6:7] = c(nX, nT)
  K_PolyMult[s, 8:9] = c(pX, pT)
  K_PolyMult[s, 10]  = sd
  K_PolyMult[s, 11:13] = Smooth$bndw[3:5]
  K_PolyMult[s, 14:16] = intPolyFunc(pX, pT, coefMat)[3:5]
}

## evaluation

K = data.frame(K_PolyMult[, ])