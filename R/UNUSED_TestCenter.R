library(mvtnorm)
library(rgl)

### General settings
n_sim = 100
seed = 42

#-------------------------------------------------------------------------------#
#                           Sine-Function                                       #
#-------------------------------------------------------------------------------#

### Set up Function
n_x = 500
n_t = 500
sigma_sq = 0.5

X = seq(from = 0, to = 1, length.out = n_x)
T = seq(from = 0, to = 1, length.out = n_t)

Y0 = matrix(NA, nrow = n_x, ncol = n_t)

for (i in 1:n_x) {
  for (j in 1:n_t) {
    Y0[i, j] =   sin(3*pi*X[i]) * sin(2*pi*T[j]) * exp(-2*X[i]*T[j])
  }
}

Y  = Y0 + matrix(sqrt(sigma_sq)*rnorm(n_x*n_t), nrow = n_x, ncol = n_t)

### Smoothing
K_sine  = matrix(NA, nrow = n_sim, ncol = 3)
set.seed(seed)

t0 = Sys.time()
for(i in 1:n_sim) {
  Y = Y0 + matrix(sqrt(sigma_sq)*rnorm(n_x*n_t), nrow = n_x, ncol = n_t)
  Smooth = DCSmooth(Y, dcsOptions = setOptions(inflPar = c(1.2, 0.5)))
  K_sine[i, 1:2] = Smooth$bndw
  K_sine[i, 3]   = mean((Smooth$M - Y0)^2)
  print(i)
}
t1 = Sys.time()
t1 - t0

plotDCS(Smooth)

#-------------------------------------------------------------------------------#
#         Normal Distribution with single Peak (Herrmann et al. 1995)           #
#-------------------------------------------------------------------------------#

### Set up Function
n_x = 400
n_t = 250
sigma_sq = 0.1

m1 = c(0.5, 0.5)
s1 = 0.2*diag(2)

X = seq(from = 0, to = 1, length.out = n_x)
T = seq(from = 0, to = 1, length.out = n_t)
XT_combine = expand.grid(X, T)

Y0 = dmvnorm(XT_combine, m1, s1)
Y0 = matrix(Y0, nrow = n_x, ncol = n_t)
Y  = Y0 + matrix(sqrt(sigma_sq)*rnorm(n_x*n_t), nrow = n_x, ncol = n_t)

### Smoothing
K_r1  = matrix(NA, nrow = n_sim, ncol = 4)
set.seed(seed)

for(i in 1:n_sim) {
  Y = Y0 + matrix(sqrt(sigma_sq)*rnorm(n_x*n_t), nrow = n_x, ncol = n_t)
  Smooth = DCSmooth(Y, dcsOptions = setOptions(inflPar = c(1.2, 0.5)))
    K_r1[i, 1:2] = Smooth$bndw
    K_r1[i, 3] = mean((Smooth$est - Y0)^2)
  print(i)
}

#-------------------------------------------------------------------------------#
#                           Saddle Function                                     #
#-------------------------------------------------------------------------------#

### Set up Function
Y.fkt = function(x, t) { 0.5*(t^3 - t) - (x^3 - x) }

sigma_sq = 0.01
n_x = 1500
n_t = 510

X = 0:(n_x - 1)/(n_x - 1)
T = 0:(n_t - 1)/(n_t - 1)

Y0 = matrix(NA, nrow = n_x, ncol = n_t)

for (i in 1:n_x) {
  for (j in 1:n_t) {
    Y0[i, j] = Y.fkt(X[i], T[j])
  }
}

### Smoothing
K_saddle = matrix(NA, nrow = n_sim, ncol = 4)
set.seed(seed)

for (i in 1:n_sim) {
  Y = Y0 + sqrt(sigma_sq) * matrix(rnorm(n_x*n_t), nrow = n_x, ncol = n_t)
  DCS_out = DCSmooth(Y)
  K_saddle[i, 1:2] = DCS_out$bndw
  K_saddle[i, 3] = mean((DCS_out$est - Y0)^2)
  print(i)
}

