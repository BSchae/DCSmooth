library(KernSmooth)
library(DCS)
library(mvtnorm)
library(rgl)

setwd("C:/Daten/Forschung/Bandwidth Selection for the DCS/Data & Pictures/Neue Daten")

### Function for Color of 3D-Plots
col_function = function(Y, n_col = 100) {
  Y_col = Y - min(Y)
  Y_col = Y_col/max(Y_col)*(n_col - 1) + 1
  #hcl.pals()
  colorlut = hcl.colors(n_col, palette = "Plasma", rev = TRUE) # height color lookup table
  col = colorlut[Y_col] # assign colors to heights for each point
  dim(col) = dim(Y_col)
  return(col)
}

### Functions for Calculation of MISE
MISE = function(h, Y, Y0, X, T, delta) {
  Smooth = DCS(Y, X, T, smoptions = list(bndw = h))$est
  range_x = (n_x*delta):(n_x*(1 - delta))
  range_t = (n_t*delta):(n_t*(1 - delta))
  MISE = mean((Smooth[range_x, range_t] - Y0[range_x, range_t])^2)
  return(MISE)
}

### General settings
n_sim = 500
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
XT_combine = expand.grid(X, T)

f_sinus = function(X) {
  sin(X[, 1]*3*pi) * sin(2*X[, 2]*pi) * exp(-2*X[, 1]*X[, 2])
}

Y0 = f_sinus(XT_combine)
Y0 = matrix(Y0, nrow = n_x, ncol = n_t)
Y  = Y0 + matrix(sqrt(sigma_sq)*rnorm(n_x*n_t), nrow = n_x, ncol = n_t)

### Smoothing
K_r2  = matrix(NA, nrow = n_sim, ncol = 4)
set.seed(seed)

for(i in 1:n_sim) {
  Y = Y0 + matrix(sqrt(sigma_sq)*rnorm(n_x*n_t), nrow = n_x, ncol = n_t)
  if(i > 485) {
    Smooth = DCS(Y, X, T, smoptions = list(type = "K", shrink_par = 0))
    K_r2[i, 1:2] = Smooth$bndw
    K_r2[i, 3] = Smooth$iterations
    K_r2[i, 4] = mean((Smooth$est - Y0)^2)
  }
  print(i)
}

# Save Values
colnames(K_r2) = c("h_x", "h_t", "Iterations", "MISE")
write.csv(K_r2, "K_r2.csv", row.names = FALSE)

### Calculation of MISE optimal bandwidths
n_h = 20
h_values = seq(from = 0.04, to = 0.06, length.out = n_h)
h_test = expand.grid(h_values, h_values)
MISE_r2 = matrix(NA, nrow = n_sim, ncol = 3)

set.seed(seed)
for(i in 1:n_sim) {
  Y = Y0 + matrix(sqrt(sigma_sq)*rnorm(n_x*n_t), nrow = n_x, ncol = n_t)
  #MISE_temp = (1:n_h)*NA
  #for(j in 1:n_h^2) {
  #    MISE_temp[j] = MISE(c(h_test[j, 1], h_test[j, 2]), Y, Y0, X, T, 0)
  #}
  #place = which(MISE_temp == min(MISE_temp, na.rm = TRUE))[1]
  optim(c())
  MISE_r2[i, 1:2] = c(h_test[place, 1], h_test[place, 2])
  MISE_r2[i, 3] = MISE_temp[place]
  print(i)
}

# Save values
colnames(MISE_r2) = c("h_x", "h_t", "MISE")
write.csv(MISE_r2, "MISE_r2_500.csv", row.names = FALSE)

### Persp Plot
DCS_K  = DCS(Y, X, T, smoptions = list(type = "K"))

persp3d(X, T, Y0, color = col_function(Y0), alpha = 1, zlab = "")

user_mat = matrix(c(0.7725729, 0.6349250, -0.001069574, 0, -0.2014114,
  0.2466731, 0.947937787, 0, 0.6021332, -0.7321358, 0.318453997, 0,
  0, 0, 0, 1), 4, 4, byrow = TRUE)
par3d(windowRect= 50 + c(0, 0, 640, 640))
view3d(userMatrix = user_mat, fov = 0)


rgl.snapshot("Sim_r2_Y0.png", fmt = "png", top = TRUE)
persp3d(X, T, Y, color = col_function(Y), alpha = 1, zlab = "")
rgl.snapshot("Sim_r2_Y.png", fmt = "png", top = TRUE)
persp3d(X, T, DCS_K$est, color = col_function(DCS_K$est), alpha = 1, zlab = "")
rgl.snapshot("Sim_r2_DCS.png", fmt = "png", top = TRUE)

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
  Smooth = DCS(Y, X, T, smoptions = list(type = "K", shrink_par = 0))
  K_r1[i, 1:2] = Smooth$bndw
  K_r1[i, 3] = Smooth$iterations
  K_r1[i, 4] = mean((Smooth$est - Y0)^2)
  print(i)
}

# Save Values
colnames(K_r1) = c("h_x", "h_t", "Iterations", "MISE")
write.csv(K_r1, "K_r1.csv", row.names = FALSE)

# quantiles of evaluation
quantile(K_r1[, 1], probs = c(0.1, 0.5, 0.9))
quantile(K_r1[, 2], probs = c(0.1, 0.5, 0.9))

### Calculation of MISE
n_h = 21
h_values = seq(from = 0.1, to = 0.13, length.out = n_h)
h_test = expand.grid(h_values, h_values)
MISE_r1 = matrix(NA, nrow = n_sim, ncol = 3)

set.seed(seed)
for(i in 1:n_sim) {
  Y = Y0 + matrix(sqrt(sigma_sq)*rnorm(n_x*n_t), nrow = n_x, ncol = n_t)
  MISE_temp = optim(c(0.12, 0.12), MISE, Y=Y, Y0=Y0, X=X, T=T, delta=0)
  MISE_r1[i, 1:2] = MISE_temp$par
  MISE_r1[i, 3] = MISE_temp$value
  print(i)
}

# Save values
colnames(MISE_r1) = c("h_x", "h_t", "MISE")
write.csv(MISE_r1, "MISE_r1_400.csv", row.names = FALSE)

### Persp Plot
DCS_K  = DCS(Y, X, T, smoptions = list(type = "K"))

persp3d(X, T, Y0, color = col_function(Y0), alpha = 1, zlab = "")
rgl.snapshot("Sim_r1_Y0.png", fmt = "png", top = TRUE)
persp3d(X, T, Y, color = col_function(Y), alpha = 1, zlab = "")
rgl.snapshot("Sim_r1_Y.png", fmt = "png", top = TRUE)
persp3d(X, T, DCS_K$est, color = col_function(DCS_K$est), alpha = 1, zlab = "")
rgl.snapshot("Sim_r1_DCS.png", fmt = "png", top = TRUE)

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
Y = Y0 + sqrt(sigma_sq) * matrix(rnorm(n_x*n_t), nrow = n_x, ncol = n_t)

### Smoothing
K_saddle = matrix(NA, nrow = n_sim, ncol = 4)
set.seed(seed)

for (i in 1:n_sim) {
  Y = Y0 + sqrt(sigma_sq) * matrix(rnorm(n_x*n_t), nrow = n_x, ncol = n_t)
  if(i > 467) {
    DCS_out = DCS(Y, X, T, smoptions = list(type = "K", shrink_par = 0))
    K_saddle[i, 1:2] = DCS_out$bndw
    K_saddle[i, 3] = DCS_out$iterations
    K_saddle[i, 4] = mean((DCS_out$est - Y0)^2)
  }
  print(i)
}

# Save values
colnames(K_saddle) = c("h_x", "h_t", "Iterations", "MISE")
write.csv(K_saddle, "K_saddle.csv", row.names = FALSE)

### Persp Plot
DCS_K  = DCS(Y, X, T, smoptions = list(type = "K"))

persp3d(X, T, Y0, color = col_function(Y0), alpha = 1, zlab = "")
rgl.snapshot("Sim_saddle_Y0.png", fmt = "png", top = TRUE)
persp3d(X, T, Y, color = col_function(Y), alpha = 1, zlab = "")
rgl.snapshot("Sim_saddle_Y.png", fmt = "png", top = TRUE)
persp3d(X, T, DCS_K$est, color = col_function(DCS_K$est), alpha = 1, zlab = "")
rgl.snapshot("Sim_saddle_DCS.png", fmt = "png", top = TRUE)

#-------------------------------------------------------------------------------#
#                               Evaluations                                     #
#-------------------------------------------------------------------------------#

par(mfrow = c(3, 2), mar = c(4, 6, 2.5, 2.5), oma = c(5, 0, 2, 0), xpd = FALSE)

### Single-Peak Function
# Read Data
K_r1 = read.csv("K_r1_complete.csv")

# Optimal Bandwidths
h_opt = c(0.11961757, 0.11961757) # r1, n_X = 400, n_T = 250

# Plot
hist(K_r1[, 1], freq = FALSE, xlab = expression(h[x]), ylab = "", main = "", ylim = c(0, 160))
lines(density(K_r1[, 1]), col = 2, lwd = 2)
abline(v = h_opt[1], col =1, lwd = 2, lty = 2)

hist(K_r1[, 2], freq = FALSE, xlab = expression(h[t]), ylab = "", main = "")
lines(density(K_r1[, 2]), col = 2, lwd = 2)
abline(v = h_opt[2], col = 1, lwd = 2, lty = 2)

### Sine Function
# Read Data
K_r2 = read.csv("K_r2_complete.csv")

# Optimal Bandwidths
h_opt = c(0.05096990, 0.07512002) # r2, n_X = n_T = 500

# Plot
hist(K_r2[, 1], freq = FALSE, xlab = expression(h[x]), ylab = "", main = "", ylim = c(0, 1000))
lines(density(K_r2[, 1]), col = 2, lwd = 2)
abline(v = h_opt[1], col =1, lwd = 2, lty = 2)

hist(K_r2[, 2], freq = FALSE, xlab = expression(h[t]), ylab = "", main = "", ylim = c(0, 600))
lines(density(K_r2[, 2]), col = 2, lwd = 2)
abline(v = h_opt[2], col = 1, lwd = 2, lty = 2)

### Saddle Function
# Read Data
K_saddle = read.csv("K_saddle_complete.csv")

# Optimal Bandwidths
h_opt = c(0.065228949, 0.092247665) # saddle, n_X = 1500, n_T = 510

# Plot
hist(K_saddle[, 1], freq = FALSE, xlab = expression(paste("log(", h[x], ")")), ylab = "", main = "", ylim = c(0, 650))
lines(density(K_saddle[, 1]), col = 2, lwd = 2)
abline(v = h_opt[1], col =1, lwd = 2, lty = 2)

hist(K_saddle[, 2], freq = FALSE, ylim = c(0, 500), xlab = expression(h[t]), ylab = "", main = "")
lines(density(K_saddle[, 2]), col = 2, lwd = 2)
abline(v = h_opt[2], col =1, lwd = 2, lty = 2)

### Title and Legend
mtext(expression(paste("Gaussian-Peak Function ", f[1])), side = 3, line = -1, outer = TRUE)
mtext(expression(paste("Sine Function ", f[2])), side = 3, line = -28, outer = TRUE)
mtext(expression(paste("Saddle-Type Function ", f[3])), side = 3, line = -56, outer = TRUE)
par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
legend("bottom", legend = c("IPI Bandwidhts", "Asymptotic Optimal"),
  horiz = TRUE, col = c("red", "black"), lty = c(1, 2), lwd = 2,
  cex = 1.5, bty = "n")
