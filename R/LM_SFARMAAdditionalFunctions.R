################################################################################
#                                                                              #
#                DCSmooth Package: additional Functions for SFARIMA              #
#                                                                              #
################################################################################

# Simulation of SFARIMA-processes

sfarima.sim <- function(n_x, n_t, model, nstart = 150, K1 = 150, K2 = 150) {
  n_x = n_x + nstart
  n_t = n_t + nstart
  eps_mat = matrix(rnorm(n_x * n_t), n_x, n_t) * model$sigma
  a1 = macoef(ar = model$ar1, ma = model$ma1, d = model$d1, K1)
  a2 = macoef(ar = model$ar2, ma = model$ma2, d = model$d2, K2)
  
  X1.sim = X2.sim = matrix(0, n_x, n_t)

    a2 = t(a2)
    eps_mat = t(eps_mat)

    for(j in 1:n_t) {
      if (j <= K2) {
        X2.sim[, j] = a2[1:j] %*% eps_mat[j:1, ]
      }
      else {
        X2.sim[, j] = a2 %*% eps_mat[j:(j - K2), ]
      }
    }

    a1 = t(a1)

    for(i in 1:n_x) {
      if (i <= K1) {
        X1.sim[i, ] = a1[1:i] %*% X2.sim[i:1, ]
      }
      else {
        X1.sim[i, ] = a1 %*% X2.sim[i:(i - K1), ]
      }

    }
    return(X1.sim[(nstart + 1):n_x, (nstart + 1):n_t])
}

#----------------------------------------------------------------#

macoef <- function(ar = 0, ma = 0, d = 0, k = 50) {
  p = length(ar[ar != 0])
  q = length(ma[ma != 0])
  if (p == 0) {
    ar = 0
  }
  if (q == 0) {
    ma = 0
  }
  ma.coef = c(1, ma, rep(0, k - q))
  arma.coef = (1:(k + 1)) * 0
  arma.coef[1] = 1
  
  if (p | q > 0) {
    for (i in 2:(k + 1)) {
      if ((i - p) < 1) {
        arma.coef[i] = sum(ar[1:(p - abs(i - p) - 1)] * arma.coef[(i - 1):1]) - ma.coef[i]
      } else {
        arma.coef[i] = sum(ar[1:p] * arma.coef[(i - 1):(i - p)]) - ma.coef[i]
      }
    }
  } else {
    arma.coef + 1
  }
  
  d.coef = choose(d, 0:k) * ((-1)^(0:k))
  
  coef.all = (1:(k + 1)) * 0
  for (j in 1:(k + 1)) {
    coef.all[j] = sum(d.coef[1:j] * arma.coef[j:1])
  }
  return(coef.all)
}



