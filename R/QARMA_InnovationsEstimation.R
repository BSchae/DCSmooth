# X = arima.sim(model = list(ar = 0.5, ma = 0.4), n = 100)
# 
# 
# phi = 0.5
# theta = 0.4
# 
# logL = function(X, phi, theta)
# {
#   n = length(X)
#   theta_vec = c(1, theta)
#   
#   p = length(phi)
#   q = length(theta)
#   m = max(p, q)
#   
#   psi = rep(0, times = n + p)
#   psi[p + 1] = 1
#   
#   for (j in 1:100)
#   {
#     psi[j + p + 1] = ifelse(j <=  q, theta[j], 0) + sum(phi * psi[(j + p):(j + 1)])
#   }
#   
#   psi = psi[(p + 1):(n + p)]
#   
#   gamma = 0:(n - 1)
#   for (i in 1:n)
#   {
#     gamma[i] = sum(psi[1:(n - i + 1)]*psi[i:n])
#   }
#   
#   kappa = matrix(0, nrow = n + 1, n + 1)
#   for (i in 1:m)
#   {
#     for (j in 1:i)
#     {
#       kappa[i, j] = gamma[i - j + 1]
#     }
#   }
#   
#   for (i in (m + 1):(2*m))
#   {
#     for (j in 1:m)
#     {
#       kappa[i, j] = gamma[i - j] - sum(phi * gamma[abs(1 - i + j)])
#     }
#   }
#   
#   for (i in (m + 1):(n + 1))
#   {
#     for (j in 1:i)
#     {
#       kappa[i, j] = ifelse((i - j) <= q, sum(theta_vec[1:(q - (i - j) + 1)] * 
#                           theta_vec[(i - j + 1):(q + 1)]), 0)
#     } 
#   }
# 
#   v = 1:(n + 1)
#   theta_mat = kappa*0
#   v[1] = kappa[1, 1]
#   
#   for (i in 1:n)
#   {
#     for (k in 0:(i - 1))
#     {
#       theta_mat[i, i - k] = 1/v[k + 1] * (
#         kappa[i + 1, k + 1] - 
#           ifelse(k > 0, sum(theta_mat[k, 1:k] * theta_mat[i, i:(i - k + 1)] * v[1:k]), 0)
#       )
#       v[i + 1] = kappa[i + 1, i + 1] - sum(theta_mat[i, i:1]^2 * v[1:i])
#     }
#   }
#   
#   # calculate predictors of X[n + 1]
#   X_pred = X*0
#   if (m > 1)
#   {
#     for (i in 1:(m - 1))
#     {
#       X_pred[i + 1] = sum(theta_mat[i, 1:i] * (X[i:1] - X_pred[i:1]))
#     }
#   }
#   
#   for (i in m:(n - 1))
#   {
#     X_pred[i + 1] = phi * X[i:(i - p + 1)] + 
#               sum(theta_mat[i, 1:q] * (X[i:(i - q + 1)] - X_pred[i:(i - q + 1)]))
#   }
# 
#   # calculate log-likelihood
#   S = sum((X - X_pred)^2/v[1:n])
#   logL = log(1/n * S) + 1/n * sum(log(v[1:n]))
# 
#   return(logL)
# }