# yw_matrix3 = function(Y, ar_order = c(1, 1), nu_vec = c(1, 1))
# {
#   # set up parameters, result matrices and vectors etc.
#   nX = dim(Y)[1]; nT = dim(Y)[2]; n = nX * nT
#   p1_ar = ar_order[1]
#   p2_ar = ar_order[2]
#   d_ar = (p1_ar + 1) * (p2_ar + 1) - 1
#   
#   acf_matrix_ar = matrix(NA, nrow = d_ar, ncol = d_ar)
#   acf_vector_ar = vector(mode = "numeric", length = d_ar)
#   # st_mat is over rows first (i, j) -> (i - 1, j), which is normal R order
#   st_mat = expand.grid(p1 = c(p1_ar:0), p2 = c(p2_ar:0)) - nu_vec
#   
#   # calculate acf_matrix and acf_vector
#   for (k in 1:d_ar)
#   {
#     for (l in 1:d_ar)
#     {
#       st_vec = st_mat[k, ] - st_mat[l, ]
#       acf_matrix_ar[k, l] = .matrix_acf(Y, st_vec)
#     }
#     acf_vector_ar[k] = .matrix_acf(Y, st_mat[k, ] - nu_vec)
#   }
#   
#   # solving Yule-Walker equations
#   yw.estimators = solve(acf_matrix_ar) %*% acf_vector_ar
#   phi_ar_matrix = matrix(c(yw.estimators, -1), nrow = (p1_ar + 1),
#                          ncol = (p2_ar + 1))
#   
#   return(phi_ar_matrix)
# }
# 
# yw_matrix2 = function(Y, ar_order = c(1, 1), ma_order = c(1, 1))
# {
#   nX0 = dim(Y)[1]; nT0 = dim(Y)[2]
#   Y1 = Y[1:as.numeric(nX0 - ma_order[1]), 1:as.numeric(nT0 - ma_order[2])]
#   Y2 = Y[as.numeric(ma_order[1] + 1):nX0, as.numeric(ma_order[2] + 1):nT0]
#   
#   # set up parameters, result matrices and vectors etc.
#   nX = dim(Y1)[1]; nT = dim(Y1)[2]; n = nX * nT
#   p1_ar = ar_order[1]
#   p2_ar = ar_order[2]
#   d_ar = (p1_ar + 1) * (p2_ar + 1) - 1
#   
#   acf_matrix_ar = matrix(NA, nrow = d_ar, ncol = d_ar)
#   acf_vector_ar = vector(mode = "numeric", length = d_ar)
#   # st_mat is over rows first (i, j) -> (i - 1, j), which is normal R order
#   st_mat = expand.grid(p1 = c(p1_ar:0), p2 = c(p2_ar:0))
#   
#   # calculate acf_matrix and acf_vector
#   for (k in 1:d_ar)
#   {
#     for (l in 1:d_ar)
#     {
#       st_vec = st_mat[k, ] - st_mat[l, ]
#       acf_matrix_ar[k, l] = .matrix_acf2(Y1, Y2, st_vec)
#     }
#     acf_vector_ar[k] = .matrix_acf2(Y1, Y2, st_mat[k, ])
#   }
#   
#   # solving Yule-Walker equations
#   yw.estimators = solve(acf_matrix_ar) %*% acf_vector_ar
#   phi_ar_matrix = matrix(c(yw.estimators, -1), nrow = (p1_ar + 1),
#                          ncol = (p2_ar + 1))
#   
#   return(phi_ar_matrix)
# }
# 
# ### auxiliary functions
# .matrix_acf2 = function(Y1, Y2, st_vec)
# {
#   s = as.numeric(st_vec[1])
#   t = as.numeric(st_vec[2])
#   n = dim(Y1)
#   
#   if ((s >= 0 && t >= 0) || (s < 0 && t < 0))
#   {
#     acf_out = .matrix_acf_positive2(Y1, Y2, abs(s), abs(t), n)
#   } else if (s < 0 && t >= 0) {
#     acf_out = .matrix_acf_negative2(Y1, Y2, s = -s, t, n)
#   } else {
#     acf_out = .matrix_acf_negative2(Y1, Y2, s, t = -t, n)
#   }
#   
#   # unbiased estimator from Ha/Newton(1993)
#   unbiased_factor = (n[1] - s)*(n[2] - t)/prod(n)
#   return(acf_out/unbiased_factor)
# }
# 
# .matrix_acf_positive2 = function(Y1, Y2, s, t, n)
# {
#   Y1a = Y1[1:(n[1] - s), 1:(n[2] - t)]
#   Y2a = Y2[(s + 1):n[1], (t + 1):n[2]]
#   acf_out = sum(Y1a * Y2a)
#   return(acf_out)
# }
# 
# .matrix_acf_negative2 = function(Y1, Y2, s, t, n)
# {
#   Y1a = Y1[1:(n[1] - s), (t + 1):n[2]]
#   Y2a = Y2[(s + 1):n[1], 1:(n[2] - t)]
#   acf_out = sum(Y1a * Y2a)
#   return(acf_out)
# }
# 
# 
# ar_mat = matrix(c(1,-0.5, 0.3, 0.1), 2, 2)
# ma_mat = matrix(c(1, 0.1, 0.5, 0.2), 2, 2)
# 
# nSim = 100
# coef_mat = matrix(NA, nrow = 100, ncol = 4)
# for (i in 1:nSim)
# {
#   Y = qarma.sim(100, 100, model = list(ar = ar_mat, ma = ma_mat, sigma = 0.5))$Y
#   coef_mat[i, ] = yw_matrix2(Y, ma_order = c(1, 1))
# }
# coef_mat
# 
# 
# ar_mat = matrix(c(1, 0.3, 0.3, 0.3), 2, 2)
# ma_mat = matrix(c(1, 0.3, 0.3, 0.3), 2, 2)
# ma_vec = expand.grid(0:2, 0:2)
# ma_vec = ma_vec[order(rowSums(ma_vec), apply(ma_vec, 1, max)), ]
# ar_vec = ma_vec[2:dim(ma_vec)[1], ]
# test_mat = matrix(NA, nrow = dim(ma_vec)[1], ncol = dim(ar_vec)[1])
# colnames(test_mat) = paste0("(", ar_vec[, 1], ", ", ar_vec[, 2], ")")
# rownames(test_mat) = paste0("(", ma_vec[, 1], ", ", ma_vec[, 2], ")")
# 
# test_arr = array(dim = c(9, 8, 100))
# 
# for(u in 1:100)
# {
# print(u)
# for (k in 1:dim(ar_vec)[1])
# {
#   for (j in 1:dim(ma_vec)[1])
#   {
#     nX0 = dim(Y)[1]; nT0 = dim(Y)[2]
#     ar = ar_vec[k, ]
#     ma = ma_vec[j, ]
#     Y = qarma.sim(100, 100, model = list(ar = ar_mat, ma = ma_mat, sigma = 0.1))$Y
#     test_mat[j, k] = yw_matrix3(Y, ar_order = as.numeric(ar),
#                                         nu_vec = as.numeric(ma))[1, 1]
#   }
# }
# test_arr[, , u] = test_mat
# }
# 
# colnames(test_arr) = paste0("(", ar_vec[, 1], ", ", ar_vec[, 2], ")")
# rownames(test_arr) = paste0("(", ma_vec[, 1], ", ", ma_vec[, 2], ")")
# apply(test_arr, c(1, 2), mean)
# apply(test_arr, c(1, 2), sd)
