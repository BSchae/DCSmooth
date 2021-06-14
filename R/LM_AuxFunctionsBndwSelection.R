################################################################################
#                                                                              #
#       DCSmooth Package: Auxiliary Functions for Long Memory Bandwidths       #
#                                                                              #
################################################################################

#----------------------Formula for optimal bandwidths--------------------------#

h.opt.LM = function(mxx, mtt, d_vec, var_coef, p_order, drv_vec, n_1, n_2,
                    kern_fcn_x, kern_fcn_t)
{
  # calculation of integrals
  i11 = sum(mxx^2)/n_sub; i22 = sum(mtt^2)/n_sub; i12 = sum(mxx * mtt)/n_sub
  
  # kernel constants (kernel Functions may also depend on p, drv)
  kernel_prop_1 = kernel.prop.LP(kern_fcn_x, p_order[1], drv_vec[1])
  kernel_prop_2 = kernel.prop.LP(kern_fcn_t, p_order[2], drv_vec[2])
  
  # compute additional values
  cb = h.opt.LM.cb(i11, i22, i12, kernel_prop_1, kernel_prop_2, drv_vec, d_vec)
  cn = n2/n1
  C1 = 4 * prod(var_coef) * kernel_prop_1$R * kernel_prop_2$R
  delta = p_order[1] + 1 - drv_vec[1]
  # var_coef should contain c_f1, c_f2, and sigma^2 as a vector
  # replace "R" by "V" later.
  
  C1A = (1 - 2*d_vec[1] + 2*drv_vec[1]) * C1 * cn^(diff(d_vec)) /
        (2 * delta * kernel_prop_1$mu * cb^(1 - 2*d_vec[2] + 2*drv_vec[2]) *
        (kernel_prop_1$mu * i11 + kernel_prop_2$mu * i12 * cb^delta))
  C2A = (1 - 2*d_vec[2] + 2*drv_vec[2]) * C1 * cn^(diff(d_vec)) /
        (2 * delta * kernel_prop_2$mu * cb^(1 - 2*d_vec[2] + 2*drv_vec[2]) *
        (kernel_prop_1$mu * i12 * cb^(-delta) + kernel_prop_2$mu * i22))
  
  b1A = (C1A / n_sub^(1 - sum(d_vec)))^(1/
                                (2*(1 - sum(d_vec) + delta + sum(drv_vec))))
  b2A = (C2A / n_sub^(1 - sum(d_vec)))^(1/
                                (2*(1 - sum(d_vec) + delta + sum(drv_vec))))
    
  return(c(b1A, b2A))
}

h.opt.LM.cb = function(i11, i22, i12, kernel_prop_1, kernel_prop_2, drv_vec,
                       d_vec)
{
  denom_value = (1 - 2*d_vec[1] + 2*drv_vec[1]) * kernel_prop_2$nu^2 * i22
  sec_term = (diff(d_vec) - diff(drv_vec)) * kernel_prop_1$mu *
              kernel_prop_2$mu * i12
  sqrt_value = sec_term^2 + denom_value *
                (1 - 2*d_vec[2] + 2*drv_vec[2]) * kernel_prop_1$mu^2 * i11
  
  delta = p_order[1] + 1 - drv_order[1]
  cb_return = 1/denom_value * (sqrt(sqrt_value) - sec_term)
  
  return(cb_return^(1/delta))
}

kernel.prop.LM = function(kernel_fcn, p, drv, d, n_int = 100)
{
  u_seq  = seq(from = -1, to = 1, length.out = n_int)
  u_grid = expand.grid(u_seq, u_seq)
  
  val_V  = gamma(1 - 2*d) * sin(pi * d) *
           sum(kernel_fcn(u_grid[1]) * kernel_fcn(u_grid[2]) *
           (abs(u_grid[, 1] - u_grid[, 2]) + diff(u_seq)[1])^(2*d - 1)) / n_int^2
  val_mu = sum((add_weights * kernel_fcn_use(u_seq, q = 1, kernel_fcn)) *
                 u_seq^(p + 1)) / (n_int * factorial(p + 1))
  
  return(list(R = val_R, mu = val_mu))
}