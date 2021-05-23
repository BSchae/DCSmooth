################################################################################
#                                                                              #
#                DCSmooth Package: additional Functions for QARMA              #
#                                                                              #
################################################################################

# Simulation of QARMA-processes, estimation of the spectral density


#----------------------------Simulation Function-------------------------------#

#' Simulation of a \eqn{QARMA(p, q)}-process
#' 
#' @description \code{qarma.sim} is used to simulate a specified QARMA-model
#'  on a lattice, with normally distributed innovations.
#' 
#' @section Details:
#' The values of a \eqn{QARMA((p_2, p_2), (q_1, q_2))}-model are simulated from
#'  given matrices \eqn{\phi} (\code{ar}) and \eqn{\theta} (\code{ma}) and 
#'  a standard deviation for the innovations \eqn{\sigma} (\code{sigma}). The
#'  simulation process is based on the model
#'  \deqn{\phi(B_{1}B_{2})X_{i,j} = \theta(B_{1}B_{2})u_{i,j}}{\phi(B1 B2)
#'  X[i,j] = \theta(B1 B2)u[i,j]}
#'  where the innovations \eqn{u} are iid. draws from a standard normal distribution
#'  with zero mean and variance \eqn{\sigma^2}
#' 
#' @param Y A numeric matrix that contains the demeaned observations of the
#'   random field or functional time-series.
#' @param model_order A list containing the orders of the QARMA model in the
#'   form \code{model_order = list(ar = c(p1, p2), ma = c(q1, q2))}. Default
#'   value is a \eqn{QARMA((1, 1), (1, 1))} model.
#' 
#' @return The function returns a list, consisting of
#'  \tabular{ll}{
#'   \code{ar} \tab a \eqn{(p_1 + 1) \times (p_2 + 1)}{(p1 + 1) x (p2 + 1)}-matrix
#'     containing the given AR-parameters \eqn{\phi} of the QARMA-process. \cr
#'   \code{ma} \tab a \eqn{(q_1 + 1) \times (q_2 + 1)}{(q1 + 1) x (q2 + 1)}-matrix
#'     containing the given MA-parameters \eqn{\theta} of the QARMA-process. \cr
#'   \code{sigma} \tab the standard deviation \eqn{\sigma} of the innovations. \cr
#'   \code{innov} \tab the simulated matrix of innovations. \cr
#' }
#' 
#' @section Details: see the vignette for further details.
#' 
#' @seealso \code{\link{qarma.est}}
#' 
#' @examples
#' ma = matrix(c(1, 0.2, 0.4, 0.1), nrow = 2, ncol = 2)
#' ar = matrix(c(1, 0.5, -0.1, 0.1), nrow = 2, ncol = 2)
#' sigma = 0.5
#' q_model = list(ar = ar, ma = ma, sigma = sigma)
#' 
#' q_sim = qarma.sim(100, 100, model = model)
#' q_sim
#' 
#' @export

qarma.sim = function(n_x, n_t, model = list(ar, ma, sigma))
{
  ar_mat = as.matrix(model$ar); ma_mat = as.matrix(model$ma)
  ar_x = dim(ar_mat)[1] - 1; ar_t = dim(ar_mat)[2] - 1
  ma_x = dim(ma_mat)[1] - 1; ma_t = dim(ma_mat)[2] - 1
  x_init = max(ar_x, ma_x) + 1
  t_init = max(ar_t, ma_t) + 1
  
  # set coefficients for zero-lags
  ma_mat[1, 1] = 1 # MA-coefficients
  ar_mat[1, 1] = 0 # AR-coefficients
  
  n_mat = floor(1.25 * c(n_x, n_t))
  error_mat = matrix(rnorm(prod(n_mat)), nrow = n_mat[1], ncol = n_mat[2]) *
    model$sigma
  arma_mat = error_mat
  
  for (i in x_init:n_mat[1])
  {
    for (j in t_init:n_mat[2])
    {
      arma_mat[i, j] = sum(-ar_mat * arma_mat[i:(i - ar_x), j:(j - ar_t)]) +
        sum(ma_mat * error_mat[i:(i - ma_x), j:(j - ma_t)])
    }
  }
  
  arma_out = arma_mat[(n_mat[1] - n_x + 1):
                        n_mat[1], (n_mat[2] - n_t + 1):n_mat[2]]
  error_out = error_mat[(n_mat[1] - n_x + 1):
                          n_mat[1], (n_mat[2] - n_t + 1):n_mat[2]]
  
  list_out = list(Y = arma_out, innov = error_out,
                 ar = ar_mat, ma = ma_mat, sigma = model$sigma)
  return(list_out)
}

#----------------------Calculation of spectral density-------------------------#

qarma.ssde = function(Y, ar, ma, st_dev)
{
  X = T = seq(from = -pi, to = pi, length.out = 100)
  y_spectrum = matrix(NA, nrow = 100, ncol = 100)
  
  for (i in 1:100)
  {
    for (j in 1:100)
    {
      omega = c(X[i], T[j])
      y_spectrum[i, j] = qarma.spec(Y, ar, ma, st_dev, omega)
    }
  }
  return(y_spectrum)
}

qarma.spectral.density = function(Y, ar, ma, st_dev, omega)
{
  lag_ar = dim(ar) - 1; lag_ma = dim(ma) - 1
  z1 = complex(argument = omega[1]); z2 = complex(argument = omega[2])
  
  g00 = sd(Y)^2
  
  ar_sum = Re(sum(ar * (z1^(lag_ar[1]:0)) %*% t(z2^(lag_ar[2]:0))) * 
               sum(ar * ((1/z1)^(lag_ar[1]:0)) %*% t((1/z2)^(lag_ar[2]:0))))
  ma_sum = Re(sum(ma * (z1^(lag_ma[1]:0)) %*% t(z2^(lag_ma[2]:0))) * 
               sum(ma * ((1/z1)^(lag_ma[1]:0)) %*% t((1/z2)^(lag_ma[2]:0))))
  
  spectral_dens_out = st_dev^2/g00 * ma_sum/ar_sum
  
  return(spectral_dens_out)
}
