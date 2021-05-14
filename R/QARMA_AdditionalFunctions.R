################################################################################
#                                                                              #
#                DCSmooth Package: additional Functions for QARMA              #
#                                                                              #
################################################################################

# Simulation of QARMA-processes, estimation of the spectral density


#----------------------------Simulation Function-------------------------------#

#' Simulation of a \eqn{QARMA(p, q)}-process on a lattice.
#' 
#' The MA- and AR-parameters of a top-left quadrant ARMA process are estimated
#' using the Hannen-Rissanen Algorithm (Hannen-Rissanen ???). The lag-orders of 
#' the \eqn{QARMA(p, q)} are given by \eqn{p = (p_1, p_2), q = (q_1, q_2)}{p =
#' (p1, p2), q = (q1, q2)}, where \eqn{p_1, q_1}{p1, q1} are the lags over the
#' rows and \eqn{p_2, q_2}{p2, q2} are the lags over the columns. The estimation
#' process is based on the model
#' \deqn{\phi(B_{1}B_{2})X_{i,j} = \theta(B_{1}B_{2})u_{i,j}}{\phi(B1 B2)
#' X[i,j] = \theta(B1 B2)u[i,j]}
#' 
#' @param Y A numeric matrix that contains the demeaned observations of the
#'   random field or functional time-series.
#' @param model_order A list containing the orders of the QARMA model in the
#'   form \code{model_order = list(ar = c(p1, p2), ma = c(q1, q2))}. Default
#'   value is a \eqn{QARMA((1, 1), (1, 1))} model.
#' 
#' @return The function returns a list including
#' 
#' \describe{
#' \item{ar}{A \eqn{(p_1 + 1) \times (p_2 + 1)}{(p1 + 1) x (p2 + 1)}-matrix
#' containing the estimated AR-parameters \eqn{\phi} of the QARMA-process. The
#' \eqn{[i, j]}th entry is the \eqn{(p_1 + 1 - i, p_2 + 1 - j)}th lag with the
#' \eqn{[p_1 + 1, p_2 + 1]}th entry being 1.}
#' \item{ma}{A \eqn{(q_1 + 1) \times (q_2 + 1)}{(q1 + 1) x (q2 + 1)}-matrix
#' containing the estimated MA-parameters \eqn{\theta} of the QARMA-process. The
#' \eqn{[i, j]}th entry is the \eqn{(q_1 + 1 - i, q_2 + 1 - j)}th lag with the 
#' \eqn{[q_1 + 1, q_2 + 1]}th entry being 1.}
#' \item{sigma}{The estimated standard deviation \eqn{\sigma} of the QARMA-model.}
#' \item{innov}{The matrix of innovations resulting from the QARMA estimation. 
#' The initial values (from \eqn{i = 1,...,p_1}{i = 1,...,p1} and
#' \eqn{j = 1,...,p_2}{j = 1,...,p2}) are set to zero.} 
#' }
#' 
#' @section Usage:
#' \code{qarma(Y, model_order = list(ar = c(1, 1), ma = c(1, 1)))}
#' 
#' @section Details:
#' ???
#' 
#' @examples
#' qarma.example1
#' qarma.example2
#' qarma.est(qarma.example1)
#' qarma.est(qarma.example2, model_order = list(ar = c(1, 1), ma = c(0, 0))
#' 
#' @export

qarma.sim = function(nX, nT, model = list(ar, ma, sigma))
{
  ar_mat = as.matrix(model$ar); ma_mat = as.matrix(model$ma)
  ar_x = dim(ar_mat)[1] - 1; ar_t = dim(ar_mat)[2] - 1
  ma_x = dim(ma_mat)[1] - 1; ma_t = dim(ma_mat)[2] - 1
  xInit = max(ar_x, ma_x) + 1
  tInit = max(ar_t, ma_t) + 1
  
  # set coefficients for zero-lags
  ma_mat[1, 1] = 1 # MA-coefficients
  ar_mat[1, 1] = 0 # AR-coefficients
  
  nMat = floor(1.25 * c(nX, nT))
  errorMat = matrix(rnorm(prod(nMat)), nrow = nMat[1], ncol = nMat[2]) *
    model$sigma
  armaMat = errorMat
  
  for (i in xInit:nMat[1])
  {
    for (j in tInit:nMat[2])
    {
      armaMat[i, j] = sum(-ar_mat * armaMat[i:(i - ar_x), j:(j - ar_t)]) +
        sum(ma_mat * errorMat[i:(i - ma_x), j:(j - ma_t)])
    }
  }
  
  armaOut = armaMat[(nMat[1] - nX + 1):nMat[1], (nMat[2] - nT + 1):nMat[2]]
  errorOut = errorMat[(nMat[1] - nX + 1):nMat[1], (nMat[2] - nT + 1):nMat[2]]
  
  listOut = list(Y = armaOut, innov = errorOut,
                 ar = ar_mat, ma = ma_mat, sigma = model$sigma)
  return(listOut)
}

#----------------------Calculation of spectral density-------------------------#

qarma.SSDE = function(Y, ar, ma, stdev)
{
  X = T = seq(from = -pi, to = pi, length.out = 100)
  ySpectrum = matrix(NA, nrow = 100, ncol = 100)
  
  for (i in 1:100)
  {
    for (j in 1:100)
    {
      omega = c(X[i], T[j])
      ySpectrum[i, j] = qarma.spec(Y, ar, ma, stdev, omega)
    }
  }
  return(ySpectrum)
}

qarma.spec = function(Y, ar, ma, stdev, omega)
{
  lagAR = dim(ar) - 1; lagMA = dim(ma) - 1
  z1 = complex(argument = omega[1]); z2 = complex(argument = omega[2])
  
  g00 = sd(Y)^2
  
  arSum = Re(sum(ar * (z1^(lagAR[1]:0)) %*% t(z2^(lagAR[2]:0))) * 
               sum(ar * ((1/z1)^(lagAR[1]:0)) %*% t((1/z2)^(lagAR[2]:0))))
  maSum = Re(sum(ma * (z1^(lagMA[1]:0)) %*% t(z2^(lagMA[2]:0))) * 
               sum(ma * ((1/z1)^(lagMA[1]:0)) %*% t((1/z2)^(lagMA[2]:0))))
  
  specDensOut = stdev^2/g00 * maSum/arSum
  
  return(specDensOut)
}
