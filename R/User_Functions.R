###############################################################################
#                                                                             #
#                     DCSmooth Package: User Functions                        #
#                                                                             #
###############################################################################


#-----------------------Set Options via Function------------------------------#

#' Set Options for the DCS procedure
#' 
#' @param kernPar kernel parameters \eqn{k} and \eqn{\mu} as vector.
#' @param pOrder order of polynomials used for smoothing with default value 1 
#' (local linear regression). Kernel regression is done via \code{pOrder = 0}.
#' Orders are assumed to be the same in both directions.
#' @param inflExp inflation exponents for bandwidth selection.
#' @param inflPar inflation parameters for bandwidth selection.
#' @param delta shrink parameters for derivatives.
#' @param constWindow use a constant window width for estimation.
#' @param varEst procedure for estimation of variance factor \eqn{c_f}.
#' 
#' @export
#' 
#' @examples
#' myOpt = setOptions(kernPar = c(4,2), pOrder = 0, inflExp = c(0.5, 0.5),
#' inflPar = c(2, 1), delta = c(0.05, 0.05), constWindow = TRUE,
#' varEst = "iid")

setOptions = function(...) # user friendly wrapper function for .setOptions
{
  .setOptions(...)
}

#--------------Function for smoothing and bandwidth estimation----------------#

#' Nonparametric Double Conditional Smoothing for 2D Surfaces
#' 
#' \code{DCSmooth} provides a double conditional nonparametric smoothing of the
#' expectation surface of a functional time series or a random field on a
#' lattice. Bandwidth selection is done via an iterative plug-in method.
#' 
#' @section Usage
#' \code{DCSmooth(Y, X = 1, T = 1, bndw = "auto", dcsOptions = setOptions())}
#' 
#' @section Details
#' Blafasel
#' 
#' @param Y A numeric matrix that contains the observations of the random field
#'   or functional time-series.
#' @param X An optional numeric vector containing the exogenous covariates
#'   with respect to the rows.
#' @param T An optional numeric vector containing the exogenous covariates
#'   with respect to the columns.
#' @param bndw Bandwidth for smoothing the observations in \code{Y}. Can be a
#'   two-valued numerical vector with bandwidths in row- and column-direction.
#'   If the value is \code{"auto"} bandwidth selection will be carried out by
#'   the iterative plug-in algorithm.
#' @param dcsOptions Additional options for the smoothing procedure. Is set
#'   via the \code{setOptions} command.
#' 
#' @return \code{DCSmooth} returns an object of class "dcs".
#' 
#' @section Details
#' The function \code{summary}
#' 
#' @example 

DCSmooth = function(Y, X = 1, T = 1, bndw = "auto", dcsOptions = setOptions())
{
  # set up vectors for X and T
  if (length(X) == 1 && length(T) == 1)
  {
    X = seq(from = 0, to = 1, length.out = dim(Y)[1])
    T = seq(from = 0, to = 1, length.out = dim(Y)[2])
  }
  
  # set kernel function (type, k, mu, nu = 0)
  nameKernFcn = paste0("MW", dcsOptions$kernPar[1], dcsOptions$kernPar[2], "0")
  kernelFcn  = kernelFcn_assign(nameKernFcn) # set kernel Function to use in optimization

  # bandwidth selection process
  if (bndw == "auto" && dcsOptions$pOrder == 0)
  {
    # kernel regression
    bndwObj = KR_bndwSelect(Y, kernelFcn, dcsOptions)
    bndw = pmin(bndwObj$bndw, c(0.45, 0.45)) # KR cannot handle larger bndws.
  } else if (bndw == "auto" && dcsOptions$pOrder > 0) {
    # local polynomial regression
    bndwObj = LP_bndwSelect(Y, kernelFcn, dcsOptions)
    bndw = bndwObj$bndw
  }
  
  # estimation of resulting surface
  if (dcsOptions$pOrder == 0) # kernel regression
  {
    DCSOut = KR_DoubleSmooth2(yMat = Y, hVec = bndw, drvVec = c(0, 0), 
      kernFcnPtrX = kernelFcn, kernFcnPtrT = kernelFcn)
  } else if (dcsOptions$pOrder > 0) {
    # local polynomial regression
    DCSOut = LP_DoubleSmooth2(yMat = Y, hVec = bndw, polyOrderVec 
      = c(dcsOptions$pOrder, dcsOptions$pOrder), drvVec = c(0,0),
      kernFcnPtr = kernelFcn)
  }
  
  DCS_out = list(X = X, T = T, Y = Y, M = DCSOut, bndw = bndw,
                 varCoef = bndwObj$varCoef, iterations = bndwObj$iterations,
                 dcsOptions = dcsOptions)
  
  # apply class to output object (not finished)
  new("DCSmooth", DCS_out)
  
  #return(DCS_out)
}

#-------------------Function for plotting smoothed surface--------------------#

plotDCS = function(DCSobj, ...)
{
  .plotDCS(DCSobj = DCSobj, ...)
}
