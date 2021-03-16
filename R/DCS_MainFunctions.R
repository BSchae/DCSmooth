###############################################################################
#                                                                             #
#                     DCSmooth Package: User Functions                        #
#                                                                             #
###############################################################################


#-----------------------Set Options via Function------------------------------#

#' Set Options for the DCS procedure
#' 
#' @param kernPar kernel parameters \eqn{k} and \eqn{\mu} as vector.
#' @param drv order for estimation of \eqn{(\nu_1, \nu_2)} to be estimated. The
#' polynomial order is selected as \eqn{(\nu_1 + 1, \nu_2 + 1)}. Default value 
#' is (0,0).
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
#' @section Usage:
#' \code{DCSmooth(Y, X = 1, T = 1, bndw = "auto", dcsOptions = setOptions())}
#' 
#' @section Details:
#' The function \code{summary}
#' 
#' @examples
#' y = y.norm1 + rnorm(100^2)
#' dcs(y)
#' 
#' @export
#' 

dcs = function(Y, X = 1, T = 1, bndw = "auto", dcsOptions = setOptions())
{
  time_start = Sys.time() # get starting time for performance measuring
  
  # check for correct inputs of data and options
  .dcsCheck_Y(Y)
  .dcsCheck_bndw(bndw, dcsOptions)
  .dcsCheck_options(dcsOptions)
  
  # set up vectors for X and T
  if (length(X) == 1 && length(T) == 1)
  {
    X = seq(from = 0, to = 1, length.out = dim(Y)[1])
    T = seq(from = 0, to = 1, length.out = dim(Y)[2])
  }
  
  # set kernel function (type, k, mu, nu = 0)
  nameKernFcn = paste0("MW", dcsOptions$kernPar[1], dcsOptions$kernPar[2], "0")
  kernelFcn  = kernelFcn_assign(nameKernFcn) # set kernel Function to use in optimization

  # check for given bandwidths
  if (bndw[1] == "auto") {
    bndwAuto = TRUE
  
    # bandwidth selection process
    if (all(dcsOptions$pOrder == 0))
    {
      # kernel regression
      bndwObj = KR_bndwSelect(Y, kernelFcn, dcsOptions)
      bndw = pmin(bndwObj$bndw, c(0.45, 0.45)) # KR cannot handle larger bndws.
    } else if (any(dcsOptions$pOrder > 0)) {
      # local polynomial regression
      bndwObj = LP_bndwSelect(Y, kernelFcn, dcsOptions)
      bndw = bndwObj$bndw
    } 
  } else {
    bndwAuto = FALSE
  }
  
  # estimation of resulting surface
  if (all(dcsOptions$pOrder == 0)) # kernel regression
  {
    DCSOut = KR_DoubleSmooth2(yMat = Y, hVec = bndw, drvVec = c(0, 0), 
      kernFcnPtrX = kernelFcn, kernFcnPtrT = kernelFcn)
  } else if (any(dcsOptions$pOrder > 0)) {
    # local polynomial regression
    DCSOut = LP_DoubleSmooth2(yMat = Y, hVec = bndw, polyOrderVec = 
            dcsOptions$pOrder, drvVec = dcsOptions$drv, kernFcnPtr = kernelFcn)
  }
  
  # calculate residuals
  R = Y - DCSOut
  
  time_end = Sys.time()
  
  if (bndwAuto == TRUE)
  {
    DCS_out = list(X = X, T = T, Y = Y, M = DCSOut, R = R,bndw = bndw,
                   cf = bndwObj$varCoef, iterations = bndwObj$iterations,
                   qarma_model = bndwObj$qarma_model, dcsOptions = dcsOptions,
                   timeUsed = difftime(time_end, time_start, units = "secs"))
    attr(DCS_out, "bndwAuto") = bndwAuto
  } else if (bndwAuto == FALSE) {
    # probably unadvised, that the dcs object differs according to the type
    # of bndw selection
    DCS_out = list(X = X, T = T, Y = Y, M = DCSOut, R = R, bndw = bndw,
                   dcsOptions = dcsOptions,
                   timeUsed = difftime(time_end, time_start, units = "secs"))
    attr(DCS_out, "bndwAuto") = bndwAuto
  }
  
  # apply class to output object
  class(DCS_out) = "dcs"

  return(DCS_out)
}

#-------------------Function for plotting smoothed surface--------------------#

plotDCS = function(DCSobj, ...)
{
  .persp3ddcs(DCSobj = DCSobj, ...)
}