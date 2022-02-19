#' Nonparametric Double Conditional Smoothing for 2D Surfaces
#' 
#' \code{spec_test}
#' 
#' @param Y A numeric matrix that contains the observations of the random field
#'   or functional time-series.
#' @param omega local point
#' @param lng language
#' 
#' @export
#' 

spec_test = function(Y, omega, lng)
{
  if (lng == "R")
  {
    return(specDens(Y, omega))
  }
  else if (lng == "Cpp")
  {
    return(specDens_cpp(Y, omega))
  }
}