% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/DCS_MainFunctions.R
\name{dcs}
\alias{dcs}
\title{Nonparametric Double Conditional Smoothing for 2D Surfaces}
\usage{
dcs(Y, dcs_options = set.options(), h = "auto", parallel = FALSE, ...)
}
\arguments{
\item{Y}{A numeric matrix that contains the observations of the random field
or functional time-series.}

\item{dcs_options}{An object of class \code{"dcs_options"}, specifying the
parameters for the smoothing and bandwidth selection procedure.}

\item{h}{Bandwidth for smoothing the observations in \code{Y}. Can be a
two-valued numerical vector with bandwidths in row- and column-direction.
If the value is \code{"auto"} (the default), bandwidth selection will be 
carried out by the iterative plug-in algorithm.}

\item{parallel}{A logical value indicating if parallel computing should be
used for faster computation. Default value is \code{parallel = FALSE}.
Parallelization seems to be efficient at above 400,000 observations.}

\item{...}{Additional arguments passed to \code{dcs}. Currently supported are
numerical vectors \code{X} and/or \code{T} containing the exogenous
covariates with respect to the rows and columns.}
}
\value{
\code{dcs} returns an object of class "dcs", including
 \tabular{ll}{
 \code{Y} \tab matrix of original observations. \cr
 \code{X, T} \tab vectors of covariates over rows (\code{X}) and columns 
  (\code{T}). \cr
 \code{M} \tab resulting matrix of smoothed values. \cr
 \code{R} \tab matrix of residuals of estimation, \eqn{Y - M}. \cr
 \code{h} \tab optimized or given bandwidths. \cr
 \code{c_f} \tab estimated variance coefficient. \cr
 \code{var_est} \tab estimated variance model. If the variance function is
  modeled by an SARMA/SFARIMA, \code{var_est} is an object of class "sarma"/
  "sfarima".\cr
 \code{dcs_options} \tab an object of class \code{cds_options} containing the
  initial options of the dcs procedure. \cr
 \code{iterations} \tab number of iterations of the IPI-procedure. \cr
 \code{time_used} \tab time spend searching for optimal bandwidths (not
  overall runtime of the function). \cr
}
}
\description{
\code{dcs} provides a double conditional nonparametric smoothing of the
expectation surface of a functional time series or a random field on a
lattice. Bandwidth selection is done via an iterative plug-in method.
}
\section{Details}{

See the vignette for a more detailed description of the function.
}

\examples{
# See vignette("DCSmooth") for examples and explanation

y <- y.norm1 + matrix(rnorm(101^2), nrow = 101, ncol = 101)
dcs(y)

}
\references{
Schäfer, B. and Feng, Y. (2021). Fast Computation and Bandwidth Selection 
 Algorithms for Smoothing Functional Time Series. Working Papers CIE 143, 
 Paderborn University.
}
\seealso{
\code{\link{set.options}}
}
