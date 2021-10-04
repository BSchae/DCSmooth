#' Single Gaussian Peak
#' 
#' Example data for using the DCSmooth functions. Data resembles a single
#'   gaussian peak on the interval \eqn{[0,1] \times [0,1]}{[0,1] x [0,1]}
#'   with maximum at \eqn{(0.5, 0.5)} and variance matrix \eqn{0.1 \cdot
#'  {\bf I}}{0.05 I}, where 
#'   \eqn{{\bf I}}{I} represents the \eqn{2 \times 2}{2 x 2} identity matrix.
#'   
#' @format A numeric matrix with 101 rows and 101 columns.
"y.norm1"

#' Double Gaussian Peak
#' 
#' Example data for using the DCSmooth functions. Data resembles two gaussian
#'   peaks on the interval \eqn{[0,1] \times [0,1]}{[0,1] x [0,1]} with maxima
#'   at \eqn{(0.5, 0.3)} with variance matrix \eqn{0.1 \cdot {\bf I}}{0.1 I}
#'   and at \eqn{(0.2, 0.8)} with variance matrix \eqn{0.05 \cdot {\bf I}}{
#'   0.05 I}, where \eqn{{\bf I}}{I} represents the \eqn{2 \times 2}{2 x 2}
#'   identity matrix.
#'   
#' @format A numeric matrix with 101 rows and 101 columns.
"y.norm2"

#' Double Gaussian Ridges
#' 
#' Example data for using the DCSmooth functions. Data resembles two gaussian
#'   ridges on the interval \eqn{[0,1] \times [0,1]}{[0,1] x [0,1]} with maxima
#'   at \eqn{(0.25, 0.75)} with variance matrix
#'   \eqn{(0.01, -0.1) \cdot {\bf I}}{(0.01, -0.1) I} and at
#'   \eqn{(0.75, 0.5)} with variance matrix
#'   \eqn{(0.01, -0.1) \cdot {\bf I}}{(0.01, -0.1) I}, where 
#'   \eqn{{\bf I}}{I} represents the \eqn{2 \times 2}{2 x 2} identity matrix.
#'   
#' @format A numeric matrix with 101 rows and 101 columns.
"y.norm3"

#' Temperatures from Yuma, AZ
#' 
#' This dataset contains the 5-minute observations of the 2020 temperature in
#' Yuma, AZ. The data is from the U.S. Climate Reference Network database at 
#' \url{https://www.ncdc.noaa.gov/crn/qcdatasets.html}. (see Diamond et. al.
#' (2013), \href{https://journals.ametsoc.org/view/journals/bams/94/4/bams-d-12-00170.1.xml}{doi:10.1175/BAMS-D-12-00170.1}).
#' The observations were adjusted matrix form for direct use with the 
#' \code{DCSmooth}-functions.
#' 
#' @format A numeric matrix with 366 rows and 288 columns containing the
#' temperatures in Celsius.
"temp.yuma"

#' Temperatures from Nunn, CO
#' 
#' This dataset contains the 5-minute observations of the 2020 temperature in
#' Nunn, CO. The data is from the U.S. Climate Reference Network database at 
#' \url{https://www.ncdc.noaa.gov/crn/qcdatasets.html}. (see Diamond et. al.
#' (2013), \href{https://journals.ametsoc.org/view/journals/bams/94/4/bams-d-12-00170.1.xml}{doi:10.1175/BAMS-D-12-00170.1}).
#' The observations were adjusted to matrix form for direct use with the 
#' \code{DCSmooth}-functions.
#' 
#' @format A numeric matrix with 366 rows and 288 columns containing the
#' temperatures in Celsius.
"temp.nunn"

#' Wind Speed from Yuma, AZ
#' 
#' This dataset contains the 5-minute observations of the 2020 wind speed in
#' Yuma, AZ. The data is from the U.S. Climate Reference Network database at 
#' \url{https://www.ncdc.noaa.gov/crn/qcdatasets.html}. (see Diamond et. al.
#' (2013), \href{https://journals.ametsoc.org/view/journals/bams/94/4/bams-d-12-00170.1.xml}{doi:10.1175/BAMS-D-12-00170.1}).
#' The observations were adjusted matrix form for direct use with the 
#' \code{DCSmooth}-functions.
#' 
#' @format A numeric matrix with 366 rows and 288 columns containing the
#' wind speeds in \eqn{m/s}.
"wind.yuma"

#' Wind Speed from Nunn, CO
#' 
#' This dataset contains the 5-minute observations of the 2020 wind speed in
#' Nunn, CO. The data is from the U.S. Climate Reference Network database at 
#' \url{https://www.ncdc.noaa.gov/crn/qcdatasets.html}. (see Diamond et. al.
#' (2013), \href{https://journals.ametsoc.org/view/journals/bams/94/4/bams-d-12-00170.1.xml}{doi:10.1175/BAMS-D-12-00170.1}).
#' The observations were adjusted to matrix form for direct use with the 
#' \code{DCSmooth}-functions.
#' 
#' @format A numeric matrix with 366 rows and 288 columns containing the
#' wind speed in \eqn{m/s}.
"wind.nunn"

#' Solar Radiation from Yuma, AZ
#' 
#' This dataset contains the 5-minute observations of the 2020 solar radiation
#' in Yuma, AZ. The data is from the U.S. Climate Reference Network database at 
#' \url{https://www.ncdc.noaa.gov/crn/qcdatasets.html}. (see Diamond et. al.
#' (2013), \href{https://journals.ametsoc.org/view/journals/bams/94/4/bams-d-12-00170.1.xml}{doi:10.1175/BAMS-D-12-00170.1}).
#' The observations were adjusted matrix form for direct use with the 
#' \code{DCSmooth}-functions.
#' 
#' @format A numeric matrix with 366 rows and 288 columns containing the
#' solar radiation in \eqn{W/m^2}.
"sun.yuma"

#' Solar Radiation from Nunn, CO
#' 
#' This dataset contains the 5-minute observations of the 2020 solar radiation
#' in Nunn, CO. The data is from the U.S. Climate Reference Network database at 
#' \url{https://www.ncdc.noaa.gov/crn/qcdatasets.html}. (see Diamond et. al.
#' (2013), \href{https://journals.ametsoc.org/view/journals/bams/94/4/bams-d-12-00170.1.xml}{doi:10.1175/BAMS-D-12-00170.1}).
#' The observations were adjusted to matrix form for direct use with the 
#' \code{DCSmooth}-functions.
#' 
#' @format A numeric matrix with 366 rows and 288 columns containing the
#' solar radiation in \eqn{W/m^2}.
"sun.nunn"

#' #' Returns of Allianz SE
#' #' 
#' #' The (log-) returns of the shares of the German insurance company Allianz SE 
#' #' from 2004-01-02 to 2014-09-30 aggregated to per-minute observations. The data
#' #' is adjusted to matrix form for direct use with the \code{DCSmooth}-functions.
#' #' 
#' #' @format A numeric matrix with 2739 rows representing the days and 510 columns
#' #' representing the intraday time points.
#' "returns.alv"
#' 
#' #' Returns of BMW AG
#' #' 
#' #' The (log-) returns of shares of the German car manufacturer BMW AG from 
#' #' 2004-01-02 to 2014-09-30 aggregated to per-minute observations. The data is 
#' #' adjusted to matrix form for direct use with the \code{DCSmooth}-functions.
#' #' 
#' #' @format A numeric matrix with 2739 rows representing the days and 510 columns
#' #' representing the intraday time points.
#' "returns.bmw"
#' 
#' #' Returns of Siemens AG
#' #' 
#' #' The (log-) returns of shares of the German technology corporation Siemens AG 
#' #' from 2004-01-02 to 2014-09-30 aggregated to per-minute observations. The data
#' #' is adjusted to matrix form for direct use with the \code{DCSmooth}-functions.
#' #' 
#' #' @format A numeric matrix with 2739 rows representing the days and 510 columns
#' #' representing the intraday time points.
#' "returns.sie"
#' 
#' #' Volumes of Allianz SE
#' #' 
#' #' The trading volumes of the shares of the German insurance company Allianz SE 
#' #' from 2004-01-02 to 2014-09-30 aggregated to per-minute observations. The data
#' #' is adjusted to matrix form for direct use with the \code{DCSmooth}-functions.
#' #' 
#' #' @format A numeric matrix with 2739 rows representing the days and 510 columns
#' #' representing the intraday time points.
#' "volumes.alv"
#' 
#' #' Volumes of BMW AG
#' #' 
#' #' The trading volumes of shares of the German car manufacturer BMW AG from 
#' #' 2004-01-02 to 2014-09-30 aggregated to per-minute observations. The data is 
#' #' adjusted to matrix form for direct use with the \code{DCSmooth}-functions.
#' #' 
#' #' @format A numeric matrix with 2739 rows representing the days and 510 columns
#' #' representing the intraday time points.
#' "volumes.bmw"
#' 
#' #' Volumes of Siemens AG
#' #' 
#' #' The trading volumes of shares of the German technology corporation Siemens AG 
#' #' from 2004-01-02 to 2014-09-30 aggregated to per-minute observations. The data
#' #' is adjusted to matrix form for direct use with the \code{DCSmooth}-functions.
#' #' 
#' #' @format A numeric matrix with 2739 rows representing the days and 510 columns
#' #' representing the intraday time points.
#' "volumes.sie"
