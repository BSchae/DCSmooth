################################################################################
#                                                                              #
#                      DCSmooth Package: User Functions                        #
#                                                                              #
################################################################################

# Includes the main functions for the DCS package

# set.options (exported)

# dcs (exported)


#------------------------Set Options via Function------------------------------#

#' Set Options for the DCS procedure
#' 
#' @param type either local polynomial regression (\code{"LP"}, the default) or 
#'  kernel regression (\code{"KR"}).
#' @param kerns a character vector of length 2 containing the identifier for the
#'  kernels to be used in kernel regression. Weighting functions in local
#'  polynomial regression are computed according to the identifier. Default value
#'  is \code{MW_220}, the Mueller-Wang kernel of order \eqn{(2, 2, 0)}.
#' @param drv A non-negative vector of length 2, containing the derivative
#'  orders to be estimated from the given data. The default is \code{c(0, 0)}. For
#'  LP-regression, polynomial order is selected as \eqn{(\nu_1 + 1, \nu_2 + 1)}.
#' @param var_model the method of estimating the variance coefficient \eqn{c_f}. 
#'  Currently available are \code{var_model = c("iid", "sarma_HR", "sarma_sep",
#'  "sarma_RSS", "sfarima_RSS")}. Replacing the argument \code{var_model}. For
#'  code using \code{var_est}, the argument is converted to \code{var_model}.
#' @param ... Additional arguments passed to \code{set.options()}. This includes
#'  \code{IPI_options}, a list containing further options used by the iterative
#'  plug-in algorithm.
#' 
#' @section Details:
#' This function is used to set the options for bandwidth selection in the 
#' \code{dcs} function.
#' Detailed information can be found in the vignette.
#' 
#' @return  An object of class \code{"dcs_options"}.
#' 
#' @seealso \code{\link{dcs}}
#' 
#' @export
#' 
#' @examples
#' # See vignette("DCSmooth") for examples and explanation
#' 
#' set.options()
#' 
#' myOpt <- set.options(type = "KR", var_model = "iid")
#' y <- y.norm1 + matrix(rnorm(101^2), nrow = 101, ncol = 101)
#' dcs(y, dcs_options = myOpt)

set.options <- function(    # inside function with default values in arguments
  # standard options
  type      = "LP",        # either "LP" for local polynomial regression or
  # "KR" for kernel regression
  kerns     = c("MW_220", "MW_220"), # choose a kernel function
  drv       = c(0, 0),
  var_model = "iid",
  
  # advanced options in ellipsis
  ...
)
{
  # get ellipsis
  add_options <- list(...)
  # include old var_est options
  if (exists("var_est", add_options))
  {
    message("Note: option \"var_est\" is deprecated, argument is converted to ",
            "\"var_model\" automatically.")
    if (add_options$var_est == "iid") { var_model = "iid" }
    if (add_options$var_est == "qarma") { var_model = "sarma_HR" }
    if (add_options$var_est == "sarma") { var_model = "sarma_sep" }
    if (add_options$var_est == "sarma2") { var_model = "sarma_RSS" }
    if (add_options$var_est == "lm") { var_model = "sfarima_RSS" }
    if (add_options$var_est %in% c("qarma_gpac", "qarma_bic"))
    {
      add_options$var_model = "sarma_HR"
      warning("For automatic order selection, use \"model_order\" in ",
              "\"dcs()\".")
    }  
  }
  if (exists("IPI_options", add_options))
  {
    IPI_options = add_options$IPI_options
  } else {
    IPI_options = list()
  }
  
  # check if inputs are vectors
  if (length(kerns) == 1) { kerns <- c(kerns, kerns) }
  if (length(drv) == 1) { drv <- c(drv, drv) }
  
  exception.check.options.input(type, kerns, drv, var_model, IPI_options)
  
  # Select options according to type ("LP", "KR")
  if (type == "LP")
  {
    p_order <- drv + 1
    if (!exists("infl_exp", IPI_options))
    {
      IPI_options$infl_exp <- c("auto", " ")
    }
    if (!exists("infl_par", IPI_options)) { IPI_options$infl_par <- c(1, 1) }
    if (!exists("delta", IPI_options)) { IPI_options$delta <- c(0.05, 0.05) }
    if (!exists("const_window", IPI_options))
    {
      IPI_options$const_window <- FALSE
    }
  } else if (type == "KR") {
    p_order <- NA
    if (!exists("infl_exp", IPI_options))
    {
      IPI_options$infl_exp <- c(0.5, 0.5)
    }
    if (!exists("infl_par", IPI_options)) { IPI_options$infl_par <- c(2, 1) }
    if (!exists("delta", IPI_options)) { IPI_options$delta <- c(0.05, 0.05) }
    if (!exists("const_window", IPI_options))
    {
      IPI_options$const_window <- FALSE
    }
  } else {
    stop("Unknown type \"", type, "\"")
  }
  
  # change IPI_options for long memory estimation
  if (!exists("infl_exp", IPI_options) && var_model == "sfarima_RSS")
  {
    IPI_options$infl_exp <- c(0.5, 0.5)
  }
  if (!exists("infl_par", IPI_options) && var_model == "sfarima_RSS")
  {
    IPI_options$infl_par <- c(3, 1)
  }
  
  options_list <- list(type = type, kerns = kerns, p_order = p_order,
                       drv = drv, var_model = var_model,
                       IPI_options = IPI_options)
  
  # apply class to output object
  class(options_list) <- "dcs_options"
  
  exception.check.options(options_list)
  
  return(options_list)
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
#' @param dcs_options An object of class \code{"dcs_options"}, specifying the
#'  parameters for the smoothing and bandwidth selection procedure.
#' @param h Bandwidth for smoothing the observations in \code{Y}. Can be a
#'  two-valued numerical vector with bandwidths in row- and column-direction.
#'  If the value is \code{"auto"} (the default), bandwidth selection will be 
#'  carried out by the iterative plug-in algorithm.
#' @param ... Additional arguments passed to \code{dcs}. These might include
#'  numerical vectors \code{X} and/or \code{T} containing the exogenous
#'  covariates with respect to the rows and columns. If the error term model in 
#'  \code{set.options()$var_model} is any of the available parametric models
#'  (SARMA, SFARIMA), \code{model_order} is either an optional list containing
#'  the two-dimensional model orders in the form
#'  \code{list(ar = c(1, 1), ma = c(1, 1)} or a character string specifying
#'  an order selection criterion by \code{c("AIC", "BIC", "gpac")}. In the case
#'  of automatic order selection, \code{order_max} is an optional list
#'  containing the two-dimensional maximum orders for the order selection in the
#'  form \code{list(ar = c(1, 1), ma = c(1, 1)}.
#' 
#' @return \code{DCSmooth} returns an object of class "dcs", including
#'  \tabular{ll}{
#'  \code{Y} \tab matrix of original observations. \cr
#'  \code{X, T} \tab vectors of covariates over rows (\code{X}) and columns 
#'   (\code{T}). \cr
#'  \code{M} \tab resulting matrix of smoothed values. \cr
#'  \code{R} \tab matrix of residuals of estimation, \eqn{Y - M}. \cr
#'  \code{h} \tab optimized or given bandwidths. \cr
#'  \code{c_f} \tab estimated variance coefficient. \cr
#'  \code{dcs_options} \tab an object of class \code{cds_options} containing the
#'   initial options of the dcs procedure. \cr
#'  \code{iterations} \tab number of iterations of the IPI-procedure. \cr
#'  \code{time_used} \tab time spend searching for optimal bandwidths (not
#'   overall runtime of the function). \cr
#'  \code{qarma} \tab optional return, if method \code{"qarma"} is chosen for 
#'   estimation of the variance factor. Omitted, if \code{"iid"} is used.
#' }
#' 
#' @section Details:
#' See the vignette for a more detailed description of the function.
#' 
#' @seealso \code{\link{set.options}}
#' 
#' @examples 
#' # See vignette("DCSmooth") for examples and explanation
#' 
#' y <- y.norm1 + matrix(rnorm(101^2), nrow = 101, ncol = 101)
#' dcs(y)
#' 
#' @export
#' 

dcs <- function(Y, dcs_options = set.options(), h = "auto", ...)
{
  #-------Front Matter-------#
  
  # check for correct inputs of data and options
  exception.check.Y(Y)
  exception.check.bndw(h, dcs_options)
  exception.check.options(dcs_options)
  n_x = dim(Y)[1]; n_t = dim(Y)[2]
  
  # process ellipsis arguments from (...)
  {
    # set up vectors for X and T if neccessary
    args_list <- list(...)
    add_options <- list() # additional options depending on type and var_model
    
    if (!exists("X", where = args_list))
    {
      X <- 0:(n_x - 1)/(n_x - 1)
    } else {
      X <- args_list$X
    }
    if (!exists("T", where = args_list))
    {
      T <- 0:(n_t - 1)/(n_t - 1)
    } else {
      T <- args_list$T
    }
    exception.check.XT(Y, X, T)
    
    # model order
    if (exists("model_order", where = args_list))
    {
      exception.check.model_order(args_list$model_order, dcs_options)
      add_options$model_order = args_list$model_order
    } else if (dcs_options$var_model %in% dcs_list_var_model &&
               dcs_options$var_model != "iid") {
      add_options$model_order = list(ar = c(1, 1), ma = c(1, 1))
    }
    if (exists("order_max", where = args_list) &&
        length(add_options$model_order) == 1 &&
        !is.list(add_options$model_order) &&
        add_options$model_order %in% c("aic", "bic", "gpac"))
    {
      exception.check.order_max(args_list$order_max)
      add_options$order_max = args_list$order_max
    } else if (length(add_options$model_order) == 1 &&
               !is.list(add_options$model_order) &&
               add_options$model_order %in% c("aic", "bic", "gpac")) {
      add_options$order_max = list(ar = c(1, 1), ma = c(1, 1))
    }
  }

  #-------Bandwidth Selection-------#
  
  # check if given bandwidths exist
  if (length(h) == 1 && h[1] == "auto")
  {
    h_select_auto <- TRUE
    
    time_start <- Sys.time() # get starting time for performance measuring
  
    # bandwidth selection process
    if (dcs_options$type == "KR")
    {
      ### kernel regression
      # TODO: allow for derivatives in kernel regression
      h_select_obj <- KR.bndw(Y, dcs_options, add_options)
      h_opt <- pmin(h_select_obj$h_opt, c(0.45, 0.45)) 
                                          # KR cannot handle larger bandwidths
    } else if (dcs_options$type == "LP") {
      ### local polynomial regression
      h_select_obj <- LP.bndw(Y, dcs_options, add_options)
      h_opt <- h_select_obj$h_opt
    }
    
    if (h_select_obj$var_model$stnry == FALSE)
    {
      warning("error model is not stationary.")
    }
    
    time_end <- Sys.time()
    
  } else {
    h_opt <- h
    h_select_auto <- FALSE
    time_end <- time_start <- Sys.time()
  }
  
  #-------Estimation of Result-------#
  
  if (dcs_options$type == "KR") # kernel regression
  {
    ### prepare kernel functions
      kernel_x <- kernel_fcn_assign(dcs_options$kerns[1])
      kernel_t <- kernel_fcn_assign(dcs_options$kerns[2])
      
    ### check bandwidth
      if (any(h_opt > c(0.45, 0.45)))
      {
        h_opt = pmin(h_opt, c(0.45, 0.45))
        warning("Bandwidth h too large for \"KR\", changed to largest ",
                "working value.")
      }
    
    Y_smth_out <- KR_dcs_const0(yMat = Y, hVec = h_opt, drvVec = c(0, 0), 
                            kernFcnPtrX = kernel_x,
                            kernFcnPtrT = kernel_t)
  } else if (dcs_options$type == "LP") {     # local polynomial regression
    ### prepare weight functions
      # set variables for weight type
      kern_type_vec = sub("^([[:alpha:]]*).*", "\\1", dcs_options$kern)
      # extract weighting type
      mu_vec = as.numeric(substr(dcs_options$kern, 
                                 nchar(dcs_options$kern) - 1,
                                 nchar(dcs_options$kern) - 1))
      # extract kernel parameter mu
      weight_x = weight_fcn_assign(kern_type_vec[1])
      weight_t = weight_fcn_assign(kern_type_vec[2])
    
    ### check bandwidth
      if (any(h_opt < (dcs_options$p_order + 2)/(c(n_x, n_t) - 1)))
      {
        h_opt = pmax(h_opt, (dcs_options$p_order + 2)/(c(n_x, n_t) - 1))
        warning("Bandwidth h too small for \"LP\", changed to smallest ",
                "working value.")
      }
      
    Y_smth_out <- LP_dcs_const0_BMod(yMat = Y, hVec = h_opt, polyOrderVec = 
                            dcs_options$p_order, drvVec = dcs_options$drv,
                            muVec = mu_vec, weightFcnPtr_x = weight_x,
                            weightFcnPtr_t = weight_t)
  }
  
  # calculate residuals and estimate model
  R <- Y - Y_smth_out
  var_model_est <- cf.estimation(R, dcs_options, add_options)
  
  if (h_select_auto == TRUE)
  {
    dcs_out <- list(Y = Y, X = X, T = T, M = Y_smth_out, R = R, h = h_opt,
                   c_f = h_select_obj$var_coef, 
                   var_est = var_model_est$model_est,
                   dcs_options = dcs_options,
                   iterations = h_select_obj$iterations,
                   time_used = difftime(time_end, time_start, units = "secs"))
    attr(dcs_out, "h_select_auto") <- h_select_auto
  } else if (h_select_auto == FALSE) {
    dcs_out <- list(X = X, T = T, Y = Y, M = Y_smth_out, R = R, h = h_opt,
                   c_f = NA, var_est = var_model_est$model_est,
                   dcs_options = dcs_options,
                   iterations = NA, time_used = NA)
    attr(dcs_out, "h_select_auto") <- h_select_auto
  }
  
  # apply class to output object
  class(dcs_out) <- "dcs"

  return(dcs_out)
}

#-------------------Function for plotting smoothed surface--------------------#

#' 3D Surface Plot of "dcs"-object or numeric matrix
#' 
#' @section Details:
#' \code{surface.dcs} uses the plotly device to plot the 3D surface of the given
#' \code{"dcs"}-object or matrix. If a "dcs"-object is passed to the function, 
#'  it can be chosen between plots of the original data (1), smoothed surface 
#'  (2) and residuals (3).
#' @param Y an object of class \code{"dcs"} or a numeric matrix that contains the
#'   values to be plotted.
#' @param plot_choice override the prompt to specify a plot, can be 
#'  \code{c(1, 2, 3)}.
#' @param ... optional arguments passed to the plot function.
#' 
#' @return \code{dcs.3d} returns an object of class "plotly" and "htmlwidget".
#' 
#' @seealso \code{\link{plot.dcs}}
#' 
#' @examples
#' # See vignette("DCSmooth") for examples and explanation
#' 
#' smth =  dcs(y.norm1 + rnorm(101^2))
#' surface.dcs(smth, plot_choice = 2)
#' 
#' @export
#' 

surface.dcs <- function(Y, plot_choice = "choice", ...)
{
  if (class(Y)[1] == "dcs")
  {
    fcn_arg <- list(...)
    
    if (plot_choice == "choice")
    {
      cat("Plot choices for dcs object:", fill = TRUE)
      choices <- c(1, 2, 3)
      choice_names <- c("original observations", "smoothed surface",
                        "residuals")
      choices_df <- data.frame(choices)
      colnames(choices_df) <- ""
      rownames(choices_df) <- choice_names
      print.data.frame(choices_df)
      plot_choice <- readline("Please enter the corresponding number: ")
      plot_choice <- as.numeric(plot_choice)
    } else if (!(plot_choice %in% 1:3)) {
      stop("Invalid value in argument \"plot_choice\". Use c(1, 2, 3).")
    }
    
    if (plot_choice == 1) {
      Y_mat <- Y$Y
    } else if (plot_choice == 2) {
      Y_mat <- Y$M
    } else if (plot_choice == 3) {
      Y_mat <- Y$R
    } else {
      stop(plot_choice, " is not a valid plot-type.")
    }
    
    .plotly.3d(Y = Y_mat, X = Y$X, T = Y$T, 
              color = c("#444C5C", "#78A5A3", "#E1B16A", "#CE5A57"), ...)
  } else {
    .plotly.3d(Y = Y, ...)
  }
}