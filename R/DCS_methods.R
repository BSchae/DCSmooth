#------------------------------Summary Methods---------------------------------#

#' Summarizing Results from Double Conditional Smoothing
#' 
#' @description \code{summary} method for class \code{"dcs"}
#' 
#' @param object an object of class "dcs", usually, a result of a call to 
#'  \code{\link{dcs}}.
#' @param digits the number of significant digits to use when printing.
#' 
#' @section Details:
#' \code{summary.dcs} strips an object of class \code{"dcs"} from all large
#'  matrices (\code{Y}, \code{X}, \code{T}, \code{M}, \code{R}), allowing
#'  for easier handling of meta-statistics of the bandwidth selection procedure.
#' 
#'  \code{print.summary_dcs} returns a list of summary statistics
#'  from the dcs procedure. The output depends on the use of the \code{dcs}-
#'  function. If automatic bandwidth selection is chosen, \code{summary.dcs}
#'  prints detailed statistics of the type of regression, the estimated 
#'  bandwidths \code{h_x}, \code{h_t}, the variance coefficient \code{c_f} and 
#'  performance statistics such as the number of iterations of the IPI-algorithm 
#'  and the time used for bandwidth selection.
#' 
#'  The method used for estimation of the variance coefficient is printed and 
#'  the results of an QARMA-estimation, if available.
#' 
#'  If bandwidths are supplied to \code{dcs}, \code{summary.dcs} only prints
#'  the given bandwidths.
#' 
#' @return The function \code{summary.dcs} returns an object of class \code{summary_dcs}
#'  including \tabular{ll}{
#'  \code{h_opt} \tab estimated optimal bandwidth from the IPI-procedure. \cr
#'  \code{c_f} \tab estimated variance factor. \cr
#'  \code{iterations} \tab number of iterations of the IPI-procedure. \cr
#'  \code{time_used} \tab time spend searching for optimal bandwidths (not overall
#'   runtime of the function). \cr
#'  \code{qarma} \tab optional return, if method \code{"qarma"} is chosen for estimation
#'   of the variance factor. Omitted, if \code{"iid"} is used. \cr
#'  \code{dcs_options} \tab an object of class \code{cds_options} containing the initial
#'   options of the dcs procedure. \cr
#' }
#' 
#' @examples
#' y <- y.norm1 + matrix(rnorm(101^2), nrow = 101, ncol = 101)
#' dcs_object <- dcs(y)
#' summary(dcs_object)
#' 
#' @name summary.dcs
#' 
#' @export
#' 

summary.dcs = function(object)
{
  if (exists("qarma", object)) {
    summary_dcs = list(h = object$h, c_f = object$c_f,
                       iterations = object$iterations, 
                       time_used = object$time_used, qarma = object$qarma,
                       dcs_options = object$dcs_options)
  } else {
    summary_dcs = list(h = object$h, c_f = object$c_f,
                       iterations = object$iterations, 
                       time_used = object$time_used,
                       dcs_options = object$dcs_options)
  }
  
  attr(summary_dcs, "h_select_auto") = attr(object, "h_select_auto")
                     
  class(summary_dcs) = "summary_dcs"
  return(summary_dcs)
}

#----------------------------Print (Summary) Methods---------------------------#


#'
#' @rdname summary.dcs
#' 
#' @export

print.summary_dcs = function(object, digits = max(3, getOption("digits") - 3))
{
  name_kern_fcn = paste0("MW", object$dcs_options$kern_par[1],
                       object$dcs_options$kern_par[2], "0")
  if (object$dcs_options$type == "KR")
  {
    reg_type = paste0("kernel regression")
  } else if (object$dcs_options$type == "LP") {
    reg_type = paste0("local polynomial regression")
  }
  
  # when automatic bandwidth selection is selected
  if (attr(object, "h_select_auto") == TRUE)
  {
    cat(class(object), "\n")
    cat("------------------------------------------", "\n")
    cat("DCS with automatic bandwidth selection:\n")
    cat("------------------------------------------", "\n")
    cat("Results of ", reg_type, ":\n", sep = "")
    cat("Estimated Bandwidths: \t h_x:\t", signif(object$h[1], 
                                                  digits = digits), "\n")
    cat("\t \t \t h_t:\t", signif(object$h[2], digits = digits), "\n")
    cat("Variance Factor: \t c_f:\t", signif(object$c_f, digits = digits), "\n")
    cat("Iterations:\t\t\t", object$iterations, "\n")
    cat("Time used (seconds):\t \t", signif(object$time_used, 
                                           digits = digits), "\n")
    cat("------------------------------------------", "\n")
    if (object$dcs_options$var_est == "iid")
    {
      cat("Variance Estimation: \t iid.\n")
    } else if (object$dcs_options$var_est == "qarma") {
      cat("Variance Estimation: \t qarma\n")
      cat("Estimated Parameters:\n")
      cat("ar:\n")
      print(signif(object$qarma$ar, digits = 5))
      cat("ma:\n")
      print(signif(object$qarma$ma, digits = 5))
    }
    cat("------------------------------------------", "\n")
    cat("See used parameter with \"$dcs_options\".\n")
    cat("------------------------------------------", "\n")
    
  # when given bandwidths are used.
  } else if (attr(object, "h_select_auto") == FALSE) {
    cat(class(object), "\n")
    cat("------------------------------------------", "\n")
    cat("DCS with given bandwidths:\n")
    cat("------------------------------------------", "\n")
    cat("Used bandwidths: \t h_x:\t", object$h[1], "\n")
    cat("\t \t \t h_t:\t", object$h[2], "\n")
    cat("------------------------------------------", "\n")
    cat("See used parameter with \"$dcs_options\".\n")
    cat("------------------------------------------", "\n")
  }
}
  
#' Print Results from Double Conditional Smoothing
#' 
#' @description \code{print} method for class \code{"dcs"}
#' 
#' @section Details:
#' \code{print.dcs} prints a short summary of an object of class \code{dcs},
#' only including bandwidths and the estimated variance coefficient (only if
#' automatic bandwidth selection is used).
#' 
#' @param object an object of class "dcs", usually, a result of a call to 
#' \code{\link{dcs}}.
#' 
#' @seealso \code{\link{plot.dcs}}, \code{\link{print.dcs_options}}
#' 
#' @examples
#' y <- y.norm1 + matrix(rnorm(101^2), nrow = 101, ncol = 101)
#' dcs_object <- dcs(y)
#' print(dcs_object)
#' dcs_object
#' 
#' @export
#' 

print.dcs = function(object)
{
  # when automatic bandwidth selection is selected
  if (attr(object, "h_select_auto") == TRUE)
  {
    cat(class(object), "\n")
    cat("--------------------------------------", "\n")
    cat("DCS with automatic bandwidth selection\n")
    cat("--------------------------------------", "\n")
    cat(" Selected Bandwidths:\n")
    cat("\t\th_x:", signif(object$h[1], digits = 5), "\n")
    cat("\t\th_t:", signif(object$h[2], digits = 5), "\n")
    cat(" Variance Factor:\n")
    cat("\t\tc_f:", signif(object$c_f, digits = 5), "\n")
    cat("--------------------------------------", "\n")
    
  # when given bandwidths are used.
  } else if (attr(object, "h_select_auto") == FALSE) {
    cat(class(object), "\n")
    cat("---------------------------", "\n")
    cat("DCS with given bandwidths\n")
    cat("---------------------------", "\n")
    cat("Used Bandwidths:\n")
    cat("\th_x:", object$h[1], "\n")
    cat("\th_t:", object$h[2], "\n")
    cat("---------------------------")
  }
}

#' Print Options for Double Conditional Smoothing
#' 
#' @description \code{print} method for class \code{"dcs_options"}
#' 
#' @section Details:
#' \code{print.dcs_options} prints the obptions used for the \code{\link{dcs}} 
#' function, which should be an object of class \code{"dcs_options"}.
#' 
#' @param object an object of class "dcs_options", usually, a result of a call to 
#' \code{\link{set_options}}.
#' 
#' @seealso \code{\link{print.dcs}}
#' 
#' @examples
#' ## Default options
#' myOpt <- set.options()
#' print(myOpt)
#' 
#' ## Use Kernel regression
#' myOpt <- set.options(type = "KR")
#' print(myOpt)
#' 
#' @export
#' 

print.dcs_options = function(object)
{
  cat(class(object), "\n")
  cat("---------------------------------------", "\n")
  cat("options for DCS \t rows \t cols \n")
  cat("---------------------------------------", "\n")
  
  # when Kernel Regression is selected
  if (object$type == "KR")
  {
    cat("type: kernel regression \n")
    cat("kernel order: \t \t ", object$kern_par[1], "\t", object$kern_par[2], "\n")
    cat("derivative: \t \t ", object$drv[1], "\t", object$drv[2], "\n")
    cat("inflation exponent: \t ", object$infl_exp[1],
        "\t", object$infl_exp[2], "\n")
    cat("inflation parameters: \t ", object$infl_par[1],
        "\t", object$infl_par[2], "\n")
    cat("boundary factor: \t ", object$delta[1], "\t", object$delta[2], "\n")
    cat("constant window: \t ", object$const_window, "\n")
    cat("---------------------------------------", "\n")
  } else if (object$type == "LP") {  
  # when Local Polynomial Regression is selected
    cat("type: local polynomial regression \n")
    cat("kernel order: \t \t ", object$kern_par[1], "\t", object$kern_par[2], "\n")
    cat("derivative: \t \t ", object$drv[1], "\t", object$drv[2], "\n")
    cat("polynomial order: \t ", object$p_order[1], "\t", object$p_order[2], "\n")
    cat("inflation exponent: \t ", object$infl_exp[1],
        "\t", object$infl_exp[2], "\n")
    cat("inflation parameters: \t ", object$infl_par[1],
        "\t", object$infl_par[2], "\n")
    cat("boundary factor: \t ", object$delta[1], "\t", object$delta[2], "\n")
    cat("constant window: \t ", object$const_window, "\n")
    cat("---------------------------------------", "\n")
  }
  if (object$var_est == "iid")
  {
    cat("variance estimation: \t  iid. \n")
  } else if (object$var_est == "qarma")
  {
    cat("variance estimation: \t  qarma \n")
    if (!(is.list(object$qarma_order)) && 
        object$qarma_order %in% c("gpac", "bic"))
    {
      cat("order selection: \t ", object$qarma_order, "\n")
    } else {
      cat("model order: \t ar: \t ", object$qarma_order$ar[1], 
          "\t", object$qarma_order$ar[2], "\n")
      cat("\t \t ma: \t ", object$qarma_order$ma[1], 
          "\t", object$qarma_order$ma[2], "\n")
    }
  }
  cat("---------------------------------------", "\n")
}

#----------------------------------Other Methods-------------------------------#

#' Residuals of "dcs"-object
#' 
#' @description Returns the residuals of an object of class \code{"dcs"}.
#' 
#' @param dcs_object an object of class \code{"dcs"}, usually the result of a call to \code{\linl{dcs}}
#' 
#' @seealso \code{\link{dcs}}
#' 
#' @example
#' y = y.norm1 + matrix(rnorm(101^2), nrow = 101, ncol = 101)
#' dcs_object = dcs(y)
#' residuals(dcs_object)
#' 

residuals.dcs = function(dcs_object)
{
  return(dcs_object$R)
}

#' Contour Plot for the Double Conditional Smoothing
#' 
#' @description \code{plot} method for class \code{"dcs"}
#' 
#' @section Details:
#' \code{plot.dcs} provides a contour plot of either the original data (1),
#'  smoothed surface (2) or residuals (3).
#' 
#' @param object an object of class "dcs_options", usually, a result of a call to 
#' \code{\link{set_options}}.
#' @param plot_choice override the prompt to specify a plot, can be 
#'  \code{c(1, 2, 3)}.
#' 
#' @seealso \code{\link{dcs.3d}} to plot the surface.
#' 
#' @examples
#' ## Contour plot of smoothed surface
#' y <- y.norm1 + matrix(rnorm(101^2), nrow = 101, ncol = 101)
#' dcs_object <- dcs(y)
#' plot(dcs_object, plot_choice = 2)
#' 
#' @export
#' 

plot.dcs = function(object, plot_choice = "choice", ...)
{
  fcn_arg = list(...)
  
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
    plot_choice <- readline(prompt="Please enter the corresponding number: ")
    plot_choice <- as.numeric(plot_choice)
  } else if (!(plot_choice %in% 1:3)) {
    stop("Invalid value in argument \"plot_choice\". Use c(1, 2, 3).")
  }
  
  if (plot_choice == 1) {
    Y = object$Y
  } else if (plot_choice == 2) {
    Y = object$M
  } else if (plot_choice == 3) {
    Y = object$R
  } else {
    stop(plot_choice, " is not a valid plot-type.")
  }
  
  if (exists("color", fcn_arg))
  {
    color = fcn_arg$color
  } else {
    color = c("#444C5C", "#78A5A3", "#E1B16A", "#CE5A57")
  }
  
  vec_XT = base::expand.grid(object$X, object$T)
  col_vec = plot.dcs.colors(Y, color = color)
  main_title = paste("Contour plot of ", choice_names[plot_choice])
  graphics::plot(vec_XT[, 1], vec_XT[, 2], pch = 15, col = col_vec,
       xlab = "T", ylab = "X", main = main_title)
}

# necessary from H. Wickhams best practices
.onUnload <- function(libpath) { library.dynam.unload("DCSmooth", libpath) }
