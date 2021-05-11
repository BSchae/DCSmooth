#' Summarizing Results from Double Conditional Smoothing
#' 
#' \code{summary} method for class \code{"dcs"}
#' 
#' @export

summary.dcs = function(object)
{
  nameKernFcn = paste0("MW", object$dcsOptions$kernPar[1],
                       object$dcsOptions$kernPar[2], "0")
  
  # when automatic bandwidth selection is selected
  if (attr(object, "bndwAuto") == TRUE)
  {
    cat(class(object), "\n")
    cat("---------------------------------------------", "\n")
    cat("DCS with automatic bandwidth selection:\n")
    cat("---------------------------------------------", "\n")
    cat("Results:\n")
    cat("Estimated Bandwidths: \t h_x:\t", object$bndw[1], "\n")
    cat("\t \t \t h_t:\t", object$bndw[2], "\n")
    cat("Variance Factor: \t c_f:\t", object$cf, "\n")
    cat("Iterations:\t\t\t", object$iterations, "\n")
    cat("---------------------------------------------", "\n")
    cat("Used Parameters:\n")
    cat("Local Polynomial Order:\t", object$dcsOptions$pOrder, "\n")
    cat("Kernel Type:\t\t", nameKernFcn, "\n")
    cat("Variance Estimation:\t", object$dcsOptions$varEst, "\n")
    cat("Inflation Exponents:\t", object$dcsOptions$inflExp, "\n")
    cat("Inflation Parameters\t", object$dcsOptions$inflPar, "\n")
    cat("Const. Window\t\t", object$dcsOptions$constWindow, "\n")
    cat("Shrink Parameter:\t", object$dcsOptions$delta, "\n")
    cat("Time used (seconds):\t", object$timeUsed, "\n")
    cat("---------------------------------------------", "\n")
  
  # when given bandwidths are used.
  } else if (attr(object, "bndwAuto") == FALSE) {
    cat(class(object), "\n")
    cat("---------------------------------------------", "\n")
    cat("DCS with given bandwidths:\n")
    cat("---------------------------------------------", "\n")
    cat("Used Parameters:\n")
    cat("Bandwidths: \t h_x:\t", object$bndw[1], "\n")
    cat("\t \t \t h_t:\t", object$bndw[2], "\n")
    cat("Local Polynomial Order:\t", object$dcsOptions$pOrder, "\n")
    cat("Kernel Type:\t\t", nameKernFcn, "\n")
    cat("Const. Window\t\t", object$dcsOptions$constWindow, "\n")
    cat("---------------------------------------------", "\n")
  }
}
  
#' Print Results from Double Conditional Smoothing
#' 
#' \code{print} method for class \code{"dcs"}
#' 
#' @export

print.dcs = function(object)
{
  # when automatic bandwidth selection is selected
  if (attr(object, "bndwAuto") == TRUE)
  {
    cat(class(object), "\n")
    cat("---------------------------------------", "\n")
    cat("DCS with automatic bandwidth selection\n")
    cat("---------------------------------------", "\n")
    cat("Selected Bandwidths:\n")
    cat("\th_x:", object$bndw[1], "\n")
    cat("\th_t:", object$bndw[2], "\n")
    cat("Variance Factor:\n")
    cat("\tc_f:", object$cf, "\n")
    cat("---------------------------------------")
    
  # when given bandwidths are used.
  } else if (attr(object, "bndwAuto") == FALSE) {
    cat(class(object), "\n")
    cat("---------------------------", "\n")
    cat("DCS with given bandwidths\n")
    cat("---------------------------", "\n")
    cat("Used Bandwidths:\n")
    cat("\th_x:", object$bndw[1], "\n")
    cat("\th_t:", object$bndw[2], "\n")
    cat("---------------------------")
  }
}

#' Contour Plot for the Double Conditional Smoothing
#' 
#' @export

plot.dcs = function(object, color = c("#444C5C", "#78A5A3", "#E1B16A",
                                      "#CE5A57"))
{
  cat("Plot choices for dcs object:", fill = TRUE)
  choices <- c(1, 2, 3)
  choice_names <- c("Original observations Y:", "Smoothed Surface M:",
                    "Residuals R:")
  choices_df <- data.frame(choices)
  colnames(choices_df) <- ""
  rownames(choices_df) <- choice_names
  print.data.frame(choices_df)
  plot_choice <- readline(prompt="Please enter the corresponding number: ")
  plot_choice <- as.numeric(plot_choice)
  
  if (plot_choice == 1) {
    Y = object$Y
  } else if (plot_choice == 2) {
    Y = object$M
  } else if (plot_choice == 3) {
    Y = object$R
  } else {
    stop(plot_choice, " is not a valid plot-type.")
  }

  vecXT = base::expand.grid(object$X, object$T)
  colVec = plotDCScolFcn(Y, myColor = color)
  mainTitle = paste0("Contour plot of", plot_choice)
  graphics::plot(vecXT[, 1], vecXT[, 2], pch = 15, col = colVec,
       xlab = "T", ylab = "X", main = "Contour plot")
}

# neccessary from H. Wickhams best practices
.onUnload <- function(libpath) { library.dynam.unload("DCSmooth", libpath) }
