summary.dcs = function(object)
  nameKernFcn = paste0("MW", object$dcsOptions$kernPar[1],
                       object$dcsOptions$kernPar[2], "0")
  cat(class(object), "\n")
  cat("--------------------------------------------------", "\n")
  cat("Results:\n")
  cat("Estimated Bandwidths: \t h_x:\t", object$bndw[1], "\n")
  cat("\t \t h_t:\t", object$bndw[2], "\n")
  cat("Variance Factor: \t c_f:\t", object$varCoef, "\n")
  cat("Iterations: ")
  
  
  cat("Used Parameters:\n")
  cat("Local Polynomial Order:\t", object$dcsOptions$pOrder, "\n")
  cat("Kernel Type:\t\t", nameKernFcn, "\n")
  cat("Variance Estimation:\t", object$dcsOptions$varEst, "\n")
  cat("--------------------------------", "\n")
  cat("Estimated Bandwidths:\n")
  cat("\t h_x:\t", object$bndw[1], "\n")
  cat("\t h_t:\t", object$bndw[2], "\n")
  cat("Variance Factor: \n")
  cat("\t c_f:\t", object$varCoef, "\n")
  cat("--------------------------------")  


print.dcs = function(object)
{
  nameKernFcn = paste0("MW", object$dcsOptions$kernPar[1],
                       object$dcsOptions$kernPar[2], "0")
  cat(class(object), "\n")
  cat("--------------------------------", "\n")
  cat("Used Parameters:\n")
  cat("Local Polynomial Order:\t", object$dcsOptions$pOrder, "\n")
  cat("Kernel Type:\t\t", nameKernFcn, "\n")
  cat("Variance Estimation:\t", object$dcsOptions$varEst, "\n")
  cat("--------------------------------", "\n")
  cat("Estimated Bandwidths:\n")
  cat("\t h_x:\t", object$bndw[1], "\n")
  cat("\t h_t:\t", object$bndw[2], "\n")
  cat("Variance Factor: \n")
  cat("\t c_f:\t", object$varCoef, "\n")
  cat("--------------------------------")
}
