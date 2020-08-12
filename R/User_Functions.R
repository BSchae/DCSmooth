###############################################################################
#                                                                             #
#                     DCSmooth Package: User Functions                        #
#                                                                             #
###############################################################################


#-----------------------Set Options via Function------------------------------#

setOptions = function(...) # user friendly wrapper function for .setOptions
{
  .setOptions(...)
}

#--------------Function for smoothing and bandwidth estimation----------------#

DCSmooth = function(Y, X = 1, T = 1, bndw = "auto", dcsOptions = setOptions())
{
  # set up vectors for X and T
  if (length(X) == 1 && length(T) == 1)
  {
    X = seq(from = 0, to = 1, length.out = dim(Y)[1])
    T = seq(from = 0, to = 1, length.out = dim(Y)[2])
  }
  
  # set kernel function
  nameKernFcn = paste0("MW", dcsOptions$kernPar[1], dcsOptions$kernPar[2], "0")
  kernelFcn  = kernelFcn_assign(nameKernFcn) # set kernel Function to use in optimization

  # bandwidth selection
  if (bndw == "auto" && dcsOptions$pOrder == 0)
  {
    bndwObj = KR_bndwSelect(Y, kernelFcn, dcsOptions)
    bndw = pmin(bndwObj$bndw, c(0.45, 0.45))
  }
  else if (bndw == "auto" && dcsOptions$pOrder > 0)
  {
    bndwObj = LP_bndwSelect(Y, kernelFcn, dcsOptions)
    bndw = pmin(bndwObj$bndw, c(0.45, 0.45))
  }
  
  if (dcsOptions$pOrder == 0)
  {
    DCSOut = KR_DoubleSmooth2(yMat = Y, hVec = bndw, drvVec = c(0, 0), 
      kernFcnPtrX = kernelFcn, kernFcnPtrT = kernelFcn)
  }
  else if (dcsOptions$pOrder > 0)
  {
    DCSOut = LP_DoubleSmooth(yMat = Y, hVec = bndw, polyOrderVec 
      = c(dcsOptions$pOrder, dcsOptions$pOrder), drvVec = c(0,0))
  }
  
  return(list(X = X, T = T, Y = Y, M = DCSOut, bndw = bndw,
              iterations = bndwObj$iterations, dcsOptions = dcsOptions))
}

#-------------------Function for plotting smoothed surface--------------------#

plotDCS = function(DCSobj, ...)
{
  .plotDCS(DCSobj = DCSobj, ...)
}
