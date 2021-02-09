LP_Smooth_test = function(yMat, hVec, polyOrderVec,
                          drvVec, kernelFcn)
{
  nX = dim(yMat)[1]; nT = dim(yMat)[2]

  if (hVec[1] > 0.45) { hVec[1] = 0.45 }
  if (hVec[2] > 0.45) { hVec[2] = 0.45 }
  
  m1 = matrix(NA, nrow = nX, ncol = nT)
  m2 = matrix(NA, nrow = nX, ncol = nT)

  for (i in 1:nX)
  {
    m1[i, ] = gsmoothCalcCpp(yMat[i, ], v = drvVec[2], p = polyOrderVec[2],
                             mu = 2, b = hVec[2], bb = 1)
  }
  for (j in 1:nT)
  {
    m2[, j] = gsmoothCalcCpp(m1[, j], v = drvVec[1], p = polyOrderVec[1],
                             mu = 2, b = hVec[1], bb = 1)
  }
  
  return(m2)
}