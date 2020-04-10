.setOptions = function(
  bandwidth = "auto",
  kernPar   = c(2, 2),
  pOrder    = 1
)
{
  return(list(bandwidth = bandwidth, kernPar = kernPar, pOrder = pOrder))
}

setOptions = function(...)
{
  .setOptions(...)
}

bndwSelect = function(Y, dcsOptions = setOptions())
{
  n = dim(Y)[1] * dim(Y)[2]
  
  hOpt = c(0.1, 0.1)
  iterate = TRUE
  s = 0
  while(iterate)
  {
    s = s + 1
    hOptTemp = hOpt
    mxx = FastDoubleSmooth(yMat = Y, hVec = hOptTemp, polyOrderVec 
                      = c(dcsOptions$pOrder + 2, dcsOptions$pOrder), drvVec = c(2, 0))
    mtt = FastDoubleSmooth(yMat = Y, hVec = hOptTemp, polyOrderVec 
      = c(dcsOptions$pOrder, dcsOptions$pOrder + 2), drvVec = c(0, 2))
    
    hOpt = hOptFunction(mxx, mtt, 1, n, dcsOptions$pOrder, list(R = 0.5, mu = 0.5))
    if(s > 5) { iterate = FALSE }
  }
  
  return(hOpt)
}
  
hOptFunction = function(mxx, mtt, varCoef, n, p, kernelProp)
{
  i0x = intCalc(mxx, mtt, n, p)
  i0t = intCalc(mtt, mxx, n, p)
  
  hxOpt = (kernelProp$R^2 * varCoef)/((p + 1) * n * kernelProp$mu^2 * i0x)
  htOpt = (kernelProp$R^2 * varCoef)/((p + 1) * n * kernelProp$mu^2 * i0x)
  hxOpt = hxOpt^(1/(2*p + 4))
  htOpt = htOpt^(1/(2*p + 4))
  
  return(c(hxOpt, htOpt))
}

intCalc = function(m11, m22, n, p)
{
  i11 = sum(m11 * m11)/n
  i22 = sum(m22 * m22)/n
  i12 = sum(m11 * m22)/n
  
  iOut = (i11/i22)^(1/(2*p + 2)) * (2*i11 + sqrt(i11/i22) * i12)
  
  return(iOut)
}
