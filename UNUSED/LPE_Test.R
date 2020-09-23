# ###############################################################################
# #                                                                             #
# #         DCSmooth Package: R-Functions for LP Bandwidth Selection            #
# #                                                                             #
# ###############################################################################
# 
# LPSmooth = function(X, Y, p, nu, h) {  #
#   n = length(X)     # 
#   ord = order(X)    # 
#   X = X[ord]        #  
#   Y = Y[ord]        # 
#   Ye.LP = (1:n)*0   # 
#   
#   for(i in 1:n) {
#     X0 = X[i]                  # 
#     Xi = X[abs(X-X0)<h] - X0   # 
#     Yi = Y[abs(X-X0)<h]        # 
#     u  = Xi/h                  # 
#     WK = kernFkt_MW220(u, 1)   # 
#     
#     # 
#     if(p == 0) { Ye.LP[i] = sum(Yi*WK)/sum(WK) }
#     if(p == 1) { Ye.LP[i] = (lm(Yi~Xi, weights=WK)$coefficients[nu+1])*factorial(nu) }
#     if(p == 2) { Ye.LP[i] = (lm(Yi~Xi+I(Xi^2), weights=WK)$coefficients[nu+1])*factorial(nu) }
#     if(p == 3) { Ye.LP[i] = (lm(Yi~Xi+I(Xi^2)+I(Xi^3), weights=WK)$coefficients[nu+1])*factorial(nu) }
#     if(p == 4) { Ye.LP[i] = (lm(Yi~Xi+I(Xi^2)+I(Xi^3)+I(Xi^4), weights=WK)$coefficients[nu+1])*factorial(nu) }
#   }
#   
#   # 
#   results=list(X=X, Y=Y, Ye.LP=Ye.LP)
#   #   
#   return(results)
# }
# 
# LPE_test = function(Y, drvVec, hVec)
# {
#   nX = dim(Y)[1]
#   nT = dim(Y)[2]
#   
#   #hVec = hVec/3
#   
#   M = M2 = Y*NA
#   
#   T = 0:(nT - 1)/(nT - 1)
#   X = 0:(nX - 1)/(nX - 1)
#   
#   for (i in 1:nX)
#   {
#     M[i, ] = LPSmooth(T, Y[i, ], p = drvVec[2] + 2, nu = drvVec[2], h = hVec[2])$Ye.LP
#   }
#   
#   for (j in 1:nT)
#   {
#     M2[, j] = LPSmooth(X, t(M[, j]), p = drvVec[1] + 2, nu = drvVec[1], h = hVec[1])$Ye.LP
#   }
#   
#   return(M2)
# }
#   
# 
# LP_bndwSelect2 = function(Y, kernelFcn, dcsOptions)
# {
#   nX = dim(Y)[1]; nT = dim(Y)[2]
#   n  = nX * nT                                  # total number of observations is needed later
#   
#   kernelProp = kernelPropFcn(kernelFcn)         # calculate properties R and mu_2 of kernel
#   
#   hOpt = c(0.1, 0.1)                            # initial values for h_0, arbitrary chosen
#   
#   iterate = TRUE                                # iteration indicator
#   iterationCount = 0
#   while(iterate)                                # loop for IPI
#   {
#     iterationCount = iterationCount + 1
#     hOptTemp   = hOpt[1:2]       # store old bandwidths for breaking condition
#     hInfl      = inflationFcnLP(hOptTemp, c(nX, nT), dcsOptions)  # inflation of bandwidths for drv estimation
#     
#     if (dcsOptions$constWindow == TRUE)
#     {
#       
#     } else {
#       # pre-smoothing of the surface function m(0,0) for better estimation of derivatives
#       YSmth = LPE_test(Y, c(0, 0), hOptTemp)
#       
#       # smoothing of derivatives m(2,0) and m(0,2)
#       mxx = LPE_test(Y, c(2, 0), hInfl$h_xx)
#       mtt = LPE_test(Y, c(0, 2), hInfl$h_tt)
#     }
#     
#     # shrink mxx, mtt from boundaries
#     if (dcsOptions$delta[1] != 0 || dcsOptions$delta[2] != 0)
#     {
#       shrinkX = ceiling(dcsOptions$delta[1]*nX):
#         floor((1 - dcsOptions$delta[1])*nX)
#       shrinkT = ceiling(dcsOptions$delta[2]*nT):
#         floor((1 - dcsOptions$delta[2])*nT)
#       
#       mxx = mxx[shrinkX, shrinkT]
#       mtt = mtt[shrinkX, shrinkT]
#       nSub = dim(mxx)[1]*dim(mxx)[2]
#     } else {
#       nSub = n
#     }
#     
#     # calculate variance factor
#     if (dcsOptions$varEst == "iid")
#     {
#       varCoef = (sd(Y - YSmth))^2
#     }
#     else if (dcsOptions$varEst == "qarma")
#     {
#       if (exists("model", dcsOptions))
#       {
#         model = dcsOptions$model
#       } else {
#         model = list(ar_x = 1, ar_t = 1, ma_x = 1, ma_t = 1)
#       }
#       
#       varCoef = QARMA.cf(Y - YSmth, model = model)
#     }
#     
#     
#     # calculate optimal bandwidths for next step
#     hOpt = hOptLP(mxx, mtt, varCoef, n, nSub, dcsOptions$pOrder, kernelProp)
#     
#     # break condition
#     if( ((hOpt[1]/hOptTemp[1] - 1 < 0.001) && (hOpt[2]/hOptTemp[2] - 1 
#       < 0.001) && (iterationCount > 2)) || iterationCount > 10)
#     {
#       iterate = FALSE
#     }
#   }
#   return(list(bndw = hOpt, iterations = iterationCount))
# }