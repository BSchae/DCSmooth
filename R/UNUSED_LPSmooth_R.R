# ####### R function for LP estimation of m & its derivatives ###
# 
# # put a comment behind every '#'
# 
# LPSmooth_R = function(X, Y, p, nu, mu.K, h) {  # start "function editor", give function a name and declare parameters
#   n = length(X)     # Number of observations  
#   ord = order(X)    # Sort observations in x in ascending order
#   X = X[ord]        # Sort X 
#   Y = Y[ord]        # Sort Y in order of X
#   Ye.LP = (1:n)*0   # Vector to store smoothed data
#   
#   for(i in 1:n) {
#     X0 = X[i]                  # The X point at which to estimate
#     Xi = X[abs(X-X0)<h] - X0   # Taking the localized and centralized values of X
#     Yi = Y[abs(X-X0)<h]        # Taking the localized values of Y
#     u  = Xi/h                  # Taking u for the kernel function
#     WK = (1 - u^2)^mu.K        # Calculating the weights (standardization not necessary)
#     
#     # Local polynomial estimation based on lm, for polynomial order p
#     if(p == 0) { Ye.LP[i] = sum(Yi*WK)/sum(WK) }
#     if(p == 1) { Ye.LP[i] = (lm(Yi~Xi, weights=WK)$coefficients[nu+1])*factorial(nu) }
#     if(p == 2) { Ye.LP[i] = (lm(Yi~Xi+I(Xi^2), weights=WK)$coefficients[nu+1])*factorial(nu) }
#     if(p == 3) { Ye.LP[i] = (lm(Yi~Xi+I(Xi^2)+I(Xi^3), weights=WK)$coefficients[nu+1])*factorial(nu) }
#     if(p == 4) { Ye.LP[i] = (lm(Yi~Xi+I(Xi^2)+I(Xi^3)+I(Xi^4), weights=WK)$coefficients[nu+1])*factorial(nu) }
#   }
#   
#   # All outputs -- X, Y ordered according to X
#   results=list(X=X, Y=Y, Ye.LP=Ye.LP)
#   # output the results  
#   return(results)
# }