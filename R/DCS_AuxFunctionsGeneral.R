################################################################################
#                                                                              #
#                DCSmooth Package: General Auxiliary Funcions                  #
#                                                                              #
################################################################################

#------------------------Set Options via Function------------------------------#

.setOptions = function(    # inside function with default values in arguments
  type      = "LP",        # either "LP" for local polynomial regression or
                           # "KR" for kernel regression
  kernPar   = c(2, 2),     # choose a kernel function with mu = 2
  drv       = c(0, 0),
  inflExp   = "auto",      # inflation exponent
  inflPar   = "auto",      # inflation parameters c (regression),
                           # d (2nd derivative)
  delta     = "auto",      # parameter for shrinking the derivatives
  constWindow = FALSE,
  varEst    = "iid",
  qarmaOrder = "auto",
  orderMax = "auto"
  
  # Options with default value "auto" are dependent on type or varEst and will
  # be selected in the following, if not given a specific value.
)
{
  # Select options according to type ("LP", "KR")
  if (type == "LP")
  {
    pOrder = drv + 1
    if (inflExp[1] == "auto") { inflExp = c(0.6, 0.6) }
    if (inflPar[1] == "auto") { inflPar = c(1.5, 1) }
    if (delta[1] == "auto")   { delta = c(0, 0) }
    
    options_list = list(type = type, kernPar = kernPar, pOrder = pOrder,
                       drv = drv, inflExp = inflExp, inflPar = inflPar,
                       delta = delta, constWindow = constWindow,
                       varEst = varEst)
  } else if (type == "KR") {
    if (inflExp[1] == "auto") { inflExp = c(0.6, 0.6) }
    if (inflPar[1] == "auto") { inflPar = c(3, 1) }
    if (delta[1] == "auto")   { delta = c(0.05, 0.05) }
    options_list = list(type = type, kernPar = kernPar, drv = drv, 
                        inflExp = inflExp, inflPar = inflPar,
                        delta = delta, constWindow = constWindow,
                        varEst = varEst)
  }
  
  # Select options according to varEst ("iid", "QARMA")
  if (varEst == "qarma")
  {
    if (modelOrder[1] == "auto") { modelOrder = list(ar = c(1, 1), ma = c(1, 1)) }
    if (any(modelOrder %in% c("gpac", "bic")))
    {
      if (orderMax == "auto") { orderMax = list(ar = c(1, 1), ma = c(1, 1)) }
      options_list = c(options_list, list(orderMax = orderMax))
    }
    
    options_list = c(options_list, list(modelOrder = modelOrder))
  }
  
  # apply class to output object
  class(options_list) = "dcsOpt"
  
  .dcsCheck_options(options_list)
  
  return(options_list)
}

#----------------------Function for 3d plots using plotly----------------------#

.plotly3d = function(Y, X = 1, T = 1,
                     color = c("#444C5C", "#78A5A3", "#E1B16A", "#CE5A57"),
                     x_lab = "X_1", t_lab = "X_2", y_lab = "Value", 
                     showaxes = TRUE)
{
  if (class(Y)[1] == "dcs")
  {
    Y_data = Y$M
    X0 = Y$X
    T0 = Y$T
  } else {
    Y_data = Y
    
    if(length(X) == 1)
    {
      X0 = seq(from = 0, to = 1, length.out = dim(Y_data)[1])
    } else {
      X0 = X
    }
    if (length(T) == 1)
    {
      T0 = seq(from = 0, to = 1, length.out = dim(Y_data)[2])
    } else {
      T0 = T
    }
  }
  
  f1 = list(family = "Computer Modern", size = 14, color = "black")
  f2 = list(family = "Courier New", size = 14, color = "black")
  
  if (showaxes == TRUE)
  {
    axx = list(titlefont = f1, tickfont = f2, title = t_lab,
              gridcolor = "lightgray", zerolinecolor = "black")
    axy = list(titlefont = f1, tickfont = f2, title = x_lab,
              gridcolor = "lightgray", zerolinecolor = "black",
              autorange = "reversed")
    axz = list(titlefont = f1, tickfont = f2, title = y_lab,
              backgroundcolor = "white",  gridcolor = "lightgray",
              showbackground = TRUE, zerolinecolor = "black",
              tickformat = ".1e")
  } else {
    axx = list(title = "", showticklabels = FALSE, gridcolor = "lightgray",
               zerolinecolor = "black")
    axy = list(title = "", showticklabels = FALSE, gridcolor = "lightgray", 
               zerolinecolor = "black", autorange = "reversed")
    axz = list(title = "", showticklabels = FALSE, backgroundcolor = "white",
               gridcolor = "lightgray", showbackground = TRUE, 
               zerolinecolor = "black")
  }
  scene = list(xaxis = axx, yaxis = axy, zaxis = axz,
             camera = list(eye = list(x = -1.8, y = 1.2, z = 1),
                           projection = list(type = "orthographic")))
  # projection can be "projection" or "orthographic"
  
  fig = plotly::plot_ly(x = T0, y = X0, z = Y_data, showscale = FALSE,
              width = 600, height = 600)
  fig = plotly::add_surface(fig, colors = color)
  fig = plotly::layout(fig, scene = scene, margin = list(l = 100))
  
  return(fig)
}