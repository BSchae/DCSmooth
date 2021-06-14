################################################################################
#                                                                              #
#                DCSmooth Package: General Auxiliary Funcions                  #
#                                                                              #
################################################################################

#------------------------Set Options via Function------------------------------#

.set.options = function(    # inside function with default values in arguments
  # standard options
  type      = "LP",         # either "LP" for local polynomial regression or
                            # "KR" for kernel regression
  drv       = c(0, 0),
  var_est   = "iid",
  qarma_order = "auto",
  order_max = "auto",
  
  # advanced options
  kerns     = c("MW220", "MW220"), # choose a kernel function
  #bound_mod = FALSE,
  infl_exp  = "auto",       # inflation exponent
  infl_par  = "auto",       # inflation parameters c (regression),
                            # d (2nd derivative)
  delta     = "auto",       # parameter for shrinking the derivatives
  const_window = FALSE     # should a constant window be maintained?
  #bound_imprv  = FALSE      # should estimation at boundaries be improved?
  
  # Options with default value "auto" are dependent on type or var_est and will
  # be selected in the following, if not given a specific value.
)
{
  # Select options according to type ("LP", "KR")
  if (type == "LP")
  {
    p_order = drv + 1
    if (infl_exp[1] == "auto") { infl_exp = c(0.6, 0.6) }
    if (infl_par[1] == "auto") { infl_par = c(1, 1) }
    if (delta[1] == "auto")    { delta = c(0.05, 0.05) }
    
    options_list = list(type = type, kerns = kerns, p_order = p_order,
                       drv = drv, infl_exp = infl_exp, infl_par = infl_par,
                       delta = delta, const_window = const_window,
                       #bound_imprv = bound_imprv,
                       var_est = var_est)
  } else if (type == "KR") {
    if (infl_exp[1] == "auto") { infl_exp = c(0.5, 0.5) }
    if (infl_par[1] == "auto") { infl_par = c(2, 1) }
    if (delta[1] == "auto")   { delta = c(0.05, 0.05) }
    options_list = list(type = type, kerns = kerns, drv = drv, 
                        infl_exp = infl_exp, infl_par = infl_par,
                        delta = delta, const_window = const_window,
                        #bound_imprv = bound_imprv,
                        var_est = var_est)
  }
  
  # Select options according to var_est ("iid", "QARMA")
  if (var_est == "qarma")
  {
    if (qarma_order[1] == "auto")
    {
      qarma_order = list(ar = c(1, 1), ma = c(1, 1))
    }
    if (any(qarma_order %in% c("gpac", "bic")))
    {
      if (order_max == "auto")
      {
        order_max = list(ar = c(1, 1), ma = c(1, 1))
      }
      options_list = c(options_list, list(order_max = order_max))
    }
    
    options_list = c(options_list, list(qarma_order = qarma_order))
  }
  
  # apply class to output object
  class(options_list) = "dcs_options"
  
  exception.check.options(options_list)
  
  return(options_list)
}

#----------------------Function for 3d plots using plotly----------------------#

.plotly.3d = function(Y, X = 1, T = 1,
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

plot.dcs.colors = function(Y, color, n_color = 100)
{
  Y_color = Y - min(Y)
  Y_color = Y_color/max(Y_color) * (n_color - 1) + 1
  color_fcn = grDevices::colorRampPalette(colors = color)
  color_grad = color_fcn(n_color)
  color_out = matrix(color_grad[Y_color])
  return(color_out)
}
