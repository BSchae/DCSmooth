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
  kerns     = c("MW_220", "MW_220"), # choose a kernel function
  drv       = c(0, 0),
  var_est   = "iid",
  
  # advanced options in ellipsis
  ...
)
{
  # get ellipsis
  IPI_options = list(...)
  
  # check if inputs are vectors
  if (length(kerns) == 1) { kerns = c(kerns, kerns) }
  if (length(drv) == 1) { drv = c(drv, drv) }
  
  exception.check.options.input(type, kerns, drv, var_est, IPI_options)
  
  if (var_est == "lm" && type == "KR")
  {
    warning("Long-Memory bandwidth selection only supported for local ", 
            "polynomial regression. Type set to \"LP\".")
    type = "LP"
  }
  if (var_est == "lm")
  {
    message("Estimation under long-memory errors (SFARIMA) is currently in ", 
            "experimantal state.")
  }
  
  # Select options according to type ("LP", "KR")
  if (type == "LP")
  {
    p_order = drv + 1
    if (!exists("infl_exp", IPI_options))
    {
      IPI_options$infl_exp = c("auto", " ")
    }
    if (!exists("infl_par", IPI_options)) { IPI_options$infl_par = c(1, 1) }
    if (!exists("delta", IPI_options)) { IPI_options$delta = c(0.05, 0.05) }
    if (!exists("const_window", IPI_options))
    {
      IPI_options$const_window = FALSE
    }
    
    options_list = list(type = type, kerns = kerns, p_order = p_order,
                       drv = drv, var_est = var_est, IPI_options = IPI_options)
  } else if (type == "KR") {
    p_order = NA
    if (!exists("infl_exp", IPI_options)) { IPI_options$infl_exp = c(0.5, 0.5) }
    if (!exists("infl_par", IPI_options)) { IPI_options$infl_par = c(2, 1) }
    if (!exists("delta", IPI_options)) { IPI_options$delta = c(0.05, 0.05) }
    if (!exists("const_window", IPI_options))
    {
      IPI_options$const_window = FALSE
    }
    
    options_list = list(type = type, kerns = kerns, p_order = p_order,
                        drv = drv, var_est = var_est, IPI_options = IPI_options)
  } else {
    stop("Unknown type \"", type, "\"")
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
    # axx = list(titlefont = f1, tickfont = f2, title = t_lab,
    #           gridcolor = "lightgray", zerolinecolor = "black")
    # axy = list(titlefont = f1, tickfont = f2, title = x_lab,
    #           gridcolor = "lightgray", zerolinecolor = "black",
    #           autorange = "reversed")
    # axz = list(titlefont = f1, tickfont = f2, title = y_lab,
    #           backgroundcolor = "white",  gridcolor = "lightgray",
    #           showbackground = TRUE, zerolinecolor = "black",
    #           tickformat = ".1e")
    axx = list(title = t_lab, gridcolor = "lightgray", zerolinecolor = "black",
               automargin = TRUE)
    axy = list(title = x_lab, gridcolor = "lightgray", zerolinecolor = "black",
               autorange = "reversed", automargin = TRUE)
    axz = list(title = y_lab, backgroundcolor = "white",  gridcolor = "lightgray",
               showbackground = TRUE, zerolinecolor = "black",
               tickformat = ".1e", automargin = TRUE) #, tickmode = "linear", dtick = 4000, tick0 = 2000)
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
                           projection = list(type = "orthographic")),
             autosize = FALSE)
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
