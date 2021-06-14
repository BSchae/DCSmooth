################################################################################
#                                                                              #
#                 DCSmooth Package: Checks for Exceptions                      #
#                                                                              #
################################################################################

# These functions check the data, the options and a given bandwidth for conform-
# ity with the requirements of the package. They stop the function and return an
# error, if something wrong.

#----------------------Check Matrix of Observations Y--------------------------#

exception.check.Y = function(Y, X, T)
{
  sys.call(-1)
  # check for missing values
  if (any(is.na(Y)))
  {
    stop("Y contains missing values (NAs)")
  }
  if (any(is.na(X)) || any(is.na(T)))
  {
    stop("X and/or T contain missing values (NAs)")
  }
  
  # check if Y, X, T is numeric
  if (!is.numeric(Y) || !is.matrix(Y))
  {
    stop("Y must be a numeric matrix")
  }
  if (!is.numeric(X) || !is.numeric(T))
  {
    stop("X and/or T be a numeric vector")
  }
  
  # check for correct dimension of Y
  if (any(dim(Y) < 3)) {
    stop("Y has to be at least of dimension 3 in each direction.")
  }
  
  # check for correct dimension of X, T
  if (dim(Y)[1] != length(X))
  {
    stop("Rows of Y not matching X")
  }
  if (dim(Y)[2] != length(T))
  {
    stop("Columns of Y not matching T")
  }
}

#-----------------------Check for usable bandwidths----------------------------#

exception.check.bndw = function(bndw, dcs_options)
{
  if (bndw[1] != "auto")
  {
    if (!(is.numeric(bndw)) || (length(bndw) != 2))
    {
      stop("bndw must be a numeric vector of length 2.")
    }
      
    if (any(bndw < 0))
    {
      stop("bndw must be positive")
    }  
      
    if (any(bndw > 0.45) && dcs_options$type == "KR")
    {
      stop("bndw must be < 0.45 for kernel regression")
    }  
      
    if (any(bndw > 0.5))
    {
      warning("bandwidth seems unusually high, computation time might be",
             "increased.")
    }  
  }
}

#-----------------------Check for correct Options------------------------------#

# This function is used in the .setOptions() function for direct check as well
# as in the dcs() function.

exception.check.options = function(dcs_opt)
{
  sys.call(-1)
  # check class
  if(!(class(dcs_opt) == "dcs_options"))
  {
    stop("Incorrect options specified, please use \"setOptions()\".")
  }
  
  # check for unknown or missing options
  unknown_name = names(dcs_opt)[which(!(names(dcs_opt) %in% dcs_list_options))]
  if (length(unknown_name) > 0)
  {
    warning("Option \"", unknown_name, "\" is unknown and will be ignored.")
  }
  # unspec_name = dcs_list_options[which(!(dcs_list_options %in%
  #                                                 names(dcs_opt)))]
  # if (length(unspec_name) > 0)
  # {
  #   stop("Option \"", unknown_name, "\" not specified.")
  # }
  
  # check kernels
  if (!(dcs_opt$kerns[1] %in% DCSmooth:::dcs_list_kernels) || 
      !(dcs_opt$kerns[1] %in% DCSmooth:::dcs_list_kernels))
  {
    stop("Unsuppored kernels specified.")
  }
  
  # check regression type
  if (!(dcs_opt$type %in% c("LP", "KR")))
  {
    stop("Unsupported regression type. Choose \"KR\" or \"LP\"")
  }
  
  ### Options for Local Polynomial Regression
  if (dcs_opt$type == "LP")
  {
    # check derivative orders
    if (any(dcs_opt$drv < 0))
    {  
      stop("Your derivative order is smaller than 0.")
    }
    
    # check inflation exponents
    if (any(dcs_opt$infl_exp != 0.6))
    {
      warning("Inflation exponents have been changed.")
    }
  
    # check inflation parameters
    if (any(dcs_opt$infl_par != c(1, 1)))
    {
      warning("Inflation parameters have been changed.")
    } 
  }
  
  ### Options for Kernel Regression
  if (dcs_opt$type == "KR")
  {
    # 
    # check inflation exponents
    if (any(dcs_opt$infl_exp != 0.5))
    {
      warning("Inflation exponents have been changed.")
    }
    
    # check inflation parameters
    if (any(dcs_opt$infl_par != c(2, 1)))
    {
      warning("Inflation parameters have been changed.")
    } 
  }
  
  ### Options for variance estimation method
  if (!(dcs_opt$var_est %in% c("iid", "qarma", "qarma_gpac", "qarma_bic")))
  {
    stop("unsupported method in varEst (use only \"iid\" and \"qarma\")")
  } else {
    if (dcs_opt$var_est == "iid")
    {
      if (exists("qarma_order", where = dcs_opt))
      {
        warning("QARMA order is no valid parameter in iid. case and will be
                 ignored.")
      }
    } else if (dcs_opt$var_est == "qarma") {
      # check for correct qarma_order
      if (!exists("qarma_order", where = dcs_opt))
      {
        stop("No model order for QARMA estimation provided.")
      } else {
        if (is.list(dcs_opt$qarma_order) &&
            exists("ar", where = dcs_opt$qarma_order) && 
            exists("ma", where = dcs_opt$qarma_order))
        {
          if (!(all(dcs_opt$qarma_order$ar %in% 0:100) &&
                length(dcs_opt$qarma_order$ar) == 2))
          {
            stop("unsupported values in AR-order")
          }
          if (!(all(dcs_opt$qarma_order$ma %in% 0:100) &&
                length(dcs_opt$qarma_order$ma) == 2))
          {
            stop("unsupported values in MA-order")
          }
        } else if (!any(dcs_opt$qarma_order %in% c("gpac", "bic"))) {
          stop("unsupported values in model order")
        }
        
        if (any(dcs_opt$qarma_order %in% c("gpac", "bic")) && 
            exists("order_max", where = dcs_opt))
        {
          if (!(is.list(dcs_opt$order_max) &&
              exists("ar", where = dcs_opt$order_max) &&
              exists("ma", where = dcs_opt$order_max)))
          {
            stop("Unsupported values for max. order of order selection.")
          } else {
            if (!(all(dcs_opt$order_max$ar %in% 0:100) &&
                  length(dcs_opt$order_max$ar) == 2))
            {
              stop("unsupported values in AR-part of max. order")
            }
            if (!(all(dcs_opt$order_max$ma %in% 0:100) && 
                  length(dcs_opt$order_max$ma) == 2))
            {
              stop("unsupported values in MA-part of max. order")
            }
          }
        }
      }
    }
  }
}