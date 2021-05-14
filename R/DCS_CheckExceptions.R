################################################################################
#                                                                              #
#                 DCSmooth Package: Checks for Exceptions                      #
#                                                                              #
################################################################################

# These functions check the data, the options and a given bandwidth for conform-
# ity with the requirements of the package. They stop the function and return an
# error, if something wrong.

#----------------------Check Matrix of Observations Y--------------------------#

.dcsCheck_Y = function(Y)
{
  # check for missing values
  if (any(is.na(Y)))
  {
    stop("Y contains missing values (NAs)")
  }
  
  # check if Y is numeric
  if (!is.numeric(Y))
  {
    stop("Y must be a numeric matrix")
  }
  
  # check for correct dimension of Y
  if (any(dim(Y) < 3)) {
    stop("Y has to be at least of dimension 3 in each direction.")
  }
}

#-----------------------Check for usable bandwidths----------------------------#

.dcsCheck_bndw = function(bndw, dcsOpt)
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
      
    if (any(bndw > 0.45) && dcsOpt$pOrder == 0)
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

.dcsCheck_options = function(dcsOpt)
{
  if(!(class(dcsOpt) == "dcsOpt"))
  {
    stop("Incorrect options specified, please use \"setOptions()\".")
  }
  
  # check regression type
  if (!(dcsOpt$type %in% c("LP", "KR")))
  {
    stop("Unsupported regression type. Choose \"KR\" or \"LP\"")
  }
  
  ### Options for Local Polynomial Regression
  if (dcsOpt$type == "LP")
  {
    # check derivative orders
    if (any(dcsOpt$drv < 0))
    {  
      stop("Your derivative order is smaller than 0.")
    }
    
    # check inflation exponents
    if (any(dcsOpt$inflExp != 0.6))
    {
      warning("Inflation exponents have been changed.")
    }
  
    # check inflation parameters
    if (any(dcsOpt$inflPar != c(1.5, 1)))
    {
      warning("Inflation parameters have been changed.")
    } 
  }
  
  ### Options for Kernel Regression
  if (dcsOpt$type == "KR")
  {
    # 
    # check inflation exponents
    if (any(dcsOpt$inflExp != 0.6))
    {
      warning("Inflation exponents have been changed.")
    }
    
    # check inflation parameters
    if (any(dcsOpt$inflPar != c(3, 1)))
    {
      warning("Inflation parameters have been changed.")
    } 
  }
  
  ### Options for variance estimation method
  if (!(dcsOpt$varEst %in% c("iid", "qarma", "qarma_gpac", "qarma_bic")))
  {
    stop("unsupported method in varEst (use only \"iid\" and \"qarma\")")
  } else {
    if (dcsOpt$varEst == "iid")
    {
      if (exists("modelOrder", where = dcsOpt))
      {
        warning("model order is no valid parameter in iid. case and will be
                 ignored.")
      }
    } else if (dcsOpt$varEst == "qarma") {
      # check for correct modelOrder
      if (!exists("modelOrder", where = dcsOpt))
      {
        stop("No model order for QARMA estimation provided.")
      } else {
        if (is.list(dcsOpt$modelOrder) &&
            exists("ar", where = dcsOpt$modelOrder) && 
            exists("ma", where = dcsOpt$modelOrder))
        {
          if (!(all(dcsOpt$modelOrder$ar %in% 0:100) &&
                length(dcsOpt$modelOrder$ar) == 2))
          {
            stop("unsupported values in AR-order")
          }
          if (!(all(dcsOpt$modelOrder$ma %in% 0:100) &&
                length(dcsOpt$modelOrder$ma) == 2))
          {
            stop("unsupported values in MA-order")
          }
        } else if (!any(dcsOpt$modelOrder %in% c("gpac", "bic"))) {
          stop("unsupported values in model order")
        }
        
        if (any(dcsOpt$modelOrder %in% c("gpac", "bic")) && 
            exists("orderMax", where = dcsOpt))
        {
          if (!(is.list(dcsOpt$orderMax) &&
              exists("ar", where = dcsOpt$orderMax) &&
              exists("ma", where = dcsOpt$orderMax)))
          {
            stop("Unsupported values for max. order of order selection.")
          } else {
            if (!(all(dcsOpt$orderMax$ar %in% 0:100) &&
                  length(dcsOpt$orderMax$ar) == 2))
            {
              stop("unsupported values in AR-part of max. order")
            }
            if (!(all(dcsOpt$orderMax$ma %in% 0:100) && 
                  length(dcsOpt$orderMax$ma) == 2))
            {
              stop("unsupported values in MA-part of max. order")
            }
          }
        }
      }
    }
  }
}