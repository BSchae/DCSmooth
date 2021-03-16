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

#----------------------Check for useable bandwidths----------------------------#

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
  # check kernel parameters
  
  # check derivative orders
  if (any(dcsOpt$pOrder < 0))
  {  
    stop("Your derivative order is smaller than 0.")
  }
  
  # check inflation exponents
  if (any(dcsOpt$inflExp != 0.5))
  {
    warning("Inflation exponents have been changed.")
  }
  
  # check inflation parameters
  if (any(dcsOpt$inflPar != c(1.5, 0.25)))
  {
    warning("Inflation parameters have been changed.")
  }
  
  # check for const window
  if (dcsOpt$inflPar == TRUE && dcsOpt$pOrder > 0)
  {
    warning("Constant window width is not yet supported for polynomial",
              "regression")
  }
  
  # check variance estimation method
  if (!(dcsOpt$varEst %in% c("iid", "qarma")))
  {
    stop("unsupported method in varEst (use only \"iid\" and \"qarma\")")
  }
}