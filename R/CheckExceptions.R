################################################################################
#                                                                              #
#                 DCSmooth Package: Checks for Exceptions                      #
#                                                                              #
################################################################################

# These functions check the data, the options and a given bandwidth for conform-
# ity with the requirements of the package. They return TRUE, if everything is
# okay.

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
  
  return(TRUE)
}

#-----------------------Check for correct Options------------------------------#

# This function is used in the .setOptions() function for direct check as well
# as in the dcs() function.

.dcsCheck_Options = function(dcsOpt)
{
  # check kernel parameters
  
  # check polynomial order
  if (dcsOpt$pOrder < 0) stop("Your polynomial order is smaller than 0.")
  if (dcsOpt$pOrder > 3) warning("Your large polynomial order might result
                                 in a longer computation time.")
  
  # check inflation exponents
  if (any(dcsOpt$inflExp != 0.5)) warning("Inflation exponents have been
                                          changed.")
  # check inflation parameters
  if (dcsOpt$inflPar != c(2, 1)) warning("Inflation parameters have been
                                          changed.")
  
  # check for const window
  if (dcsOpt$inflPar == TRUE && dcsOpt$pOrder > 0) warning("Constant window
                        width is not yet supported for polynomial regression")
  
  # check variance estimation method
  if (!(dcsOpt$varEst %in% c("iid", "qarma"))) stop("unsupported method in
                                    varEst (use only \"iid\" and \"qarma\")")
}

set


