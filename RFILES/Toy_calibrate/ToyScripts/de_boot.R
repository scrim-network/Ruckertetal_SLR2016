###################################################################################
#
#  -file = "de_boot.R"   Code written March 2014
#  - Author: Kelsey Ruckert (klr324@psu.edu)
#
#  -This function sources a simple linear model to be sourced in the DEoptim
#       R function. Deoptim is used in the bootstrap and MCMC codes as described
#       in Ruckert et al. (2016) to find good initial values for each of the
#       parameters.
#
# THIS CODE IS PROVIDED AS-IS WITH NO WARRANTY (NEITHER EXPLICIT
# NOT IMPLICIT).  I SHARE THIS CODE IN HOPES THAT IT IS USEFUL, 
# BUT I AM NOT LIABLE FOR THE BEHAVIOR OF THIS CODE IN YOUR OWN
# APPLICATION.  YOU ARE FREE TO SHARE THIS CODE SO LONG AS THE
# AUTHOR(S) AND VERSION HISTORY REMAIN INTACT.
#
# To use this function, simply source this file:
#   source("de_boot.R")
#
###################################################################################

model = function(p){ # p represents the parameters in a vector
  
  obser=p[1]+p[2]*x # estimate observations using a linear model
  #p[1] = a
  #p[2] = b
  return(obser)
}

#################################### END ##########################################
