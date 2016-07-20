#################################################################################
#
#  -file = "min_res_bootstrapped.R"   Code written March 2014
#  - Author: Kelsey Ruckert (klr324@psu.edu)
#
#  -This function find the sum of the absolute residuals estimated from the model
#       function in "Deoptim_rahm_model.R". Finding the sum is used with the DEoptim
#       R function to minimize the residuals and find the best values for the
#       parameters. The best parameter values are used as initial values for the
#       bootstrap and MCMC codes as described in Ruckert et al. (GRL 2016).
#
# THIS CODE IS PROVIDED AS-IS WITH NO WARRANTY (NEITHER EXPLICIT
# NOT IMPLICIT).  I SHARE THIS CODE IN HOPES THAT IT IS USEFUL,
# BUT I AM NOT LIABLE FOR THE BEHAVIOR OF THIS CODE IN YOUR OWN
# APPLICATION.  YOU ARE FREE TO SHARE THIS CODE SO LONG AS THE
# AUTHOR(S) AND VERSION HISTORY REMAIN INTACT.
#
# To use this function, simply source this file:
#   source("min_res_bootstrapped.R")
#
# INPUTS:
#   vector of parameter values:
#       p[1]: alpha; sensitivity of sea-level to temperature changes
#       p[2]: T_0; temperature when sea-level anomaly is zero
#       p[3]: H_0; initial sea-level anomaly
#
##################################################################################

min_res = function(p){
    # Return the sum of the absolute values
    # from sea-level minus bootstrap simulated model values
    sum(abs(slr.boot.sim[i,] - model(p)))
}

################################## END ###########################################