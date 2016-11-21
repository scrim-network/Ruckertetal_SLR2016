###################################################################################
#
#  -file = "DEoptim_rahm_model.R"   Code written March 2014
#  - Author: Kelsey Ruckert (klr324@psu.edu)
#
#  -This function sources the Rahmstorf(2007) model to be sourced in the DEoptim
#       R function. Deoptim is used in the bootstrap and MCMC codes as described
#       in Ruckert et al. (accepted) to find good initial values for each of the
#       parameters.
#
# THIS CODE IS PROVIDED AS-IS WITH NO WARRANTY (NEITHER EXPLICIT
# NOT IMPLICIT).  I SHARE THIS CODE IN HOPES THAT IT IS USEFUL,
# BUT I AM NOT LIABLE FOR THE BEHAVIOR OF THIS CODE IN YOUR OWN
# APPLICATION.  YOU ARE FREE TO SHARE THIS CODE SO LONG AS THE
# AUTHOR(S) AND VERSION HISTORY REMAIN INTACT.
#
# To use this function, simply source this file:
#   source("DEoptim_rahm_model.R")
#
# INPUTS:
#   vector of parameter values:
#       p[1]: alpha; sensitivity of sea-level to temperature changes
#       p[2]: T_0; temperature when sea-level anomaly is zero
#       p[3]: H_0; initial sea-level anomaly
#
###################################################################################

model = function(p){# p represents the parameters in a vector
    
    # Estimate the rate, dH, of sea-level change each year, Equation (S17)
    dH = p[1]*(hist.temp - p[2])
    #p[1] = sensitivity of sea-level to temperature changes
    #p[2] = temperature when sea-level is zero

    # Set up empty vector for sea level anomalies.
    H_1 = rep(NA, to)
    H_1[1] = p[3]     # sea-level in 1880
    
    # Run a forward euler to estimate sea-level over time
    for (i in from:to){
        H_1[i] = H_1[i-1] + dH[i-1]*timestep
    }
    
    # Return sea-level anomalies
    return(H_1)
}

#################################### END ##########################################
