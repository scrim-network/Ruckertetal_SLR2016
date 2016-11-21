###################################################################################
#
#  -file = "sealevel_rahm_model.R"   Code written March 2014
#  - Author: Kelsey Ruckert (klr324@psu.edu)
#
#  -This function sources the Rahmstorf(2007) model to be sourced into the uncertainty methods as described
#       in Ruckert et al. (2016). For further
#       description and references, please read the paper.
#
# THIS CODE IS PROVIDED AS-IS WITH NO WARRANTY (NEITHER EXPLICIT
# NOT IMPLICIT).  I SHARE THIS CODE IN HOPES THAT IT IS USEFUL,
# BUT I AM NOT LIABLE FOR THE BEHAVIOR OF THIS CODE IN YOUR OWN
# APPLICATION.  YOU ARE FREE TO SHARE THIS CODE SO LONG AS THE
# AUTHOR(S) AND VERSION HISTORY REMAIN INTACT.
#
# To use this function, simply source this file:
#   source("sealevel_rahm_model.R")
#
# INPUTS:
#   vector of parameter values:
#       parameters[1]: alpha; sensitivity of sea-level to temperature changes
#       parameters[2]: T_0; temperature when sea-level anomaly is zero
#       parameters[3]: H_0; initial sea-level anomaly
#
#   vector of temperatures:
#
###################################################################################

rahmfunction = function(parameters, Temp){ #inputs are parameters and temperature
    
    # Determine number of model parameters
    model.p = length(parameters)
    
    # Extract parameter values
    a = parameters[1]             # sensitivity of sea-level to temperature changes
    Ti = parameters[2]            # temperature when sea-level is zero
    initialvalue = parameters[3]  # initial value of sea-level in 1880
    
    # Estimate the rate of sea-level change each year, Equation (S17)
    rate = a*(Temp - Ti)
    
    # Set up empty vector for sea level anomalies.
    values = rep(NA, to)
    values[1] = initialvalue
    
    # Run a forward euler to estimate sea-level over time
    for(i in from:to){
        values[i] = values[i-1]+rate[i-1]*timestep
    }
    
    # Return sea-level anomalies, sea-level rates, and number of parameters.
    return(list(sle = values, slrate = rate, model.p = model.p))
}

################################### END ##########################################
