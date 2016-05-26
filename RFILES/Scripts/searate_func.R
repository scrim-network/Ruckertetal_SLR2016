###################################################################################
#
#  -file = "searate_func.R"   Code written March 2014
#  - Author: Kelsey Ruckert (klr324@psu.edu)
#
#  -This function sources the Rahmstorf(2007) model to find the rate of change in
#       sea-level per year as described in Ruckert et al. (GRL 2016). For further
#       description and references, please read the paper.
#
# THIS CODE IS PROVIDED AS-IS WITH NO WARRANTY (NEITHER EXPLICIT
# NOT IMPLICIT).  I SHARE THIS CODE IN HOPES THAT IT IS USEFUL,
# BUT I AM NOT LIABLE FOR THE BEHAVIOR OF THIS CODE IN YOUR OWN
# APPLICATION.  YOU ARE FREE TO SHARE THIS CODE SO LONG AS THE
# AUTHOR(S) AND VERSION HISTORY REMAIN INTACT.
#
# To use this function, simply source this file:
#   source("searate_func.R")
#
###################################################################################

thermal = function(parameters, Temp){ #inputs are parameters and temperature
    a = parameters[1] # sensitivity of sea-level to temperature changes
    Ti = parameters[2] # temperature when sea-level is zero
    therm_rate=a*(Temp-Ti)
    list(therm_rate) #therm_rate returns the rates of sea-level change each year
}

#################################### END #########################################