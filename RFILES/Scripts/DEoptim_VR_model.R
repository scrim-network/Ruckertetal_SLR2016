###################################################################################
#
#  -file = "DEoptim_VR_model.R"   Code written Sept. 16 2016
#  - Author: Kelsey Ruckert (klr324@psu.edu)
#
#  -This function sources the Vermeer and Rahmstorf(2009) model to be sourced in the DEoptim
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
#   source("DEoptim_VR_model.R")
#
###################################################################################

VR_model = function(p){ # p represents the parameters in a vector
  #p[1] = sensitivity of sea-level to temperature changes
  #p[2] = equilibrium Sea-level anomaly [m] when the temperature anomaly = 0
  #p[3] = initial value of sea-level in 1880 [m]
  #p[4] = rapid response of sea-level term
  
  dS     = rep(NA, to) 
  S_1    = rep(NA, to)
  dT     = rep(NA, to)
  
  # Estimate the rate of temperature change with smoothing
  temp.dif = diff(hist.temp)
  
  for(i in 2:(length(hist.temp)-1)){
    dT[i] = 0.5*(temp.dif[i-1] + temp.dif[i])
  }
  dT[1]  = temp.dif[1] - (0.5*(temp.dif[2] - temp.dif[1]))
  dT[length(hist.temp)] = temp.dif[length(hist.temp)-1] + (0.5*(temp.dif[length(hist.temp)-1] - temp.dif[length(hist.temp)-2]))
  
  # Estimate the rate of sea-level change each year
  dS = p[1]*(hist.temp - p[2]) + p[4]*dT
  
  # sea-level values
  S_1[1] = p[3] # sea-level in 1880
  
  for(i in from:to){
    # Run a forward euler to estimate sea-level over time
    S_1[i] = S_1[i-1] + dS[i-1]*timestep
  }
  return(S_1)
}

#################################### END ##########################################