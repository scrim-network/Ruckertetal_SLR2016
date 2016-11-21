###################################################################################
#
#  -file = "sealevel_grinsted_model.R"   Code written Sept. 16 2016
#  - Author: Kelsey Ruckert (klr324@psu.edu)
#
#  -This function sources the Grinsted et al. (2010) model to be sourced into the uncertainty methods as described
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
#   source("sealevel_grinsted_model.R")
#
###################################################################################

grinsted_sealevel = function(parameters, Temp){ #inputs are parameters and temperature
  model.p = length(parameters)    # number of parameters in the model
  a = parameters[1]               # sensitivity of sea-level to temperature changes [m/C]
  b = parameters[2]               # equilibrium Sea-level anomaly [m] when the temperature anomaly = 0
  initialvalue = parameters[3]    # initial value of sea-level in 1880 [m]
  Tau = parameters[4]             # response time (1/time-scale of exponential decay (e-folding time) [years^-1]m)
  
  # Equilibrium sea level
  Seq = a*Temp + b
  
  sea_rate = rep(NA,to) 
  values = rep(NA,to)
  values[1] = initialvalue
  sea_rate[1] = (Seq[1] - values[1])*Tau
  
  for(i in from:to){
    # Estimate the rate of sea-level change each year
    sea_rate[i] = (Seq[i-1] - values[i-1])*Tau
    
    # Run a forward euler to estimate sea-level over time
    values[i] = values[i-1] + sea_rate[i-1]*timestep
  }
  #return sea-level, sea-level rates, and number of parameters
  return(list(sle = values, slrate = sea_rate, model.p = model.p))
}

################################### END ##########################################