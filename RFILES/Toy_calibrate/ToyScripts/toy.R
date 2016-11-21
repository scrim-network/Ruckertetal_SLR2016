##################################################################################
#
#  -file = "toy.R"   Code written February 2015
#  - Author: Kelsey Ruckert (klr324@psu.edu)
#
#  -This function sources the toy model used to test where the surprise index fails
#       as described in Ruckert et al. (2016). For further description and references,
#       please read the paper and appendix.
#
# THIS CODE IS PROVIDED AS-IS WITH NO WARRANTY (NEITHER EXPLICIT
# NOT IMPLICIT).  I SHARE THIS CODE IN HOPES THAT IT IS USEFUL, 
# BUT I AM NOT LIABLE FOR THE BEHAVIOR OF THIS CODE IN YOUR OWN
# APPLICATION.  YOU ARE FREE TO SHARE THIS CODE SO LONG AS THE
# AUTHOR(S) AND VERSION HISTORY REMAIN INTACT.
#
# To use this function, simply source this file:
#   source("toy.R")
#
###################################################################################

model<-function(parm,x){ # Inouts are parameters and length of data
    model.p=length(parm) # number of parameters in the model
    
  a=parm[1]
  b=parm[2]
  
  y.mod <- a+b*x #This linear equation represent the model
  #Return model observations and number of parameters
  return(list(mod.obs = y.mod, model.p = model.p))
}
#################################### END #############################################
