#################################################################################
#
#  -file = "ar.R"   Code written July 2014
#  - Author: Yawen Guan (yig5031@psu.edu)
#
#  -This function finds the log likelihood for a zero-mean AR1 process and the
#       simulated lag-1 autocorrelation coefficient for the model. Both outputs are
#       used in the MCMC likelihood programs as described in Ruckert et al. (GRL 2016). For
#       further description and references, please read the paper and the appendix.
#
#   -NOTE: The innovation variance = sigma^2 and the lag-1 autocorrelation
#       coefficient = rho1. In this function the initial values are ignored to
#       reduce bias and to make consistent with a VAR code. Description of a VAR
#       code can be found in the R package in review "VAR1"
#
# THIS CODE IS PROVIDED AS-IS WITH NO WARRANTY (NEITHER EXPLICIT
# NOT IMPLICIT).  I SHARE THIS CODE IN HOPES THAT IT IS USEFUL,
# BUT I AM NOT LIABLE FOR THE BEHAVIOR OF THIS CODE IN YOUR OWN
# APPLICATION.  YOU ARE FREE TO SHARE THIS CODE SO LONG AS THE
# AUTHOR(S) AND VERSION HISTORY REMAIN INTACT.
#
# To use this function, simply source this file:
#   source("ar.R")
#
###################################################################################

# Estimate the log likelihood of the AR1 process
logl.ar1 = function(r,sigma1,rho1,eps1=y.meas.err) # default obs error is 0
{
    n = length(r) # r is the residuals
    
  #sigma.proc = sigma1/sqrt(1-rho1^2) # stationary process variance sigma.proc^2
  #logl = dnorm(r[1],sd=sigma.proc+eps1[1],log=TRUE)
  
  logl=0
  if(n>1) {
    w = r[2:n] - rho1*r[1:(n-1)] # this process whitens the residuals
    logl = logl + sum(dnorm(w,sd=sqrt((sigma1)^2+(eps1[c(-1)])^2),log=TRUE)) # add in the sum of
           # density of the whitened residuals with a standard deviation of the
           # variance and the obs. errors
  }
  return(logl)
}

# Simulate the lag-1 autocorrelation coefficient
ar1.sim = function(N,rho1,sigma,eps1=y.meas.err) {
	x = rep(NA,N)
  
	x[1] = sigma/sqrt(1-rho1^2)
	if(N>length(eps1)) eps1=c(eps1,rep(eps1[n],N-n))
    for(i in 2:N)
		x[i] = rho1*x[i-1] + rnorm(1,sd=sqrt(sigma^2+eps1[i]^2))
        return(x) # x returns the lag-1 autocorrelation coefficients
}

################################## END ############################################