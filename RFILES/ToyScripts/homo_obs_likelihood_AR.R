#################################################################################
#
#  -file = "homo_obs_likelihood_AR.R"   Code written July 2014
#  - Author: Yawen Guan (yig5031@psu.edu)
#  - Edits by: Kelsey Ruckert (klr324@psu.edu) February 2015
#
#  -This function computes the log likelihood for a zero-mean AR1 process from
#       observations. The observations include i.i.d. and AR(1) assuming
#       homoskedastic or heteroskedastic. The tests using this program are
#       described in the appendix of Ruckert et al. (2016). For further
#       description and references, please read the paper
#       and the appendix.
#
#   -NOTE: Before sourcing this code make sure the code is set to run for either
#       Ar(1) observations or iid observations. Descriptions of how to use this
#       for other observations and models can be found in the R package in review
#       "VAR1"
#
# THIS CODE IS PROVIDED AS-IS WITH NO WARRANTY (NEITHER EXPLICIT
# NOT IMPLICIT).  I SHARE THIS CODE IN HOPES THAT IT IS USEFUL, 
# BUT I AM NOT LIABLE FOR THE BEHAVIOR OF THIS CODE IN YOUR OWN
# APPLICATION.  YOU ARE FREE TO SHARE THIS CODE SO LONG AS THE
# AUTHOR(S) AND VERSION HISTORY REMAIN INTACT.
#
# To use this function, simply source this file:
#   source("toy_obs_likelihood_AR.R")
#
################################################################################

# source in the function that finds the log likelihood and the lag-1
#   autocorrelation coefficient

# Use either the approximation method. In the case of small errors the approximation method produces the same results.
source("ToyScripts/ar.R")

#- OR the less stable large matrix method. NOTE: This method takes significantly longer.
#source("ToyScripts/logl.ar1hetero.r")

log.lik = function(p) # model.p is the dimension of model parameters 
{ 
  par=p[1:model.p]
  
  sigma.y = p[model.p+1] # sigma is the innovation variance
  rho.y = p[model.p+2]   # rho1 is the lag-1 autocorrelation coefficient
  
  model.out = model(par, x) # run the model and name to output y.mod
  y.mod = model.out$mod.obs
  
  llik.y  = 0
  resid.y = obs-y.mod # estimate the residuals from the model simulation

  ##---------------------------- STOP ---------------------------------------##
  #Set the likelihood for AR(1) observations
  llik.y  = logl.ar1(resid.y, sigma.y, rho.y,y.meas.err) # AR(1) observations
  ##-------------------------------------------------------------------------##
  
  llik = llik.y # assume residuals are independent
  llik  
}
# (log) prior distribution:
log.pri = function(p)
{
  par=p[1:model.p]
  
  sigma.y = p[model.p+1]
  rho.y = p[model.p+2]
  
  in.range = all(p > bound.lower) & all(p < bound.upper)
  
  if(in.range) {
    lpri=0
  } else {
    lpri = -Inf
  }
  
  lpri
}

# (log) posterior distribution:  posterior ~ likelihood * prior
log.post = function(p)
{  
  lpri = log.pri(p)
  if(is.finite(lpri)) { # evaluate likelihood if nonzero prior probability
    lpost = log.lik(p) + lpri
  } else {
    lpost = -Inf  
  }
  lpost
}


# calculate a scale matrix to transform a iid normal distribution,
# which defines a multivariate normal proposal distribution for MCMC
# (the proposal distribution is proportional to the covariance of
# a preliminary Markov chain of the posterior distribution to sample,
# tuned to be optimally scaled if the posterior is multivariate normal,
# plus a small nugget term to ensure ergodicity)
#
# from Gareth O. Roberts and Jeffrey S. Rosenthal,
# "Examples of adaptive MCMC", unpublished
proposal.matrix = function(prechain, mult=1, beta=0.05)
{
  # mult = overall scale factor to adjust all step sizes
  # beta = relative influence of nugget term
  
  p = ncol(prechain)
  precov = cov(prechain)
  
  propcov = (1-beta)*2.38^2*precov/p + beta*0.1^2*diag(p)/p
  propcov = mult*propcov
  
  mat = t(chol(propcov))
}

################################## END ###########################################
