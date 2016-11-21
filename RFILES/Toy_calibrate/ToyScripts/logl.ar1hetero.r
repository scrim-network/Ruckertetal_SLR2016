# ### create matrix for AR1 covariance matrix
# times <- 1:5
# rho <- 0.5
# sigma <- 2
# sigma^2/(1-rho^2)
# ###############
# H <- abs(outer(times, times, "-"))
# V <- sigma^2/(1-rho^2) * rho^H
# V

logl.ar1 <-
  function(r,sigma1,rho1,eps1) # default obs error is 0. sigma1 is standard error. 
  {
    n = length(r)
    sigma.proc = sigma1/sqrt(1-rho1^2) # stationary process variance sigma.proc^2
    if(all(eps1==0)){
      logl = dnorm(r[1],sd=sigma.proc,log=TRUE)
      if(n>1) {
        w = r[2:n] - rho1*r[1:(n-1)] # whitened residuals
        # logl = logl + sum(dnorm(w,sd=sigma1,log=TRUE))
        # This is what we had before to make the computation faster. 
        # This approximation should not change the result, but it is worth trying 
        logl = logl + sum(dnorm(w,sd=sqrt(sigma1^2+eps1[-1]^2),log=TRUE))
      }
    }else{
      H <- abs(outer(1:n, 1:n, "-"))
      v = sigma.proc^2*rho1^H
      v = v+diag(eps1^2)
      # Need R package "mvtnorm"       
      logl = dmvnorm(r,sigma=v,log=TRUE)
    }
    return(logl)
  }
