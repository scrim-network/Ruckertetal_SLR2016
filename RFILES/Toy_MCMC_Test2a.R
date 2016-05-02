#################################################################################
#
#  -file = "Toy_MCMC_Test2a.R"   Code written February 2015
#  - Author: Kelsey Ruckert (klr324@psu.edu)
#
#  -The toy scripts analyze three hypotheses for failure in the surprise index
#       as described in Ruckert et al. (2016). The three hypotheses include number
#       of observations, AR(1) observations, and observations assuming heteroskedastic
#       errors. It also produces the last three supplementary figures shown in that
#       paper. The graphs will be saved as pdf files in the current working directory.
#       For further description and references, please read the paper and appendix.
#
# THIS CODE IS PROVIDED AS-IS WITH NO WARRANTY (NEITHER EXPLICIT
# NOT IMPLICIT).  I SHARE THIS CODE IN HOPES THAT IT IS USEFUL,
# BUT I AM NOT LIABLE FOR THE BEHAVIOR OF THIS CODE IN YOUR OWN
# APPLICATION.  YOU ARE FREE TO SHARE THIS CODE SO LONG AS THE
# AUTHOR(S) AND VERSION HISTORY REMAIN INTACT.
#
#   -NOTE: This program applies a Markov Chain Monte Carlo method with varying
#       assumptions about the observations including:
#
#       -varies number of observations
#       -preserves the iid structure of the observations
#       -creates random realizations of "natural variability"
#       -preserves the AR(1) structure of the observations
#       -Markov Chain Monte Carlo AR(1) Heteroskedastic errors
#       -Markov Chain Monte Carlo AR(1) Homoskedastic errors
#
#  -The output from this file once ran will be saved as:
#       -Workspace file: "testiidall.RData"
#
#  -The spreadsheet: "testhypoth_surprise.csv" is a previously calculated surprise index for
#       all the following tests.
#
###################################################################################
rm(list =ls()) #Clear global environment
library(compiler)
enableJIT(3)

#source("Scripts/put_fig_letter.R")
#----------------------- Test #1: Length of Observations ----------------------------#
load("Workspace/testiidall.RData") #Load results from test #1

#------------ Test #2: Adding in Autocorrelation & Test 3: Heteroskedastic Errors -------------#
################# Step 1: Set up the true model & Generate generate AR 1 errors
#datalength = 1
#datalength = 5
#datalength = 20
#datalength = 50
datalength = 200
x = 1:datalength

true_par = c(1,0.2) #a=1 b=0.2
source("ToyScripts/toy.R") # this is the model: y = a+b*x
true.obs= model(true_par,x)

rho=0.3  # define your correlation
sigma.ar1=0.3 # define your ar1 white noise inovation variance
sigma.inov = sigma.ar1/(1-0.3^2) # this is your stationary variance for the AR1 process . i.e.the variance for time 0

t.s0 = rnorm(1,0,sd = sqrt(sigma.inov))
t.s<-c()
t.s[1]<-t.s0
for ( i in (2:datalength)){ # write a loop to simulate AR1 residuals
  t.s[i]<-rho*t.s[i-1]+rnorm(1,mean=0,sd=sqrt(sigma.ar1))
}

obs = true.obs$mod.obs + t.s #Add the observations and the residuals

################# Step 3: Find initial parameter guesses
fn = function(parameters, x, obs){
  #set parameter (a,b)
  a=parameters[1]
  b=parameters[2]
  #set up equation: a+bx= observations  
  data = a+b*x
  #set up equation to find residuals
  resid=obs-data
  #Find root mean square error
  rmse = sqrt(mean(resid^2))
  # return the root mean square error
  return(rmse)
}

result = optim(true_par, fn, gr=NULL, x, obs)

################# Step 4:  Simulate the observations
source("ToyScripts/toy.R")
parameter = c(result$par[1],result$par[2]) #Use the initial parameter guesses
y.obs= model(parameter,x)

#Plot the observations
#plot(x,obs, pch=20, xlab="x", ylab="Data")
#lines(x, y.obs$mod.obs, col="red", lwd=3)
#lines(x,err_pos, col="slateblue", lwd=3)
#lines(x,err_neg, col="slateblue", lwd=3)

res = obs-y.obs$mod.obs #Estimate the residuals

# apply auto-correlation to determine iid correlation coefficients
rho=rep(NA,3)
ac=acf(res, lag.max=5, plot=TRUE, main="")
rho[1]=ac$acf[1]
rho[2]=ac$acf[2]
rho[3]=ac$acf[3]
rho[4]=ac$acf[4]
rho[5]=ac$acf[5]

################# Step 5:  Set up and Run MCMC
##---------------------- STOP ---------------------------------------##
# Make sure the likelihood is set for AR(1) observations! in "toy_obs_likelihood_AR.R"
##-------------------------------------------------------------------##
######################### Set up and Run MCMC #######################
bound.lower = c(-10,-10,0,-0.99)
bound.upper = c(10,10,1,0.99)
y.meas.err=rep(0,length(y.obs)) #Error for homoskedastic
model.p=2
parnames=c("a","b","sigma.y","phi11")
source("ToyScripts/toy.R")
source("ToyScripts/toy_obs_likelihood_AR_1.R")

p = c(result$par[1], result$par[2], sd(res), rho[2])
p0 = c(1,0.3,0.6, 0.5)
p0 = optim(p0, function(p) -log.post(p))$par
print(round(p0,4))
library(mcmc)

step = c(0.1,0.001,0.01,0.01)
NI = 5E6 #number of iterations
burnin = seq(1,0.01*NI,1)

mcmc.out1 = metrop(log.post, p0, nbatch=NI, scale=step)
prechain1 = mcmc.out1$batch
mcmc.out1$accept
# Calculate the parameter acceptance rate
acceptrate = mcmc.out1$accept * 100
#Print the acceptance rate as a percent. Should be higher than 10% and no greater than 60%
cat("Accept rate =", acceptrate, "%\n")

new = c(median(prechain1[-burnin,1]), median(prechain1[-burnin,2]))
print(new) # Print MCMC best estimated parameters
new.est = model(new, x)

### Check trace plots for convergence
par(mfrow=c(2,2))
plot(prechain1[-burnin,1], type="l", ylab="a", xlab="Number of Runs", main="")
plot(prechain1[-burnin,2], type="l", ylab="b", xlab="Number of Runs", main="")
plot(prechain1[-burnin,3], type="l", ylab="Sigma", xlab="Number of Runs", main="")
plot(prechain1[-burnin,4], type="l", ylab="Rho", xlab="Number of Runs", main="")

################# Step 6:  Estimate the Hindcast
n=length(prechain1[,1])
ssprechain1 = prechain1[seq(length(burnin),n,500),] #Only a subsequent number is neccessary (~10,000)
h = length(ssprechain1[,1])

#Calculate the simulated model hindcasts
new.mcmc_hind=mat.or.vec(h,datalength)
for(i in 1:datalength){
  new.mcmc_hind[,i]=ssprechain1[,1]+ssprechain1[,2]*x[i]
}
# Calculate the residuals with the MCMC estimates AR(1) coefficient and the sigma
res.mcmc = mat.or.vec(h, datalength) #(nr,nc)
hindcast = res.mcmc
for(k in 1:h) {
  for(i in 2:datalength) {
    res.mcmc[k,i] = ssprechain1[k,4]*res.mcmc[k,i-1] + rnorm(1,mean=0,sd=ssprechain1[k,3]) ## For AR1
  }
}
for(i in 1:h) {
  hindcast[i,]=new.mcmc_hind[i,]+res.mcmc[i,] #Add the model simulations and the residuals
}

################# Step 7:  Analyze the Results
# Save the probability density functions of the parameters
pdfa = density(prechain1[,1])
pdfb = density(prechain1[,2])

#Calculate the 90% Confidence Interval
#five = quantile(hindcast,0.05) #for one data point
#ninef = quantile(hindcast,0.95)
five <-
  ninetyfive <-rep(NA,datalength) # For multiple observations
for(i in 1:datalength){
  five[i] = quantile(hindcast[,i],0.05)
  ninetyfive[i] = quantile(hindcast[,i],0.95)
}
range_x=c(five, rev(ninetyfive)); range_y=c(x, rev(x))

#Calculate the Surprise index by observation the ratio of outliers to the expected # of outliers
par(mfrow=c(1,1), mgp=c(1.5,.5,0), mar=c(4, 4, 3, 1))
plot(x, obs, type="l",xlab="x",ylab="Observations") #, ylim=c(-2,6))
polygon(range_y, range_x, col="red")
points(x, obs, pch=20)
legend("topleft", c("90% CI","Observations"), pch=c(15,20),bty="n", col=c("red","black"))

################# Step 8:  Save the Homoskedastic Results
# Save the pdf of the parameters
homtwoh_a = density(prechain1[,1]) ; homtwoh_b = density(prechain1[,2])
#Save the 90% confidence interval
homtwoh_x = range_x ; homtwoh_y = range_y

#Again save the workspace:
save.image(file = "Workspace/testiidall.RData")
