#################################################################################
#
#  -file = "Toy_MCMC_IID.R"   Code written February 2015
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
#       -Workspace file: "test_IID_data.RData"
#
#  -The spreadsheet: "testhypoth_surprise.csv" is a previously calculated surprise index for
#       all the following tests.
#
###################################################################################
rm(list =ls()) #Clear global environment
library(compiler)
enableJIT(3)
# set.seed(1234)
# set.seed(111)

#source("Scripts/put_fig_letter.R")
#----------------------- Test #1: Length of Observations ----------------------------#
################# Step 1: Choose the number of observations
#datalength = 1
#datalength = 5
#datalength = 20
#datalength = 50
#datalength = 100
datalength = 200
x = 1:datalength

################# Step 2: Set up the true model & Generate Random iid numbers
true_par = c(1,0.2) #a=1 b=0.2
source("ToyScripts/toy.R") # this is the model: y = a+b*x
true.obs= model(true_par,x)

residuals = rnorm(datalength,mean=0,sd=0.5) #Calculate the residuals
# observations are: obs =  a + b*x + residuals
obs = true.obs$mod.obs + residuals

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

res = obs-y.obs$mod.obs #Estimate the residuals

# apply auto-correlation to determine iid correlation coefficients
rho=rep(NA,3)
ac=acf(res, lag.max=5, plot=TRUE, main="")
rho[1]=ac$acf[1]
rho[2]=ac$acf[2]
rho[3]=ac$acf[3]
rho[4]=ac$acf[4]
rho[5]=ac$acf[5]

################# Step 5:  Set up and Run MCMC for iid observations
bound.lower = c(-10,-10,0)
bound.upper = c(10,10,1)
y.meas.err=rep(0,length(y.obs)) # assuming homoskedastic error

model.p=2
parnames=c("a","b","sigma.y")
source("ToyScripts/toy.R")
source("ToyScripts/iid_obs_likelihood.R")

p = c(result$par[1], result$par[2], sd(res))
p0 = c(1,0.3,0.6) # Initial parameter guesses
p0 = optim(p0, function(p) -log.post(p))$par
print(round(p0,4))
library(mcmc)

step = c(0.1,0.001,0.1) # for datasets larger than 20 observations
#step = c(0.1,0.1,0.1) # for datasets equal to or less than 20 observations
NI = 1E6 #number of iterations
burnin = seq(1,0.01*NI,1) # 1% burnin

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

## Check trace plots for convergence ##
par(mfrow=c(2,2))
plot(prechain1[-burnin,1], type="l", ylab="a", xlab="Number of Runs", main="")
plot(prechain1[-burnin,2], type="l", ylab="b", xlab="Number of Runs", main="")
plot(prechain1[-burnin,3], type="l", ylab="Sigma", xlab="Number of Runs", main="")

################# Step 6:  Estimate the Hindcast
n=length(prechain1[,1])
ssprechain1 = prechain1[seq(length(burnin),n,100),] #Only a subsequent number is neccessary (~10,000)
h = length(ssprechain1[,1])

#If the dataset length is one observation uncomment the following line
# new.mcmc_hind=prechain1[,1]+prechain1[,2]*x ## For one point

new.mcmc_hind=mat.or.vec(h,datalength) #More than 1 observation
for(i in 1:datalength){
  new.mcmc_hind[,i]=ssprechain1[,1]+ssprechain1[,2]*x[i]
}

#Superimpose the residuals on the hindcast fits #
hindcast=new.mcmc_hind+rnorm(datalength,mean=0,sd=ssprechain1[,3]) ## For IID

################# Step 7:  Analyze the Results
# Save the probability density functions of the parameters
pdfa = density(prechain1[-burnin,1])
pdfb = density(prechain1[-burnin,2])
twoh_a = density(prechain1[-burnin,1]) ; twoh_b = density(prechain1[-burnin,2])

#Calculate the 90% Confidence Interval
# five = quantile(hindcast,0.05) #for one data point
# ninef = quantile(hindcast,0.95)
# range_x=c(five, rev(ninef)); range_y=c(x, rev(x))

five <-
  ninetyfive <-rep(NA,datalength) # For multiple observations
for(i in 1:datalength){
  five[i] = quantile(hindcast[,i],0.05)
  ninetyfive[i] = quantile(hindcast[,i],0.95)
}
range_x=c(five, rev(ninetyfive)); range_y=c(x, rev(x))

#Calculate the Surprise index by observation the ratio of outliers to the expected # of outliers
# The confidence intervals in the paper are:
# Percentages: 10%, 20%, 30%, 40%, 50%, 60%, 70%, 80%, 90%, 92%, 95%,  96%, 97%,  98%, 99%,  100%
# low numbers: 0.45,0.40,0.35,0.30,0.25,0.20,0.15,0.10,0.05,0.04,0.025,0.02,0.015,0.01,0.005,0
# high numbers:0.55,0.60,0.65,0.70,0.75,0.80,0.85,0.90,0.95,0.96,0.975,0.98,0.985,0.99,0.995,1

par(mfrow=c(1,1), mgp=c(1.5,.5,0), mar=c(4, 4, 3, 1))
plot(x, obs, type="l",xlab="x",ylab="Observations") #, ylim=c(-2,6))
polygon(range_y, range_x, col="red")
points(x, obs, pch=20, cex=0.5)
legend("topleft", c("90% CI","Observations"), pch=c(15,20),bty="n", col=c("red","black"))

#IID
IID.005 <-
  IID.995 <- rep(NA,datalength) # 99%
IID.45 <- IID.995; IID.55 <- IID.995 # 10%
IID.40 <- IID.995; IID.60 <- IID.995 # 20%
IID.35 <- IID.995; IID.65 <- IID.995 # 30%
IID.30 <- IID.995; IID.70 <- IID.995 # 40%
IID.25 <- IID.995; IID.75 <- IID.995 # 50%
IID.20 <- IID.995; IID.80 <- IID.995 # 60%
IID.15 <- IID.995; IID.85 <- IID.995 # 70%
IID.10 <- IID.995; IID.90 <- IID.995 # 80%
IID.5 <- IID.995; IID.95 <- IID.995  # 90%
IID.4 <- IID.995; IID.96 <- IID.995  # 92%
IID.025 <- IID.995; IID.975 <- IID.995 # 95%
IID.02 <- IID.995; IID.98 <- IID.995   # 96%
IID.015 <- IID.995; IID.985 <- IID.995 # 97%
IID.01 <- IID.995; IID.99 <- IID.995   # 98%
IID.0 <- IID.995; IID.100 <- IID.995   # 100%

for(i in 1:datalength){
  IID.45[i] = quantile(hindcast[,i],0.45); IID.55[i] = quantile(hindcast[,i],0.55)
  IID.40[i] = quantile(hindcast[,i],0.40); IID.60[i] = quantile(hindcast[,i],0.60)
  IID.35[i] = quantile(hindcast[,i],0.35); IID.65[i] = quantile(hindcast[,i],0.65)
  IID.30[i] = quantile(hindcast[,i],0.30); IID.70[i] = quantile(hindcast[,i],0.70)
  IID.25[i] = quantile(hindcast[,i],0.25); IID.75[i] = quantile(hindcast[,i],0.75)
  IID.20[i] = quantile(hindcast[,i],0.20); IID.80[i] = quantile(hindcast[,i],0.80)
  IID.15[i] = quantile(hindcast[,i],0.15); IID.85[i] = quantile(hindcast[,i],0.85)
  IID.10[i] = quantile(hindcast[,i],0.10); IID.90[i] = quantile(hindcast[,i],0.90)
  IID.5[i] = quantile(hindcast[,i],0.05); IID.95[i] = quantile(hindcast[,i],0.95)
  IID.4[i] = quantile(hindcast[,i],0.04); IID.96[i] = quantile(hindcast[,i],0.96)
  IID.025[i] = quantile(hindcast[,i],0.025); IID.975[i] = quantile(hindcast[,i],0.975)
  IID.02[i] = quantile(hindcast[,i],0.02); IID.98[i] = quantile(hindcast[,i],0.98)
  IID.015[i] = quantile(hindcast[,i],0.015); IID.985[i] = quantile(hindcast[,i],0.985)
  IID.01[i] = quantile(hindcast[,i],0.01); IID.99[i] = quantile(hindcast[,i],0.99)
  IID.005[i] = quantile(hindcast[,i],0.005); IID.995[i] = quantile(hindcast[,i],0.995)
  IID.0[i] = quantile(hindcast[,i],0); IID.100[i] = quantile(hindcast[,i],1)
}

range_IID_10=c(IID.45, rev(IID.55)); range_IID_20=c(IID.40, rev(IID.60))
range_IID_30=c(IID.35, rev(IID.65)); range_IID_40=c(IID.30, rev(IID.70))
range_IID_50=c(IID.25, rev(IID.75)); range_IID_60=c(IID.20, rev(IID.80))
range_IID_70=c(IID.15, rev(IID.85)); range_IID_80=c(IID.10, rev(IID.90))
range_IID_90=c(IID.5, rev(IID.95)); range_IID_92=c(IID.4, rev(IID.96))
range_IID_95=c(IID.025, rev(IID.975)); range_IID_96=c(IID.02, rev(IID.98))
range_IID_97=c(IID.015, rev(IID.985)); range_IID_98=c(IID.01, rev(IID.99))
range_IID_99=c(IID.005, rev(IID.995)); range_IID_100=c(IID.0, rev(IID.100))

plot(x, obs, type="l",xlab="x",ylab="Observations") #, ylim=c(-2,6))
points(x, obs, pch=20, cex=0.85)
polygon(range_y, range_IID_10, col="red")
polygon(range_y, range_IID_20, col="red")
polygon(range_y, range_IID_30, col="red")
polygon(range_y, range_IID_40, col="red")
polygon(range_y, range_IID_50, col="red")
polygon(range_y, range_IID_60, col="red")
polygon(range_y, range_IID_70, col="red")
polygon(range_y, range_IID_80, col="red")
polygon(range_y, range_IID_90, col="red")
polygon(range_y, range_IID_92, col="red")
polygon(range_y, range_IID_95, col="red")
polygon(range_y, range_IID_96, col="red")
polygon(range_y, range_IID_97, col="red")
polygon(range_y, range_IID_98, col="red")
polygon(range_y, range_IID_99, col="red")
polygon(range_y, range_IID_100, col="red")

################# Step 8:  Save the results from using mulitple datalengths:
#After running each analysis varying the observations save the pdf
one_a = density(prechain1[-burnin,1]) ; one_b = density(prechain1[-burnin,2])
five_a = density(prechain1[-burnin,1]) ; five_b = density(prechain1[-burnin,2])
twent_a = density(prechain1[-burnin,1]) ; twent_b = density(prechain1[-burnin,2])
fity_a = density(prechain1[-burnin,1]) ; fity_b = density(prechain1[-burnin,2])
hund_a = density(prechain1[-burnin,1]) ; hund_b = density(prechain1[-burnin,2])
twoh_a = density(prechain1[-burnin,1]) ; twoh_b = density(prechain1[-burnin,2])

#After running each analysis varying the observations save the 90% Confidence Interval
one_x = range_x ; one_y = range_y
five_x = range_x ; five_y = range_y
twent_x = range_x ; twent_y = range_y
fity_x = range_x ; fity_y = range_y
hund_x = range_x ; hund_y = range_y
twoh_x = range_x ; twoh_y = range_y

#Save the workspace:
# save.image(file = "Workspace/testiidall.RData")
save.image(file = "Workspace/test_IID_data.RData")
