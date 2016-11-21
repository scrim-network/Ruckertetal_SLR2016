#################################################################################
#
#  -file = "testAR1vsIID.R"   Code written February 2015
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
#       -Workspace file: "test_AR1vsIID.RData"
#
#  -The spreadsheet: "testhypoth_surprise.csv" is a previously calculated surprise index for
#       all the following tests.
#
###################################################################################
rm(list =ls()) #Clear global environment
library(compiler)
enableJIT(3)
set.seed(1234)
# set.seed(111)
source("Scripts/put_fig_letter.R")

# 1: Set up the true model & Generate generate AR 1 errors
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

# 2: Find initial parameter guesses
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

# 3:  Simulate the observations
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

####### Set up and Run MCMC for AR(1) homoskedastic observations assuming UNKNOWN RHO ---------------------------------
bound.lower = c(-10,-10,0,-0.99)
bound.upper = c(10,10,1,0.99)
y.meas.err=rep(0,length(y.obs)) #Error for homoskedastic
model.p=2
parnames=c("a","b","sigma.y","phi11")
source("ToyScripts/toy.R")
source("ToyScripts/homo_obs_likelihood_AR.R")

p = c(result$par[1], result$par[2], sd(res), rho[2])
p0 = c(1,0.3,0.6, 0.5)
p0 = optim(p0, function(p) -log.post(p))$par
print(round(p0,4))
library(mcmc)

step = c(0.1,0.001,0.01,0.01)
NI = 5E6 #number of iterations
burnin = seq(1,0.01*NI,1)

AR1.mcmc = metrop(log.post, p0, nbatch=NI, scale=step)
AR1.mcmc.chain = AR1.mcmc$batch
AR1.mcmc$accept
# Calculate the parameter acceptance rate
acceptrate = AR1.mcmc$accept * 100
#Print the acceptance rate as a percent. Should be higher than 10% and no greater than 60%
cat("Accept rate =", acceptrate, "%\n")

AR1.mcmc.parameters = c(median(AR1.mcmc.chain[-burnin,1]), median(AR1.mcmc.chain[-burnin,2]))
print(AR1.mcmc.parameters) # Print MCMC best estimated parameters
median.AR1.fit = model(AR1.mcmc.parameters, x)

### Check trace plots for convergence
par(mfrow=c(2,2))
plot(AR1.mcmc.chain[-burnin,1], type="l", ylab="a", xlab="Number of Runs", main="")
plot(AR1.mcmc.chain[-burnin,2], type="l", ylab="b", xlab="Number of Runs", main="")
plot(AR1.mcmc.chain[-burnin,3], type="l", ylab="Sigma", xlab="Number of Runs", main="")
plot(AR1.mcmc.chain[-burnin,4], type="l", ylab="Rho", xlab="Number of Runs", main="")

# 5:  Estimate the Hindcast
n=length(AR1.mcmc.chain[,1])
AR1.subset = AR1.mcmc.chain[seq(length(burnin),n,500),] #Only a subsequent number is neccessary (~10,000)
h = length(AR1.subset[,1])

#Calculate the simulated model hindcasts
AR1.mcmc_hind=mat.or.vec(h,datalength)
for(i in 1:datalength){
  AR1.mcmc_hind[,i]=AR1.subset[,1]+AR1.subset[,2]*x[i]
}
# Calculate the residuals with the MCMC estimates AR(1) coefficient and the sigma
res.mcmc.AR1 = mat.or.vec(h, datalength) #(nr,nc)
hindcast.AR1 = res.mcmc.AR1
for(k in 1:h) {
  for(i in 2:datalength) {
    res.mcmc.AR1[k,i] = AR1.subset[k,4]*res.mcmc.AR1[k,i-1] + rnorm(1,mean=0,sd=AR1.subset[k,3]) ## For AR1
  }
}
for(i in 1:h) {
  hindcast.AR1[i,]=AR1.mcmc_hind[i,]+res.mcmc.AR1[i,] #Add the model simulations and the residuals
}

# 6:  Analyze the Results
# Save the probability density functions of the parameters
pdfa.AR1 = density(AR1.mcmc.chain[-burnin,1])
pdfb.AR1 = density(AR1.mcmc.chain[-burnin,2])
pdfsigma.AR1 = density(AR1.mcmc.chain[-burnin,3])
pdfrho.AR1 = density(AR1.mcmc.chain[-burnin,4])

#Calculate the 90% Confidence Interval
#five = quantile(hindcast,0.05) #for one data point
#ninef = quantile(hindcast,0.95)
five.AR1 <-
  ninetyfive.AR1 <-rep(NA,datalength) # For multiple observations
for(i in 1:datalength){
  five.AR1[i] = quantile(hindcast.AR1[,i],0.05)
  ninetyfive.AR1[i] = quantile(hindcast.AR1[,i],0.95)
}
range_x_AR1=c(five.AR1, rev(ninetyfive.AR1)); range_y_AR1=c(x, rev(x))

####### Set up and Run MCMC for AR(1) homoskedastic observations assuming KNOWN RHO ---------------------------------
bound.lower = c(-10,-10,0)
bound.upper = c(10,10,1)
y.meas.err=rep(0,length(y.obs)) #Error for homoskedastic
known.rho = rho[2]

model.p=2
parnames=c("a","b","sigma.y")
source("ToyScripts/toy.R")
source("ToyScripts/KnownRho_homo_obs_likelihood_AR.R")

p = c(result$par[1], result$par[2], sd(res))
p0 = c(1,0.3,0.6)
p0 = optim(p0, function(p) -log.post(p))$par
print(round(p0,4))
library(mcmc)

step = c(0.1,0.001,0.01)
NI = 5E6 #number of iterations
burnin = seq(1,0.01*NI,1)

rho.mcmc = metrop(log.post, p0, nbatch=NI, scale=step)
rho.mcmc.chain = rho.mcmc$batch
rho.mcmc$accept
# Calculate the parameter acceptance rate
acceptrate = rho.mcmc$accept * 100
cat("Accept rate =", acceptrate, "%\n")

rho.mcmc.parameters = c(median(rho.mcmc.chain[-burnin,1]), median(rho.mcmc.chain[-burnin,2]))
print(rho.mcmc.parameters) # Print MCMC best estimated parameters
median.rho.fit = model(rho.mcmc.parameters, x)

### Check trace plots for convergence
par(mfrow=c(2,2))
plot(rho.mcmc.chain[-burnin,1], type="l", ylab="a", xlab="Number of Runs", main="")
plot(rho.mcmc.chain[-burnin,2], type="l", ylab="b", xlab="Number of Runs", main="")
plot(rho.mcmc.chain[-burnin,3], type="l", ylab="Sigma", xlab="Number of Runs", main="")
# plot(IID.mcmc.chain[-burnin,4], type="l", ylab="Rho", xlab="Number of Runs", main="")

# 8:  Estimate the Hindcast
n=length(rho.mcmc.chain[,1])
rho.subset = rho.mcmc.chain[seq(length(burnin),n,500),] #Only a subsequent number is neccessary (~10,000)
h = length(rho.subset[,1])

#Calculate the simulated model hindcasts
rho.mcmc_hind=mat.or.vec(h,datalength)
for(i in 1:datalength){
  rho.mcmc_hind[,i]=rho.subset[,1]+rho.subset[,2]*x[i]
}
# Calculate the residuals with the MCMC estimates AR(1) coefficient and the sigma
res.mcmc.rho = mat.or.vec(h, datalength) #(nr,nc)
rho.hindcast = res.mcmc.rho
for(k in 1:h) {
  for(i in 2:datalength) {
    res.mcmc.rho[k,i] = known.rho*res.mcmc.rho[k,i-1] + rnorm(1,mean=0,sd=rho.subset[k,3]) ## For AR1
  }
}
#Superimpose residuals on the hindcast fits #
for(i in 1:h) {
  rho.hindcast[i,]=rho.mcmc_hind[i,]+res.mcmc.rho[i,] #Add the model simulations and the residuals
}

# 9:  Analyze the Results
# Save the probability density functions of the parameters
pdfa.rho = density(rho.mcmc.chain[-burnin,1])
pdfb.rho = density(rho.mcmc.chain[-burnin,2])
pdfsigma.rho = density(rho.mcmc.chain[-burnin,3])

five.rho <-
  ninetyfive.rho <-rep(NA,datalength) # For multiple observations
for(i in 1:datalength){
  five.rho[i] = quantile(rho.hindcast[,i],0.05)
  ninetyfive.rho[i] = quantile(rho.hindcast[,i],0.95)
}
range_x_rho=c(five.rho, rev(ninetyfive.rho)); range_y_rho=c(x, rev(x))

######### Set up and Run MCMC assuming IID and hence uncorrelated data similiar to white noise-------------------------------
bound.lower = c(-10,-10,0)
bound.upper = c(10,10,1)
y.meas.err=rep(0,length(y.obs)) #Error for homoskedastic
model.p=2
parnames=c("a","b","sigma.y")
source("ToyScripts/toy.R")
source("ToyScripts/iid_obs_likelihood.R")

p = c(result$par[1], result$par[2], sd(res))
p0 = c(1,0.3,0.6)
p0 = optim(p0, function(p) -log.post(p))$par
print(round(p0,4))
library(mcmc)

step = c(0.1,0.001,0.01)
NI = 5E6 #number of iterations
burnin = seq(1,0.01*NI,1)

IID.mcmc = metrop(log.post, p0, nbatch=NI, scale=step)
IID.mcmc.chain = IID.mcmc$batch
IID.mcmc$accept
# Calculate the parameter acceptance rate
acceptrate = IID.mcmc$accept * 100
cat("Accept rate =", acceptrate, "%\n")

IID.mcmc.parameters = c(median(IID.mcmc.chain[-burnin,1]), median(IID.mcmc.chain[-burnin,2]))
print(IID.mcmc.parameters) # Print MCMC best estimated parameters
median.IID.fit = model(IID.mcmc.parameters, x)

### Check trace plots for convergence
par(mfrow=c(2,2))
plot(IID.mcmc.chain[-burnin,1], type="l", ylab="a", xlab="Number of Runs", main="")
plot(IID.mcmc.chain[-burnin,2], type="l", ylab="b", xlab="Number of Runs", main="")
plot(IID.mcmc.chain[-burnin,3], type="l", ylab="Sigma", xlab="Number of Runs", main="")
# plot(IID.mcmc.chain[-burnin,4], type="l", ylab="Rho", xlab="Number of Runs", main="")

# 11:  Estimate the Hindcast
n=length(IID.mcmc.chain[,1])
IID.subset = IID.mcmc.chain[seq(length(burnin),n,500),] #Only a subsequent number is neccessary (~10,000)
h = length(IID.subset[,1])

IID.mcmc_hind=mat.or.vec(h,datalength) #More than 1 observation
for(i in 1:datalength){
  IID.mcmc_hind[,i]=IID.subset[,1]+IID.subset[,2]*x[i]
}

#Superimpose noise on the hindcast fits #
IID.hindcast=IID.mcmc_hind+rnorm(datalength,mean=0,sd=IID.subset[,3]) ## For IID

# 12:  Analyze the Results
# Save the probability density functions of the parameters
pdfa.IID = density(IID.mcmc.chain[-burnin,1])
pdfb.IID = density(IID.mcmc.chain[-burnin,2])
pdfsigma.IID = density(IID.mcmc.chain[-burnin,3])

five.IID <-
  ninetyfive.IID <-rep(NA,datalength) # For multiple observations
for(i in 1:datalength){
  five.IID[i] = quantile(IID.hindcast[,i],0.05)
  ninetyfive.IID[i] = quantile(IID.hindcast[,i],0.95)
}
range_x_IID=c(five.IID, rev(ninetyfive.IID)); range_y_IID=c(x, rev(x))

# ------------------------------------------------------------------------- #
#Save the workspace:
# save.image(file = "Workspace/test_AR1vsIID.RData")
load("Workspace/test_AR1vsIID.RData")
# ------------------------------------------------------------------------- #
# 12: Plot the Results
library(RColorBrewer)
test.colors = brewer.pal(9, "YlGnBu")
IID.color = test.colors[3]
rho.color = test.colors[5]
AR1.color = test.colors[7]

makeTransparent<- function(somecolor, alpha=100){
  someColor = someColor
  newColor<-col2rgb(someColor)
  apply(newColor,2 ,
        function(curcoldata)
        {rgb(red=curcoldata[1],
             green=curcoldata[2],
             blue=curcoldata[3], alpha=alpha,
             maxColorValue=255)})
}
# FOR PDF
someColor = IID.color
trans_con_color = makeTransparent(someColor,150)

mm_TO_inches = function(mm){
  mm * 0.039370
}

single_column = mm_TO_inches(84)
double_column = mm_TO_inches(174)
maximum_width = mm_TO_inches(234)
column_height=2.7

# pdf(file="SuppFigures/sFig_AR1vsIID.pdf", family="Helvetica", height=column_height*2, width=double_column, pointsize=11)
setEPS()
postscript(file="SuppFigures/sFig_AR1vsIID_fixed.eps", horizontal = FALSE, onefile = FALSE, paper = "special", family="Helvetica",
         width=double_column, height=column_height*3, pointsize=13)
par(mfrow=c(3,2), mgp=c(1.5,.5,0),  mar=c(3.5,4,1,1)) # set figure dimensions
plot(x, obs, type="l",xlab="x",ylab="Observations") #, ylim=c(-2,6))
polygon(range_y_AR1, range_x_AR1, col=AR1.color, border=NA)
polygon(range_y_rho, range_x_rho, col=rho.color, border=NA)
polygon(range_y_IID, range_x_IID, col=IID.color, border=NA)

points(x, obs, pch=20, cex=0.5)
put.fig.letter("a.",font=2)
legend("topleft", c("AR1 (uncertain autocorrelation)", "AR1 (known autocorrelation)", "IID (assume non-correlated)","Observations",
                    "True parameter"), cex=0.75, bty="n",lty=c(1,1,1,NA,2),lwd=2, pch=c(NA,NA,NA,20,NA), 
       col=c(AR1.color,rho.color,IID.color,"black","black"))

plot(pdfa.IID, xlab="a", ylab="Probability density", main="",
     xlim=c(-1,2), lwd=2, yaxt="n", col=IID.color) #iid
lines(pdfa.AR1, col=AR1.color,lwd=2)
lines(pdfa.rho, col=rho.color,lwd=2)
abline(v=true_par[1], lty=2, lwd=2)
put.fig.letter("b.",font=2)

plot(pdfb.IID, xlab="b", ylab="Probability density", main="",lwd=2,
     xlim=c(0.19,0.21), yaxt="n", col=IID.color)
lines(pdfb.AR1, col=AR1.color,lwd=2)
lines(pdfb.rho, col=rho.color,lwd=2)
abline(v=true_par[2], lty=2, lwd=2)
put.fig.letter("c.",font=2)

plot(pdfrho.AR1, xlab=expression(paste(rho, " (Residual autocorrelation coeffient)")), 
     ylab="Probability density", main="",lwd=2,
     xlim=c(-0.99,0.99), yaxt="n", col=AR1.color) 
abline(v=known.rho, col=rho.color, lwd=2)
abline(v=0, col=IID.color, lwd=2)
abline(v=0.3, lty=2, lwd=2)
put.fig.letter("d.",font=2)

plot(pdfsigma.AR1, xlab=expression(paste(sigma, " (standard error)")), 
     ylab="Probability density", main="",lwd=2,
     xlim=c(0,1), yaxt="n", col=AR1.color) 
lines(pdfsigma.IID, col=IID.color,lwd=2)
lines(pdfsigma.rho, col=rho.color,lwd=2)
abline(v=sqrt(sigma.ar1), lty=2, lwd=2)
put.fig.letter("e.",font=2)
dev.off()

####################################### END ################################################################
