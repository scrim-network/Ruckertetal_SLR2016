#################################################################################
#
#  -file = "toy_bootstrap.R"   Code written February 2015
#  - Author: Kelsey Ruckert (klr324@psu.edu)
#
#  -This program analyses the bootstrap method with 200 hundred generated AR(1) residuals.
#       The bootstrap method is applied to a linear model as described in the appendix of Ruckert et al.
#        (2016). The graphs will be saved as pdf files in the current working directory.
#       For further description and references, please read the paper and appendix.
#
#   -NOTE: This program applies a Bootstrap method that:
#
#       -creates random realizations of "natural variability"
#       -Assumes Homoskedastic errors
#
# THIS CODE IS PROVIDED AS-IS WITH NO WARRANTY (NEITHER EXPLICIT
# NOT IMPLICIT).  I SHARE THIS CODE IN HOPES THAT IT IS USEFUL, 
# BUT I AM NOT LIABLE FOR THE BEHAVIOR OF THIS CODE IN YOUR OWN
# APPLICATION.  YOU ARE FREE TO SHARE THIS CODE SO LONG AS THE
# AUTHOR(S) AND VERSION HISTORY REMAIN INTACT.
#
#  -The spreadsheet: "bootstrapsurprise.csv" is a previously calculated surprise index for
#       the following test.
#
###################################################################################
rm(list =ls()) #Clear global environment
install.packages("DEoptim")
library(DEoptim)
library(compiler)
enableJIT(3)

source("../Scripts/put_fig_letter.R")

#----------------------- Step 1: Choose the number of observations ----------------------------#
datalength = 200
x = 1:datalength
################# Step 1: Set up the true model & Generate generate AR 1 errors
true_par = c(1,0.2) #a=1 b=0.2
source("ToyScripts/toy.R") # this is the model: y = a+b*x
true.obs= model(true_par,x)

rho=0.3  # define your correlation
sigma.ar1=0.3 # define your ar1 white noise inovation variance
sigma.inov = sigma.ar1/(1-0.3^2) # this is your stationary variance for the AR1 process. 
# i.e.the variance for time 0

t.s0 = rnorm(1,0,sd = sqrt(sigma.inov))
t.s<-c()
t.s[1]<-t.s0
for ( i in (2:datalength)){ # write a loop to simulate AR1 residuals
  t.s[i]<-rho*t.s[i-1]+rnorm(1,mean=0,sd=sqrt(sigma.ar1))
}

obs = true.obs$mod.obs + t.s #Add the observations and the residuals

#------------------------- Find Initial Parameters -------------------------------
#Run DEoptim in R to find good initial parameters
source("ToyScripts/de_boot.R") # source the model
source("ToyScripts/min_res.R") # find the minimum residuals
lower=c(-10,-10)  
upper=c(10,10)
iter=1000  # specify number of iterations
outDEoptim <- DEoptim(min_res, lower, upper, 
                      DEoptim.control(itermax=iter,
                                      trace=FALSE))
print(outDEoptim$optim$bestmem)# find best initial parameters
parms = c(outDEoptim$optim$bestmem[1], outDEoptim$optim$bestmem[2])

#Run the model with the initial parameters to create the best simulation of the observations
source("ToyScripts/toy.R") #linear model
y.obs= model(parms,x)
# Project the linear model to 300
plength = 300
project = 1:plength
p.obs = model(parms,project)

#Plot the observations and the best fit
pdf(file="../ToyFigures/plot1.pdf")
par(mfrow=c(1,1))
plot(x, obs, pch=20, xlab="x", ylab="Observations")
lines(x, y.obs$mod.obs, lwd=2, col="blue")
legend("topleft", c("Observations","Best fit"), pch=c(20,NA), lty=c(NA,1),
       col=c("black","blue"))
dev.off()
#---------- Calculate the Residuals ------------------
### Calculate Residuals during observed time series (data - polynomial fit)  ###
res <- obs - y.obs$mod.obs

#---------- Bootstrap the Residuals & Save AR(1) Coefficient and Sigma ---------
### We need to retain the auto-correlated structure of the residuals ###
### determine AR coefficients of the original residuals (data - polynomial fit)  ###
pdf(file="../ToyFigures/plot1a.pdf")  # write to pdf, define a pdf file to write to
rho=rep(NA,3)
ac <- acf(res, lag.max=5, plot=TRUE, main="Autocorrelated Residuals")  # apply auto-correlation to determine correlation coefficientsrho[1] <- ac$acf[ 1 ]
rho[2] <- ac$acf[ 2 ]
rho[3] <- ac$acf[ 3 ]
rho[4] <- ac$acf[ 4 ]
rho[5] <- ac$acf[ 5 ]
dev.off()

# Find the standard deviation (sigma)
N=5000  # Number of bootstrap samples
white.boot = mat.or.vec(N, datalength) # create matrix (nr,nc)
white.boot_sd = rep(NA,N)
for(i in 1:N) {
  white.boot[i,1:datalength] = sample(res,size=datalength,replace=TRUE)
  white.boot_sd[i] = sd(white.boot[i,]) #this estimtes stationary variance
}

### Calculate the bootstrapped residuals with an lag-1 autocorrelation coefficient: rho[2]
# IMPORTANT: estimate the white noice variance (sigma) from the stationary variance to avoid accounting
# for autocorrelation twice
sigma = sqrt((white.boot_sd^2)*(1-(rho[2]^2))) #what is this?
res.boot=mat.or.vec(N, datalength) #(nr,nc)
for(n in 1:N) {
  for(i in 2:datalength) {
    res.boot[n,i] = rho[2]*res.boot[n,i-1] + rnorm(1,mean=0,sd=sigma[n])
  }
}
### Superimpose residuals on the hindcasts from the linear model to make N bootstrap samples###
lin.boot=mat.or.vec(N, datalength) #(nr,nc)
for(i in 1:N) {
  lin.boot[i,]=y.obs$mod.obs+res.boot[i,]
}

#Plot the best hindcast with fits (N bootstrap samples)
pdf(file="../ToyFigures/plot2.pdf")  # write to pdf, define a pdf file to write to
plot(x, lin.boot[1,], type="l",main="Hindcast Fit + Noise",
     xlab="x",ylab="observations",col="black",lwd=1,ylim=c(0,40))
for(i in 1:N) {
  lines(x, lin.boot[i,], col="gold", lwd=1)
}
lines(x, y.obs$mod.obs, col="black", lwd=1)
points(x,obs, pch=20, col="black", lwd=2)
legend("topleft", c("Observations","Best fit + Noise","Best fit"), lwd=2, 
       lty=c(NA,1,1), pch=c(20,NA,NA), col=c("black","gold","black"))
dev.off()

###IMPORTANT: calculate polynomial coefficients for the bootstrapped samples###
########## NOTE: THIS MAY TAKE A FEW MINUTES! ###
boot.fit_coef=mat.or.vec(N, 2)
for(i in 1:N) {
  source("ToyScripts/de_boot.R") # source the model
  source("ToyScripts/bootstrapped_min_res.R") # find the minimum residuals with the Many bootstrap samples
  lower=c(-10,-10)  
  upper=c(10,10)
  iter=100  # specify number of iterations
  outDEoptim <- DEoptim(min_res, lower, upper, 
                        DEoptim.control(itermax=iter,
                                        trace=FALSE))
  boot.fit_coef[i,1]=outDEoptim$optim$bestmem[1]
  boot.fit_coef[i,2]=outDEoptim$optim$bestmem[2]
}
#---------------------------- Hincast the model with uncertainty ------------------------------------
# Calculate the smooth hindcast fits
boot.fit=mat.or.vec(N, datalength)
for(n in 1:N) {
  for(i in 1:datalength) {
    boot.fit[n,i]=boot.fit_coef[n,1] + boot.fit_coef[n,2]*x[i]
  }
}
# Add the AR(1) residual noise to the smooth hindcast fits
boot.fit_prob=mat.or.vec(N, datalength) #(nr,nc)
for(i in 1:N) {
boot.fit_prob[i,] = boot.fit[i,]+res.boot[i,]
}

#---------------------------- Project the model with uncertainty ------------------------------------
#Calculate smooth projections
boot.fit_proj=mat.or.vec(N, plength)
for(n in 1:N) {
for(i in 1:plength) {
  boot.fit_proj[n,i]=boot.fit_coef[n,1] + boot.fit_coef[n,2]*project[i]
}
}
#----------------------------------------------------------------
#Plot the smoothed projected fits
#Some observations should be outside the range of the fits
pdf(file="../ToyFigures/plot3.pdf")  # write to pdf, define a pdf file to write to
plot(project, boot.fit_proj[1,], type="l",main="Projected Smooth Fits",
     xlab="x",ylab="observations",col="black",lwd=1,ylim=c(0,60))
for(i in 1:N) {
  lines(project, boot.fit_proj[i,], col="red", lwd=1)
}
lines(project, p.obs$mod.obs, col="black", lwd=1)
points(x,obs, pch=20, col="black", lwd=2)
legend("topleft", c("Observations","Boot smoothed fits","Best fit"), lwd=2, 
       lty=c(NA,1,1), pch=c(20,NA,NA), col=c("black","red","black"))
dev.off()
#---------------------------------------------------------------------------

### calculate projected AR(1) residual noise ###
res.boot_proj=mat.or.vec(N, plength) #(nr,nc)
boot_proj=res.boot_proj
for(n in 1:N) {
  for(i in 2:plength) {
    res.boot_proj[n,i] = rho[2]*res.boot_proj[n,i-1] + rnorm(1,mean=0,sd=sigma[n])
  }
}

### superimpose residuals on polynomial smooth fits to 
### calculate Probabilistic bootstrap samples###
boot.proj=mat.or.vec(N, plength) #(nr,nc)
for(i in 1:N) {
  boot.proj[i,]=boot.fit_proj[i,]+res.boot_proj[i,]
}
#------------------------ Analyze the Results ------------------------------------
#Plot the projected fits with the added noise
# The range should encompass all of the observations
pdf(file="../ToyFigures/plot4.pdf")  # write to pdf, define a pdf file to write to
plot(project, boot.proj[1,], type="l",main="Projected Fits + Noise",
     xlab="x",ylab="observations",col="black",lwd=1,ylim=c(0,60))
for(i in 1:N) {
  lines(project, boot.proj[i,], col="red", lwd=1)
}
lines(project, p.obs$mod.obs, col="black", lwd=1)
points(x,obs, pch=20, col="black", lwd=2)
legend("topleft", c("Observations","Boot fits + noise","Best fit"), lwd=2, 
       lty=c(NA,1,1), pch=c(20,NA,NA), col=c("black","red","black"))
dev.off()

### Plot the 90% hindcast ###
# Calculate the 90% confidence interval
# And Calculate the surprise index
boot_f <-
  boot_nf <-rep(NA,plength)
for(i in 1:plength){
  boot_f[i] = quantile(boot.proj[,i],0.05)
  boot_nf[i] = quantile(boot.proj[,i],0.95)
}
boot_x=c(boot_f, rev(boot_nf)); boot_y=c(project, rev(project))

#Calculate the Surprise index by observation the ratio of outliers to the expected # of outliers
# The confidence intervals in the paper are:
# Percentages: 10%, 20%, 30%, 40%, 50%, 60%, 70%, 80%, 90%, 92%, 95%,  96%, 97%,  98%, 99%,  100%
# low numbers: 0.45,0.40,0.35,0.30,0.25,0.20,0.15,0.10,0.05,0.04,0.025,0.02,0.015,0.01,0.005,0
# high numbers:0.55,0.60,0.65,0.70,0.75,0.80,0.85,0.90,0.95,0.96,0.975,0.98,0.985,0.99,0.995,100

pdf(file="../ToyFigures/plot5.pdf")
plot(x, y.obs$mod.obs, type="l",xlab="x",ylab="Observations", 
     ylim=c(0,40), xlim=c(0,200))
#Plot the 90% confidence interval envelope
polygon(boot_y, boot_x, col="red")
points(x, obs, pch=20, cex=0.5, col="black")
legend("topleft", c("90% CI Boot Fits + noise", "Observations"), 
       pch=c(15,20),bty="n", col=c("red","black"))
dev.off()

# Plot the surpise index
### NOTE: This has been previously calculated
surprise = read.csv("ToyScripts/bootstrapsurprise.csv")
percent = surprise[1:16,1]
boot.sur = surprise[1:16,3]*100
boot.dev = surprise[1:16,5]
pdf(file="ToyFigures/plot6.pdf", family="Helvetica", width=3.5, height=8.1, pointsize=13)
par(mfrow=c(3,1), mgp=c(1.5,.5,0),mar=c(2.5,4,4,2))
plot(percent, percent, typ="l", lty=2, lwd=2, xlab="Confidence Interval [%]", 
     ylab="Percent Covered [%]")
points(percent,boot.sur, pch=8, col="lightseagreen", lwd=2)
put.fig.letter("a",font=2)
text(50,65, "underconfident", srt=35)
text(55,45, "overconfident", srt=35)

# set up box in the upper corner to show where panel b is located
lines(c(88,102),c(102,102))
lines(c(88,88),c(88,102))
lines(c(102,102),c(88,102))
lines(c(88,102),c(88,88))
text(100,92, font=2, "b")

legend("bottomright", c("Perfect 1:1 line","Bootstrap SI"), 
       pch=c(NA,8), lty=c(2,NA), lwd=2,
       col=c("black","lightseagreen"))

# zoomed plot of the section above 90%
par(mgp=c(1.5,.5,0), mar=c(3.75,4,2.75,2)) 
plot(percent, percent, type="l", lty=2, ylab="Percent covered [%]", 
     xlab="Confidence interval [%]", lwd=2, xlim=c(90,100), ylim=c(90,100))
points(percent, boot.sur, pch=8, col="lightseagreen", lwd=2)
put.fig.letter("b",font=2)
text(94,95, "underconfident", srt=35)
text(95,94, "overconfident", srt=35)

par(mgp=c(1.5,.5,0), mar=c(5,4,2.5,2))
plot(percent, boot.dev, typ="l", col="lightseagreen", lwd=2, ylim=c(-0.06,0.06),
     ylab="Surprise index [Dimensionless] ", xlab="Confidence interval [%]")
points(percent, boot.dev, pch=8, col="lightseagreen")
abline(h=0, lty=2, lwd=2)
text(35,0.03, "underconfident")
text(35,-0.04, "overconfident")

put.fig.letter("c",font=2)
dev.off()

#Calculate the density of each parameter
pdfa = density(boot.fit_coef[,1])
pdfb = density(boot.fit_coef[,2])

pdf(file="../ToyFigures/plot7.pdf", width=3.5, height=5.4)
par(mfrow=c(2,1), mgp=c(1.5,.5,0),mar=c(2.5,4,4,2))
plot(pdfa, col="blue", xlab="a", ylab="Probability Density", main="", yaxt="n")
put.fig.letter("a",font=2)

par(mgp=c(1.5,.5,0), mar=c(5,4,2.5,2))
plot(pdfb, col="red", xlab="b", ylab="Probability Density", main="", yaxt="n")
put.fig.letter("b",font=2)
dev.off()

#################### END ####################################
