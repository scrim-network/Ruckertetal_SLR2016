#################################################################################
#
#  -file = "Toy_MCMC_Test2b.R"   Code written February 2015
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
#       -Workspace file: "test_observation_num_hetData.RData"
#
#  -The spreadsheet: "testhypoth_surprise.csv" is a previously calculated surprise index for
#       all the following tests.
#
###################################################################################
rm(list =ls()) #Clear global environment
#install.packages('mvtnorm')
library(mvtnorm)
library(compiler)
enableJIT(3)

source("Scripts/put_fig_letter.R")
#----------------------- Test #1: Length of Observations ----------------------------#
#load("Workspace/testiidall.RData") #Load results from test #1 and test 2a

#------------ Test #2: Adding in Autocorrelation & Test 3: Heteroskedastic Errors -------------#
################# Step 1: Set up the true model & Generate generate AR 1 errors
################# Step 1a: Choose the number of observations
#datalength = 2
#datalength = 5
#datalength = 20
#datalength = 50
#datalength = 100
#datalength = 150
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

################# Step 2: Generate time varying errors for heteroskedastic assumption
#Generate varying errors for Test #3
# NOTE: 2.9-0.9 are arbitruary numbers
err.1 = sigma.inov/5
err.end = err.1*0.25

errors = seq(err.1,err.end, length.out=datalength)
err_pos = obs+(1.65*errors) # For ploting purposes 90% confidence level
err_neg = obs-(1.65*errors)
########################
#generate hetero error:
H.err = rnorm(datalength, mean=0, sd=errors)
H.obs = true.obs$mod.obs + t.s + H.err #Add the observations and the residuals

################# Step 3: Find initial parameter guesses
fn = function(parameters, x, H.obs){
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

result = optim(true_par, fn, gr=NULL, x, H.obs)

################# Step 4:  Simulate the observations
source("ToyScripts/toy.R")
parameter = c(result$par[1],result$par[2]) #Use the initial parameter guesses
H_y.obs= model(parameter,x)

#Plot the observations
#plot(x,obs, pch=20, xlab="x", ylab="Data")
#lines(x, y.obs$mod.obs, col="red", lwd=3)
#lines(x,err_pos, col="slateblue", lwd=3)
#lines(x,err_neg, col="slateblue", lwd=3)

res = H.obs-H_y.obs$mod.obs #Estimate the residuals

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
## AND
# Make sure the observations are set to use the heteroskedastic observations (H.obs)
##-------------------------------------------------------------------##
######################### Set up and Run MCMC #######################
bound.lower = c(-10,-10,0,-0.99)
bound.upper = c(10,10,1,0.99)
y.meas.err=errors # heteroskedastic errors
model.p=2
parnames=c("a","b","sigma.y","phi11")
source("ToyScripts/toy.R")
source("ToyScripts/toy_obs_likelihood_AR.R")

p = c(result$par[1], result$par[2], sd(res), rho[2])
p0 = c(1,0.3,0.6, 0.5)
p0 = optim(p0, function(p) -log.post(p))$par
print(round(p0,4))
library(mcmc)

step = c(0.1,0.001,0.01,0.01) # for datasets larger than 20 observations
#step = c(0.2,0.1,0.1,0.1) # for datasets equal to or less than 20 observations
NI = 1E6 #number of iterations
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
ssprechain1 = prechain1[seq(length(burnin),n,1000),] #Only a subsequent number is neccessary (~10,000)
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
pdfsig = density(prechain1[,3])
pdfrho = density(prechain1[,4])

#Calculate the Surprise index by observation the ratio of outliers to the expected # of outliers
# The confidence intervals in the paper are:
# Percentages: 10%, 20%, 30%, 40%, 50%, 60%, 70%, 80%, 90%, 92%, 95%,  96%, 97%,  98%, 99%,  100%
# low numbers: 0.45,0.40,0.35,0.30,0.25,0.20,0.15,0.10,0.05,0.04,0.025,0.02,0.015,0.01,0.005,0
# high numbers:0.55,0.60,0.65,0.70,0.75,0.80,0.85,0.90,0.95,0.96,0.975,0.98,0.985,0.99,0.995,1
five <-
  ninetyfive <-rep(NA,datalength) # For multiple observations
for(i in 1:datalength){
  five[i] = quantile(hindcast[,i],0.05)
  ninetyfive[i] = quantile(hindcast[,i],0.95)
}
range_x=c(five, rev(ninetyfive)); range_y=c(x, rev(x))

#Calculate the Surprise index by observation the ratio of outliers to the expected # of outliers
par(mfrow=c(1,1), mgp=c(1.5,.5,0), mar=c(4, 4, 3, 1))
plot(x, H.obs, type="l",xlab="x",ylab="Observations")#, ylim=c(-5,20))
points(x, H.obs, pch=20)
polygon(range_y, range_x, col="red")
points(x, H.obs, pch=20)
lines(x,err_pos, col="slateblue", lwd=3)
lines(x,err_neg, col="slateblue", lwd=3)
legend("topleft", c("90% CI","Observations"), pch=c(15,20),bty="n", col=c("red","black"))

################# Step 8:  Save the Heteroskedastic Results (for 200 observations for test 2)
# Save the pdf of the parameters
hettwoh_a = density(prechain1[,1]) ; hettwoh_b = density(prechain1[,2])
#Save the 90% confidence interval
hettwoh_x = range_x ; hettwoh_y = range_y

#Again save the workspace: (for 200 observations for test 2)
#save.image(file = "Workspace/testiidall_hettest.RData")

################# Step 8:  Save the varying number of observation results (test 1)
#After running each analysis varying the observations save the pdf

two_a = density(prechain1[,1]) ; two_b = density(prechain1[,2]); two_s = density(prechain1[,3]) ; two_r = density(prechain1[,4])
five_a = density(prechain1[,1]) ; five_b = density(prechain1[,2]); five_s = density(prechain1[,3]) ; five_r = density(prechain1[,4])
twent_a = density(prechain1[,1]) ; twent_b = density(prechain1[,2]); twent_s = density(prechain1[,3]) ; twent_r = density(prechain1[,4])
fity_a = density(prechain1[,1]) ; fity_b = density(prechain1[,2]); fity_s = density(prechain1[,3]) ; fity_r = density(prechain1[,4])
hund_a = density(prechain1[,1]) ; hund_b = density(prechain1[,2]); hund_s = density(prechain1[,3]) ; hund_r = density(prechain1[,4])
hund50_a = density(prechain1[,1]) ; hund50_b = density(prechain1[,2]); hund50_s = density(prechain1[,3]) ; hund50_r = density(prechain1[,4])
twoh_a = density(prechain1[,1]) ; twoh_b = density(prechain1[,2]); twoh_s = density(prechain1[,3]) ; twoh_r = density(prechain1[,4])


#After running each analysis varying the observations save the 90% Confidence Interval
two_x = range_x ; two_y = range_y
five_x = range_x ; five_y = range_y
twent_x = range_x ; twent_y = range_y
fity_x = range_x ; fity_y = range_y
hund_x = range_x ; hund_y = range_y
hund50_x = range_x ; hund50_y = range_y
twoh_x = range_x ; twoh_y = range_y

#Save the workspace:
save.image(file = "Workspace/test_observation_num_hetData.RData")
#load("Workspace/test_observation_num_hetData.RData")

#------------------------------- Sup. Figure 7 -----------------------------------#
# Compare the surprise index and the pdfs for when the number of observations changes:
# pdf(file="SuppFigures/nRuckertetal_sup8.pdf", family="Helvetica", height=5.4, width=6.7,pointsize=11)
# par(mfrow=c(2,2),mgp=c(1.5,.5,0), mar=c(3.5,4,4,1))
library(RColorBrewer)
test.colors = brewer.pal(7, "BrBG")
test.colors[4] = "gray"
circlesize = seq(from=1, to=4, length.out=7)

setEPS()
postscript(file="SuppFigures/sfigure7a_7d.eps", horizontal = FALSE, onefile = FALSE, paper = "special", family="Helvetica", width=6.7, height=5.4, pointsize=11)
par(mfrow=c(2,2), mgp=c(1.5,.5,0),  mar=c(3.5,4,1,1)) # set figure dimensions

# 1) PDF
plot(two_a, xlab="a", ylab="Probability density", main="",
ylim=c(0,3.2), lwd=3, yaxt="n", col=test.colors[1]) #iid
lines(five_a, col=test.colors[2],lwd=3)
lines(twoh_a, col=test.colors[7],lwd=3)
lines(hund50_a, col=test.colors[6],lwd=3)
lines(hund_a, col=test.colors[5],lwd=3)

lines(fity_a, col=test.colors[4],lwd=3)
lines(twent_a, col=test.colors[3],lwd=3)
lines(c(-10,-10),c(0,0.18), lwd=2, lty=3)
lines(c(10,10),c(0,0.18), lwd=2, lty=3)
lines(c(-10,10),c(0.18,0.18), lwd=2, lty=3)
put.fig.letter("a.",font=2)
legend("topright", c("Prior","1:1 line","2","5","20","50","100","150", "200"), cex=0.75, pt.cex=c(NA,NA,circlesize[1:7]), y.intersp = 1.5,
lty=c(3,2,1,1,1,1,1,1,1), pch=c(NA,NA,20,20,20,20,20,20,20),lwd=2, col=c("black","black",test.colors[1:7]))

# 2) Surprise index
surprise = read.csv("ToyScripts/test_num_obs_hetdata.csv", skip=1)
percent = surprise[1:16,1]
twohund = surprise[1:16,15]
one50hund = surprise[1:16,13]
onehund = surprise[1:16,11]
fifty = surprise[1:16,9]
twenty = surprise[1:16,7]
five = surprise[1:16,5]
two = surprise[1:16,3]

#par(mgp=c(1.5,.5,0), mar=c(3.5, 3, 4, 2))
plot(percent, percent, typ="l", lty=2, ylab="Percent covered [%]",
xlab="Credible interval [%]", lwd=2, ylim=c(0,100))
# points(percent, two*100, pch=20, col=test.colors[1], cex=circlesize[1])
# points(percent, five*100, pch=20, col=test.colors[2], cex=circlesize[2])
# points(percent, twenty*100, pch=20, col=test.colors[3], cex=circlesize[3])
# points(percent, fifty*100, pch=20, col=test.colors[4], cex=circlesize[4])
# points(percent, onehund*100, pch=20, col=test.colors[5], cex=circlesize[5])
# points(percent, one50hund*100, pch=20, col=test.colors[6], cex=circlesize[6])
# points(percent, twohund*100, pch=20, col=test.colors[7], cex=circlesize[7])
#-----
points(percent, twohund*100, pch=20, col=test.colors[7], cex=circlesize[7])
points(percent, one50hund*100, pch=20, col=test.colors[6], cex=circlesize[6])
points(percent, onehund*100, pch=20, col=test.colors[5], cex=circlesize[5])
points(percent, fifty*100, pch=20, col=test.colors[4], cex=circlesize[4])
points(percent, twenty*100, pch=20, col=test.colors[3], cex=circlesize[3])
points(percent, five*100, pch=20, col=test.colors[2], cex=circlesize[2])
points(percent, two*100, pch=20, col=test.colors[1], cex=circlesize[1])
put.fig.letter("c.",font=2)
lines(c(88,102),c(102,102))
lines(c(88,88),c(88,102))
lines(c(102,102),c(88,102))
lines(c(88,102),c(88,88))
text(100,92, font=2, "d.")
text(45,55, "underconfident", srt=35)
text(50,35, "overconfident", srt=35)

#par(mgp=c(1.5,.5,0), mar=c(5,4,2.5,1))
plot(two_b, xlab="b", ylab="Probability density", main="",lwd=3,
ylim=c(0,350), yaxt="n", col=test.colors[1], xlim=c(-10,10)) #iid
lines(five_b,lwd=3, col=test.colors[2])
lines(twoh_b, lwd=3, col=test.colors[7])
lines(hund50_b, lwd=3, col=test.colors[6])
lines(hund_b, lwd=3, col=test.colors[5])

lines(fity_b, lwd=3, col=test.colors[4])
lines(twent_b, lwd=3, col=test.colors[3])
lines(c(-10,-10),c(0,50), lwd=2, lty=3)
lines(c(10,10),c(0,50), lwd=2, lty=3)
lines(c(-10,10),c(50,50), lwd=2, lty=3)
put.fig.letter("b.",font=2)

#par(mgp=c(1.5,.5,0), mar=c(5, 3, 2.5, 2))
plot(percent, percent, type="l", lty=2, ylab="Percent covered [%]",
xlab="Credible interval [%]", lwd=2, xlim=c(90,100), ylim=c(90,100))
# points(percent, two*100, pch=20, col=test.colors[1], cex=circlesize[1])
# points(percent, five*100, pch=20, col=test.colors[2], cex=circlesize[2])
# points(percent, twenty*100, pch=20, col=test.colors[3], cex=circlesize[3])
# points(percent, fifty*100, pch=20, col=test.colors[4], cex=circlesize[4])
# points(percent, onehund*100, pch=20, col=test.colors[5], cex=circlesize[5])
# points(percent, one50hund*100, pch=20, col=test.colors[6], cex=circlesize[6])
# points(percent, twohund*100, pch=20, col=test.colors[7], cex=circlesize[7])
#-----
points(percent, twohund*100, pch=20, col=test.colors[7], cex=circlesize[7])
points(percent, one50hund*100, pch=20, col=test.colors[6], cex=circlesize[6])
points(percent, onehund*100, pch=20, col=test.colors[5], cex=circlesize[5])
points(percent, fifty*100, pch=20, col=test.colors[4], cex=circlesize[4])
points(percent, twenty*100, pch=20, col=test.colors[3], cex=circlesize[3])
points(percent, five*100, pch=20, col=test.colors[2], cex=circlesize[2])
points(percent, two*100, pch=20, col=test.colors[1], cex=circlesize[1])
put.fig.letter("d.",font=2)
text(93,94.5, "underconfident", srt=35)
text(94,93, "overconfident", srt=35)
# legend("bottomright", c("2","5","20","50","100","150", "200","1:1 line"),
# pch=c(20,20,20,20,20,20,20,NA), lty=c(NA,NA,NA,NA,NA,NA,NA,2), ncol=2, lwd=2, cex=c(circlesize[1:7], 1),
# col=c(test.colors[1:7],"black"))
dev.off()


################################################################################################