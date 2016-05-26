#################################################################################
#
#  -file = "Rcali_heter_model_AR.R"   Code written August 2014
#  - Author: Kelsey Ruckert (klr324@psu.edu)
#
#  -This program runs a Markov Chain Monte Carlo analysis of global sea-level
#       assuming heteroskedastic errors as described in Ruckert et al. (2016).
#       For further description and references, please read the paper.
#
# THIS CODE IS PROVIDED AS-IS WITH NO WARRANTY (NEITHER EXPLICIT
# NOT IMPLICIT).  I SHARE THIS CODE IN HOPES THAT IT IS USEFUL,
# BUT I AM NOT LIABLE FOR THE BEHAVIOR OF THIS CODE IN YOUR OWN
# APPLICATION.  YOU ARE FREE TO SHARE THIS CODE SO LONG AS THE
# AUTHOR(S) AND VERSION HISTORY REMAIN INTACT.
#
#   -NOTE: This program applies the Markov Chain Monte Carlo method assuming
#       heteroskedastic errors to the Rahmstorf (2007) semi-empirical sea-level model:
#
#       -preserving the AR(1) structure of the observations
#       -creating random realizations of "natural variability"
#       -projecting new realizations to 2300
#
#----------------------------- DATA DETAILS ------------------------------------
#  -Model can be found in Rahmstorf_Science_2007
#
#  -RCP8.5 is used to create temperature simulations to 2300
#  -RCP8.5 simulates "Business as usual" and is similar to the A2 scenario
#  -The RCP8.5 temperatures are from the CNRM-CM5 model
#       (Centre National de Recherches Meteorologiques)
#  -These RCP8.5 temperatures closly resemble historical temperatures
#       from the National Oceanic and Atmospheric Administration (NOAA)
#
#  -Annual global land and ocean temperature anomalies (C)
#  -Anomalies with respect to the 20th century average
#  -http://www.ncdc.noaa.gov/cag/time-series/global
#
#  -Tide guage data from Church & White_GRL_2006
#  -http://www.psmsl.org/products/reconstructions/church.php
#
#  -This program will be sourced into the file: "Mega_Rahmstorf.R"
#
###################################################################################

#rm(list =ls()) #Clear all previous variables
#library(ash)
#library(aqfig)
#library(DEoptim)
#library(compiler)
#enableJIT(3)

# Set the seed
# set.seed(1780) #seed #3
# set.seed(1)    #seed #2
set.seed(111)  #seed #1
# set.seed(1014) #seed #4
# set.seed(1234) #seed #5

# Read in the sea-level rise observations, observation errors, years, and historic and emission temps.
source("Data/temp_sea_2300.R")
hindcast_length=122 # there are 122 years from 1880 to 2002
projection_length=421 # from 1880 to 2300

#------------------------------ Find Initial Parameter & Initial Hindcast ------------------------
#parm = c(.34, -0.5) #cm/year per C #the sensitivity of SLR to temperature change
#Ti=-0.5 #baseline temp (C) at which sea level is zero
timestep=1 # timesteps are 1 year
from=2 # start from the second year since the first year is an uncertain parameter
to=hindcast_length #122

#Run DEoptim in R to find good initial parameters
source("Scripts/Deoptim_rahm_model.R") # source the model
source("Scripts/minimize_residuals.R") # find the minimum residuals
lower=c(0,-3,err_neg[1])
upper=c(2,2,err_pos[1])
iter=1000  # specify number of iterations
outDEoptim <- DEoptim(min_res, lower, upper, 
                      DEoptim.control(itermax=iter,
                                      trace=FALSE))
print(outDEoptim$optim$bestmem)# find best initial parameters
parms = c(outDEoptim$optim$bestmem[1], outDEoptim$optim$bestmem[2], outDEoptim$optim$bestmem[3])

#Run the model with the initial parameters to create a simulation of the observations
source("Scripts/sealevel_rahm_model.R") #sealevel_rahm_model.R is the model equation
slr.est = rahmfunction(parms, hist.temp)
to=projection_length  #number of years in the projection
proj.est = rahmfunction(parms, rcp85)
to=hindcast_length

#Check to make sure true parameters fit the observations
#plot(year, slr/100, pch=20, ylab="Sea-level Anomaly [m]", xlab="Year")
#lines(year, slr.est$sle/100, col="blue", lwd=2)

#------------------------ Calculate the Residuals  & AR(1) Coefficient --------------------------
#Calculate Residuals
res=slr-slr.est$sle
nyears.obs=length(year) #number of years in observational time series

### Estimate and save the lag-1 autocorrelation coefficient (rho[2])
#pdf(file="Rheter19.pdf", family="Helvetica", pointsize=11, height=4.5, width=4.5)
rho=rep(NA,3)
ac=acf(res, lag.max=5, plot=TRUE, main="")# apply  auto-correlation to determine correlation coefficients
rho[1]=ac$acf[1]
rho[2]=ac$acf[2]
rho[3]=ac$acf[3]
rho[4]=ac$acf[4]
rho[5]=ac$acf[5]
#dev.off()

#--------------------------------- Run MCMC Calibration ----------------------------------------
# Step 1: Set up the prior ranges for the MCMC
bound.lower = c(0,-3,err_neg[1],0,-0.99)
bound.upper = c(2,2,err_pos[1],1,0.99)
y.meas.err=err.obs # measurement error, this changes over time making it heteroskedastic

# Step 2: Define number of model parameters
model.p=3
parnames=c("alpha","base temp","initialvalue", "sigma.y", "phi11")
# Step 3: Source the physical model and statistical model
source("Scripts/sealevel_rahm_model.R")
source("Scripts/Robs_likelihood_AR.R")

# Step 4: Set up the initial parameters from the DEoptim best guess parameters
p = c(outDEoptim$optim$bestmem[1], outDEoptim$optim$bestmem[2], 
      outDEoptim$optim$bestmem[3], sd(res), rho[2]) 
p0 = c(0.34,-0.5,slr[1],0.6, 0.5) # Rahmstorf estimated best guess parameters
p0 = optim(p0, function(p) -log.post(p))$par
print(round(p0,4))
library(mcmc)

# Step 5: Set up the step size, burnin, and number of iterations to run
step = c(0.02,0.02,0.1,0.01,0.01)
NI = 2.5e7 #number of iterations
burnin = seq(1,0.01*NI,1) # 1% burnin

#Run the MCMC chain
mcmc.out1 = metrop(log.post, p0, nbatch=NI, scale=step)
prechain1 = mcmc.out1$batch
mcmc.out1$accept
# Calculate the parameter acceptance rate
acceptrate = mcmc.out1$accept * 100
#Print the acceptance rate as a percent. Should be ~ 25%
cat("Accept rate =", acceptrate, "%\n")

#-------------------------- Estimate Parameter PDFs & Best Estimates ----------------------------
# Find the probability density function for each of the estimated parameters
pdfa <- density(prechain1[length(burnin):NI,1])
pdfTo <- density(prechain1[length(burnin):NI,2])
pdfinitialval <- density(prechain1[length(burnin):NI,3])
pdfsigma <- density(prechain1[length(burnin):NI,4])
pdfrho <- density(prechain1[length(burnin):NI,5])

hetChainBurnin <- prechain1[length(burnin):NI,]

# Find median estimated parameters from  with the median value
new = c(median(prechain1[-burnin,1]), median(prechain1[-burnin,2]), median(prechain1[-burnin,3]))
print(new)

to=hindcast_length #122
new.est = rahmfunction(new, hist.temp) #best fit hindcast

to = projection_length #421 #best fit projection
new.proj.est = rahmfunction(new, rcp85) #~4C increase from 1990-2100

# Estimating with all 20 million runs is not neccasary if the chains have
# converged. So we will take a subset of every 1237th number in the data set from
# the burnin to the 20 millionth run
ssprechain1 = prechain1[seq(length(burnin),NI,1237),]
h = length(ssprechain1[,1]) #length of the subset ~20,000

## To check if subset is sufficient & for convergence: 
# save.image(file = "Workspace/heter1780.RData") # seed 3
# save.image(file = "Workspace/heter1234.RData") # seed 5
# save.image(file = "Workspace/heter1014.RData") # seed 4
# save.image(file = "Workspace/heter111.RData")  # seed 1
# save.image(file = "Workspace/heter1.RData")    # seed 2
#They should be roughly similiar. The range needs to be the same
#par(mfrow=c(3,2))
#plot(density(prechain1[length(burnin):(NI/2),1]), main="alpha", xlab="")
#lines(density(prechain1[(NI/2):NI,1]), col="blue")
#lines(density(ssprechain1[,1]), col="red")
#plot(density(prechain1[1:(NI/2),2]), main="base temp", xlab="")
#lines(density(prechain1[(NI/2):NI,2]), col="blue")
#lines(density(ssprechain1[,2]), col="red")
#plot(density(prechain1[length(burnin):(NI/2),3]), main="initialvalue", xlab="")
#lines(density(prechain1[(NI/2):NI,3]), col="blue")
#lines(density(ssprechain1[,3]), col="red")
#plot(density(prechain1[1:(NI/2),4]), main="sigma.y", xlab="")
#lines(density(prechain1[(NI/2):NI,4]), col="blue")
#lines(density(ssprechain1[,4]), col="red")
#plot(density(prechain1[1:(NI/2),5]), main="phi11", xlab="")
#lines(density(prechain1[(NI/2):NI,5]), col="blue")
#lines(density(ssprechain1[,5]), col="red")

#------------- Hindcast Sea-level Rise with Uncertainty & Find Plausible Parameters ------------------
# Calculate all possible hindcasts from the subset parameter estimates.
new.rate=mat.or.vec(h, nyears.obs)
mcmc.fit=mat.or.vec(h, nyears.obs) # mcmc.fit is the hindcast SLR simulation without noise
par=mat.or.vec(h, 2)
source("Scripts/searate_func.R") #searate function finds the hindcast rates for SLR
for(i in 1:h) {
    to=hindcast_length
    par[i,1]=ssprechain1[i,1] # alpha parameter
    par[i,2]=ssprechain1[i,2] # T0 parameter
    new.rate[i,] = thermal(par[i,], hist.temp)[[1]]
    mcmc.fit[i,1] = ssprechain1[i,3]  # Initial value
    for (n in from:to){
        mcmc.fit[i,n]=mcmc.fit[i,n-1]+new.rate[i,n-1]*timestep
    }
}

### Calculate hindcast residuals with the lag-1 autocorrelation coefficient estimates: ssprechain1[n,5]
### and the standard deviation (sigma) estimates: ssprechain1[n,4]
res.mcmc_hind=mat.or.vec(h, nyears.obs) #(nr,nc)
slr.mcmc_hind=res.mcmc_hind
for(n in 1:h) {
    for(i in 2:nyears.obs) {
        res.mcmc_hind[n,i] = ssprechain1[n,5]*res.mcmc_hind[n,i-1] + 
          rnorm(1,mean=0,sd=ssprechain1[n,4]) # add in the AR(1) noise
    }
}
### superimpose residuals on the hindcasts from the SLR model ###
for(i in 1:h) {
    slr.mcmc_hind[i,]=mcmc.fit[i,]+res.mcmc_hind[i,]
}

#----------------------------- Project Sea-level Rise with Uncertainty --------------------------------
### Project SLR with uncertainty using the parameters generated from MCMC assuming heteroskedastic
### errors and RCP8.5 temp. emission
years.mod=(alltime) # all time represent the years from 1880 to 2300
nyears.mod=length(years.mod)
pred.rate=mat.or.vec(h, nyears.mod) #(nr,nc)
fit.mcmc_proj=mat.or.vec(h, nyears.mod) #(nr,nc) # fit.mcmc_proj is SLR projections without noise
for(n in 1:h) {
    to=projection_length #421
    source("Scripts/searate_func.R") #searate function finds the projected rates for SLR to 2300
    pred.rate[n,] = thermal(par[n,], rcp85)[[1]]
    fit.mcmc_proj[n,1] = ssprechain1[n,3] # initial value in 1880
    for (i in from:to){
        fit.mcmc_proj[n,i]=fit.mcmc_proj[n,i-1]+pred.rate[n,i-1]*timestep
    }
}

###Calculate projection residuals with the lag-1 autocorrelation coefficient estimates: ssprechain1[n,5]
### and the standard deviation (sigma) estimates: ssprechain1[n,4]
res.mcmc_proj=mat.or.vec(h, nyears.mod) #(nr,nc)
slr.mcmc_proj=res.mcmc_proj
for(n in 1:h) {
    for(i in 2:nyears.mod) {
        res.mcmc_proj[n,i] = ssprechain1[n,5]*res.mcmc_proj[n,i-1] + 
          rnorm(1,mean=0,sd=ssprechain1[n,4])
    }
}
### superimpose residuals on the projections from the SLR model ###
for(i in 1:h) {
    slr.mcmc_proj[i,]=fit.mcmc_proj[i,]+res.mcmc_proj[i,]
}

#--------------------- Estimate PDF, CDF, and SF of SLR in 2100 & 2050 --------------------------
# Source survival function Function"
source("Scripts/plot_sf.r")

# Find the probability density function of sea-level estimates in 2100
prob_proj2100=mat.or.vec(h,1)
prob_proj2100=slr.mcmc_proj[,221] #The year 2100 is the 221 number in the sequence
pdf2100 <- density(prob_proj2100/100)

cdf2100 = ecdf(prob_proj2100/100) # Find the cumulative density function of SLR 2100
survival2050 <- plot.sf(slr.mcmc_proj[,171]/100, make.plot=F) # Finds the survival function

# Find the probability density function of sea-level estimates in 2050
prob_proj2050=mat.or.vec(h,1)
prob_proj2050=slr.mcmc_proj[,171] #The year 2050 is the 171 number in the sequence
pdf2050 <- density(prob_proj2050/100)

cdf2050 = ecdf(prob_proj2050/100) # Find the cumulative density function of SLR 2050
survival2100 <- plot.sf(slr.mcmc_proj[,221]/100, make.plot=F) # Finds the survival function

#save.image(file = "mega_R_methods_workspace_hetcon.RData")
############################################ END #############################################
