#################################################################################
#
#  -file = "Rcali_homo_model_AR.R"   Code written August 2014
#  - Author: Kelsey Ruckert (klr324@psu.edu)
#
#  -This program runs a Markov Chain Monte Carlo analysis of global sea-level
#       assuming homoskedastic errors as described in Ruckert et al. (2016).
#       For further description and references, please read the paper.
#
# THIS CODE IS PROVIDED AS-IS WITH NO WARRANTY (NEITHER EXPLICIT
# NOT IMPLICIT).  I SHARE THIS CODE IN HOPES THAT IT IS USEFUL,
# BUT I AM NOT LIABLE FOR THE BEHAVIOR OF THIS CODE IN YOUR OWN
# APPLICATION.  YOU ARE FREE TO SHARE THIS CODE SO LONG AS THE
# AUTHOR(S) AND VERSION HISTORY REMAIN INTACT.
#
#   -NOTE: This program applies the Markov Chain Monte Carlo method assuming
#       homoskedastic errors to the Rahmstorf (2007) semi-empirical sea-level model:
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
hindcast_length = 122 # there are 122 years from 1880 to 2002
projection_length = 421 # from 1880 to 2300

#------------------------------ Find Initial Parameter & Initial Hindcast ------------------------
#parm = c(.34, -0.5) #cm/year per C #the sensitivity of SLR to temperature change
#Ti=-0.5 #baseline temp (C) at which sea level is zero
timestep=1  # timesteps are 1 year
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
#plot(year, slr, pch=20, ylab="Sea-level Anomaly [mm]", xlab="Year")
#lines(year, slr.est$sle, col="blue", lwd=2)

#------------------------ Calculate the Residuals  & AR(1) Coefficient --------------------------
#Calculate Residuals
res=slr-slr.est$sle
nyears.obs=length(year) #number of years in observational time series

### Estimate and save the lag-1 autocorrelation coefficient (rho[2])
#pdf(file="Rboot3.pdf", family="Helvetica", pointsize=11, height=4.5, width=4.5)
rho=rep(NA,3)
ac=acf(res, lag.max=5, plot=TRUE, main="")# apply auto-correlation to determine correlation coefficients
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
y.meas.err=rep(0, length(slr)) # measurement error, this is constant over time making it homoskedastic

# Step 2: Define number of model parameters
model.p=3
parnames=c("alpha","base temp","initialvalue", "sigma.y", "phi11")
# Step 3: Source the physical model and statistical model
source("Scripts/sealevel_rahm_model.R")
source("Scripts/Robs_likelihood_AR.R")

# Step 4: Set up the initial parameters from the DEoptim best guess parameters
p = c(outDEoptim$optim$bestmem[1], outDEoptim$optim$bestmem[2], 
      outDEoptim$optim$bestmem[3], sd(res), rho[2]) 
p0 = c(0.34,-0.5,slr[1],0.6, 0.5) #Rahmstorf estimated best guess parameters
p0 = optim(p0, function(p) -log.post(p))$par
print(round(p0,4))
#install.packages("mcmc")
library(mcmc)

# Step 5: Set up the step size, burnin, and number of iterations to run
step = c(0.01,0.01,0.1,0.01,0.01)
NI = 1e7 #number of iterations
burnin = seq(1,0.01*NI,1) # 1% burnin

#Run the MCMC chain
mcmc.out.homo = metrop(log.post, p0, nbatch=NI, scale=step)
homoskchain = mcmc.out.homo$batch
mcmc.out.homo$accept
# Calculate the parameter acceptance rate
acceptrate = mcmc.out.homo$accept * 100
#Print the acceptance rate as a percent. Should be ~ 25%
cat("Accept rate =", acceptrate, "%\n")

#-------------------------- Estimate Parameter PDFs & Best Estimates ----------------------------
# Find the probability density function for each of the estimated parameters
h.pdfa <- density(homoskchain[length(burnin):NI,1])
h.pdfTo <- density(homoskchain[length(burnin):NI,2])
h.pdfinitialval <- density(homoskchain[length(burnin):NI,3])
h.pdfsigma <- density(homoskchain[length(burnin):NI,4])
h.pdfrho <- density(homoskchain[length(burnin):NI,5])

homChainBurnin <- homoskchain[length(burnin):NI,]

# Find best estimated parameters from bootstrapping with the median value
h.new = c(median(homoskchain[-burnin,1]), median(homoskchain[-burnin,2]), median(homoskchain[-burnin,3]))
print(h.new)

to=hindcast_length #122
h.new.est = rahmfunction(h.new, hist.temp) #best fit hindcast
to = projection_length
h.proj.est = rahmfunction(h.new, rcp85) #421 #best fit projection

# Estimating with all 10 million runs is not neccasary if the chains have
# converged. A subset of every 495th number should be sufficient.
sschain = homoskchain[seq(length(burnin),NI,495),]
k = length(sschain[,1])

## To check if subset is sufficient & for convergence: 
# save.image(file = "Workspace/homo1780.RData") # seed 3
# save.image(file = "Workspace/homo1234.RData") # seed 5
# save.image(file = "Workspace/homo1014.RData") # seed 4
# save.image(file = "Workspace/homo111.RData")  # seed 1
# save.image(file = "Workspace/homo1.RData")    # seed 2 
#They should be roughly similiar
#par(mfrow=c(3,2))
#plot(density(homoskchain[1:(NI/2),1]), main="alpha", xlab="")
#lines(density(homoskchain[(NI/2):NI,1]), col="blue")
#lines(density(sschain[,1]), col="red")
#plot(density(homoskchain[1:(NI/2),2]), main="base temp", xlab="")
#lines(density(homoskchain[(NI/2):NI,2]), col="blue")
#lines(density(sschain[,2]), col="red")
#plot(density(homoskchain[1:(NI/2),3]), main="initialvalue", xlab="")
#lines(density(homoskchain[(NI/2):NI,3]), col="blue")
#lines(density(sschain[,3]), col="red")
#plot(density(homoskchain[1:(NI/2),4]), main="sigma.y", xlab="")
#lines(density(homoskchain[(NI/2):NI,4]), col="blue")
#lines(density(sschain[,4]), col="red")
#plot(density(homoskchain[1:(NI/2),5]), main="phi11", xlab="")
#lines(density(homoskchain[(NI/2):NI,5]), col="blue")
#lines(density(sschain[,5]), col="red")

#------------- Hindcast Sea-level Rise with Uncertainty & Find Plausible Parameters ------------------
# Calculate all possible hindcasts from the subset parameter estimates.
homo.new.rate=mat.or.vec(k, nyears.obs)
homo.mcmc.fit=mat.or.vec(k, nyears.obs) # homo.mcmc.fit is the hindcast SLR simulation without noise
source("Scripts/searate_func.R") #searate function finds the hindcast rates for SLR
par=mat.or.vec(k, 2)
for(i in 1:k) {
  to=hindcast_length
  par[i,1]=sschain[i,1] # alpha parameter
  par[i,2]=sschain[i,2] # T0 parameter
  homo.new.rate[i,] = thermal(par[i,], hist.temp)[[1]]
  homo.mcmc.fit[i,1] = sschain[i,3]  # Initial value
  for (n in from:to){
    homo.mcmc.fit[i,n]=homo.mcmc.fit[i,n-1]+homo.new.rate[i,n-1]*timestep
  }
}

### Calculate hindcast residuals with the lag-1 autocorrelation coefficient estimates: ssprechain1[n,5]
### and the standard deviation (sigma) estimates: ssprechain1[n,4]
h.res.mcmc_hind=mat.or.vec(k, nyears.obs) #(nr,nc)
h.slr.mcmc_hind=h.res.mcmc_hind
for(n in 1:k) {
  for(i in 2:nyears.obs) {
    h.res.mcmc_hind[n,i] = sschain[n,5]*h.res.mcmc_hind[n,i-1] + 
      rnorm(1,mean=0,sd=sschain[n,4]) # add in the AR(1) noise
  }
}
### superimpose residuals on the hindcasts from the SLR model ###
for(i in 1:k) {
  h.slr.mcmc_hind[i,]=homo.mcmc.fit[i,]+h.res.mcmc_hind[i,]
}

#----------------------------- Project Sea-level Rise with Uncertainty --------------------------------
### Project SLR with uncertainty using the parameters generated from MCMC assuming homoskedastic
### errors and RCP8.5 temp. emission
years.mod=(alltime) # all time represent the years from 1880 to 2300
nyears.mod=length(years.mod)
homo.pred.rate=mat.or.vec(k, nyears.mod) #(nr,nc)
h.fit.mcmc_proj=mat.or.vec(k, nyears.mod) #(nr,nc) # h.fit.mcmc_proj is SLR projections without noise
for(n in 1:k) {
  to=projection_length #421
  source("Scripts/searate_func.R") #searate function finds the projected rates for SLR to 2300
  homo.pred.rate[n,] = thermal(par[n,], rcp85)[[1]]
  h.fit.mcmc_proj[n,1] = sschain[n,3] # initial value in 1880
  for (i in from:to){
    h.fit.mcmc_proj[n,i]=h.fit.mcmc_proj[n,i-1]+homo.pred.rate[n,i-1]*timestep
  }
}

###Calculate projection residuals with the lag-1 autocorrelation coefficient estimates: ssprechain1[n,5]
### and the standard deviation (sigma) estimates: ssprechain1[n,4]
h.res.mcmc_proj=mat.or.vec(k, nyears.mod) #(nr,nc)
h.slr.mcmc_proj=h.res.mcmc_proj
for(n in 1:k) {
  for(i in 2:nyears.mod) {
    h.res.mcmc_proj[n,i] = sschain[n,5]*h.res.mcmc_proj[n,i-1] + 
      rnorm(1,mean=0,sd=sschain[n,4])
  }
}
### superimpose residuals on the projections from the SLR model ###
for(i in 1:k) {
  h.slr.mcmc_proj[i,]=h.fit.mcmc_proj[i,]+h.res.mcmc_proj[i,]
}

#--------------------- Estimate PDF, CDF, and SF of SLR in 2100 & 2050 --------------------------
# Source survival function Function"
source("Scripts/plot_sf.r")

# Find the probability density function of sea-level estimates in 2100
h.prob_proj2100=mat.or.vec(k,1)
h.prob_proj2100=h.slr.mcmc_proj[,221] #The year 2100 is the 223 number in the sequence
h.pdf2100 <- density(h.prob_proj2100/100)

h.cdf2100 = ecdf(h.prob_proj2100/100) # Find the cumulative density function of SLR 2100
h.survival2050 <- plot.sf(h.slr.mcmc_proj[,171]/100, make.plot=F)  # Finds the survival function

# Find the probability density function of sea-level estimates in 2050
h.prob_proj2050=mat.or.vec(k,1)
h.prob_proj2050=h.slr.mcmc_proj[,171] #The year 2050 is the 171 number in the sequence
h.pdf2050 <- density(h.prob_proj2050/100)

h.cdf2050 = ecdf(h.prob_proj2050/100) # Find the cumulative density function of SLR 2050
h.survival2100 <- plot.sf(h.slr.mcmc_proj[,221]/100, make.plot=F) # Finds the survival function

#save.image(file = "mega_R_methods_workspace_homcon.RData")
######################################### END ###################################################
