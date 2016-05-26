#################################################################################
#
#  -file = "bootstrap_Rahm_new.R"   Code written March 2014, Updated July 2014
#  - Author: Kelsey Ruckert (klr324@psu.edu)
#
#  -This program runs a bootstrap analysis of global sea-level
#       as described in Ruckert et al. (2016). For further
#       description and references, please read the paper.
#
# THIS CODE IS PROVIDED AS-IS WITH NO WARRANTY (NEITHER EXPLICIT
# NOT IMPLICIT).  I SHARE THIS CODE IN HOPES THAT IT IS USEFUL,
# BUT I AM NOT LIABLE FOR THE BEHAVIOR OF THIS CODE IN YOUR OWN
# APPLICATION.  YOU ARE FREE TO SHARE THIS CODE SO LONG AS THE
# AUTHOR(S) AND VERSION HISTORY REMAIN INTACT.
#
#   -NOTE: This program applies the bootstrap method to the Rahmstorf (2007)
#       semi-empirical sea-level model:
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

#rm(list = ls()) #clear all previous variables
#install packages
#install.packages("DEoptim")
#install.packages("ash")
#install.packages("aqfig")
#library(ash)
#library(aqfig)
#library(DEoptim)

# Set the seed
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

#------------------------------ Find Initial Parameter ------------------------
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
parms = c(outDEoptim$optim$bestmem[1], outDEoptim$optim$bestmem[2],outDEoptim$optim$bestmem[3])

#Run the model with the initial parameters to create best simulation of the observations
source("Scripts/sealevel_rahm_model.R") #sealevel_rahm_model.R is the semi-emipirical model
slr.est = rahmfunction(parms, hist.temp)
to=projection_length #number of years in the projection
proj.est = rahmfunction(parms, rcp85)
to=hindcast_length

#---------- Calculate the Residuals ------------------
### Calculate Residuals during the observed time series (data - model)  ###
res=slr-slr.est$sle
nyears.obs=length(year) #number of years in observational time series

#---------- Bootstrap the Residuals & Save AR(1) Coefficient and Sigma ---------
### We need to retain the auto-correlated structure of the residuals ###
### determine AR coefficients of the original residuals (data - polynomial fit)  ###
#pdf(file="autocorboot.pdf", family="Helvetica", pointsize=11, height=4.5, width=4.5)
rho=rep(NA,3)
ac=acf(res, lag.max=5, plot=TRUE, main="")# apply auto-correlation to determine correlation coefficients
rho[1]=ac$acf[1]
rho[2]=ac$acf[2] #AR1 coefficient
rho[3]=ac$acf[3]
rho[4]=ac$acf[4]
rho[5]=ac$acf[5]
#dev.off()

# Find the standard deviation (sigma)
N=20000 #Number of bootstrap samples, twenty thousand
rahm.boot=mat.or.vec(N, nyears.obs) #create matrix (# of rows, # of columns)
sigma.boot_sd=rep(NA,N)
for(i in 1:N){
  rahm.boot[i,1:nyears.obs]=sample(res, size=nyears.obs, replace=TRUE)
  sigma.boot_sd[i] = sd(rahm.boot[i,]) # Estimates the stationary variance
}
# IMPORTANT: estimate the white noise variance from the stationary vaiance to avoid accounting
# for autocorrelation twice
sigma = sqrt((sigma.boot_sd^2)*(1-(rho[2]^2)))

### Calculate the bootstrapped residuals with an lag-1 autocorrelation coefficient: rho[2]
res.boot_hind=mat.or.vec(N, nyears.obs) #(nr,nc)
slr.boot_hind=res.boot_hind
for(n in 1:N) {
    for(i in 2:nyears.obs) {
        res.boot_hind[n,i] = rho[2]*res.boot_hind[n,i-1] + rnorm(1,mean=0,sd=sigma[n]) # add in the AR(1) noise
    }
}
### superimpose residuals on the best hindcast from the SLR model to make N bootstrap samples###
for(i in 1:N) {
    slr.boot_hind[i,]=slr.est$sle+res.boot_hind[i,]
}


###IMPORTANT: calculate polynomial coefficients for the bootstrapped samples###
#------------- Find Plausible Parameters & Calculate Hindcast Sea-level Rise with Uncertainty ------------------
### calculate polynomial coefficients to obtain "N" number of simulated observation vectors###
from=2
to=hindcast_length #122
parm=mat.or.vec(N, 2) # parm is a vector of the parameter a and T0
initialvalue=rep(NA, N) # vector of the parameter of sea-level rise in 1880
for(i in 1:N) {
  timestep=1
  source("Scripts/Deoptim_rahm_model.R")
  source("Scripts/min_res_bootstrapped.R")
  lower=c(0,-3,err_neg[1])  
  upper=c(2,2,err_pos[1])
  iter=100  # specify number of iterations
  outDE <- DEoptim(min_res, lower, upper, 
                        DEoptim.control(itermax=iter,trace=FALSE))
  parm[i,1]=outDE$optim$bestmem[1]
  parm[i,2]=outDE$optim$bestmem[2]
  initialvalue[i]=outDE$optim$bestmem[3]
}
## Check for convergence: 
# save.image(file = "Workspace/bootstrap_seed1780.RData") # seed 3
# save.image(file = "Workspace/bootstrap_seed1234.RData") # seed 5
# save.image(file = "Workspace/bootstrap_seed1014.RData") # seed 4
# save.image(file = "Workspace/bootstrap_seed111.RData")  # seed 1
# save.image(file = "Workspace/bootstrap_seed1.RData")    # seed 2

# Calculate the smooth hindcast fits
boot.rate=mat.or.vec(N, nyears.obs)
boot.fit=mat.or.vec(N, nyears.obs) # boot.fit is the hindcast SLR simulation without noise
source("Scripts/searate_func.R")  #searate function finds the hindcast rates for SLR
for(n in 1:N) {
  boot.rate[n,] = thermal(parm[n,], hist.temp)[[1]]
  boot.fit[n,1] = initialvalue[n]
  for (i in from:to){
    boot.fit[n,i]=boot.fit[n,i-1]+boot.rate[n,i-1]*timestep
}
}
# Add the AR(1) residual noise to the smooth hindcast fits
prob.boot_hind=mat.or.vec(N, nyears.obs) #(nr,nc)
for(i in 1:N) {
    prob.boot_hind[i,]=boot.fit[i,]+res.boot_hind[i,]
}

#----------------------------- Project Sea-level Rise with Uncertainty --------------------------------
#Calculate smooth projections from parameters and RCP8.5 temp. emission
years.mod=(alltime) # all time represent the years from 1880 to 2300
nyears.mod=length(years.mod)
b.pred.rate=mat.or.vec(N, nyears.mod) #(nr,nc)
fit.boot_proj=mat.or.vec(N, nyears.mod) #(nr,nc) # fit.boot_proj is the projected SLR simulation without noise
to=projection_length #421
for(n in 1:N) {
  b.pred.rate[n,] = thermal(parm[n,], rcp85)[[1]] #find projected rates for SLR to 2300
  fit.boot_proj[n,1] = initialvalue[n]
  for (i in from:to){
    fit.boot_proj[n,i]=fit.boot_proj[n,i-1]+b.pred.rate[n,i-1]*timestep
}
}

### Calculate projected residuals with an lag-1 autocorrelation coefficient: rho[2]
res.boot_proj=mat.or.vec(N, nyears.mod) #(nr,nc)
slr.boot_proj=res.boot_proj
for(n in 1:N) {
  for(i in 2:nyears.mod) {
    res.boot_proj[n,i] = rho[2]*res.boot_proj[n,i-1] + rnorm(1,mean=0,sd=sigma[n])
  }
}

### superimpose residuals on polynomial smooth fits to
### calculate Probabilistic bootstrap samples###
for(i in 1:N) {
  slr.boot_proj[i,]=fit.boot_proj[i,]+res.boot_proj[i,]
}

#-------------------------- Estimate Parameter PDFs & Best Estimates ----------------------------
# Find the probability density function for each of the estimated parameters
b.pdfalpha <- density(parm[,1])
b.pdfbasetemp <- density(parm[,2])
b.pdfinitialval <- density(initialvalue)

# Find best estimated parameters from bootstrapping with the median value
b.new = c(median(parm[,1]), median(parm[,2]), median(initialvalue))
print(b.new)

# Calculate the best fit hindcast and projection
source("Scripts/sealevel_rahm_model.R")
to=hindcast_length #122
b.new.slr = rahmfunction(b.new, hist.temp) #best fit hindcast

to=projection_length #421 #best fit projection
b.new.rcp85 = rahmfunction(b.new, rcp85) #~4C increase from 1990-2100

#--------------------- Estimate PDF, CDF, and SF of SLR in 2100 & 2050 --------------------------
# Source survival function Function"
source("Scripts/plot_sf.r")

# Find the probability density function of sea-level estimates in 2100
b.prob_proj2100=mat.or.vec(N,1)
b.prob_proj2100=slr.boot_proj[,221] #The year 2100 is the 223 number in the sequence
b.pdf2100 <- density(b.prob_proj2100/100)

b.cdf2100 = ecdf(b.prob_proj2100/100)# Find the cumulative density function of SLR 2100
b.survival2050 <- plot.sf(slr.boot_proj[,171]/100, make.plot=F)  #Finds the survival function

# Find the probability density function of sea-level estimates in 2050
b.prob_proj2050=mat.or.vec(N,1)
b.prob_proj2050=slr.boot_proj[,171] #The year 2050 is the 171 number in the sequence
b.pdf2050 <- density(b.prob_proj2050/100)

b.cdf2050 = ecdf(b.prob_proj2050/100) # Find the cumulative density function of SLR 2050
b.survival2100 <- plot.sf(slr.boot_proj[,221]/100, make.plot=F)  #Finds the survival function

##################################### END ####################################################
#save.image(file = "mega_R_methods_workspace_bootstrap.RData")
