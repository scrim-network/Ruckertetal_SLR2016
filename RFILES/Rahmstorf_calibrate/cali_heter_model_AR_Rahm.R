#################################################################################
#
#  -file = "cali_heter_model_AR_Rahm.R"   Code written August 2014
#  - Author: Kelsey Ruckert (klr324@psu.edu)
#
#  -This program runs a Markov Chain Monte Carlo analysis of global sea-level
#       assuming heteroskedastic errors as described in Ruckert et al. (accepted).
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
library(mcmc)

# Set the seed
set.seed(111)

## Run multiple times with different seeds to check for convergence &
## robustness of the results:
# set.seed(1780)
# set.seed(1)
# set.seed(1014)
# set.seed(1234)

# Read in the sea-level rise observations, observation errors, years, and historic and emission temps.
source("../Data/temp_sea_2300.R")
hindcast_length=122 # there are 122 years from 1880 to 2002
projection_length=421 # from 1880 to 2300

#------------------------------ Find Initial Parameter & Initial Hindcast ------------------------
# Physical model parameters
# [1] alpha =  .34      sensitivity of SLR to temperature change (cm/year/C)
# [2] T_0   = -0.5      baseline temp at which the sea level anomaly is zero (C)
# [3] H_0   = -15      initial sea-level anomaly (cm)

timestep = 1 # timesteps are 1 year
from = 2 # start from the second year since the first year is an uncertain parameter, H_0
to=hindcast_length #122

# Run differential evolution optimization to find initial starting values for
# parameters to use in optimization of the likelihood function.
source("../Scripts/Deoptim_rahm_model.R") # physical model
source("../Scripts/minimize_residuals.R") # function to minimize the residuals

lower=c(0, -3, err_neg[1])
upper=c(2,  2, err_pos[1])
iter=1000  # specify number of iterations
outDEoptim <- DEoptim(min_res, lower, upper, 
                      DEoptim.control(itermax=iter,
                                      trace=FALSE))
print(outDEoptim$optim$bestmem) # print best initial parameters
deoptim.parameters = c(outDEoptim$optim$bestmem[1], outDEoptim$optim$bestmem[2], outDEoptim$optim$bestmem[3])

#------------------------ Calculate the Residuals  & AR(1) Coefficient --------------------------
# Load the physical sea-level model converted to R from the equaitons in Rahmstorf (2007).
source("../Scripts/sealevel_rahm_model.R")

# Use the optimized parameters to generate a fit to the data.
slr.est = rahmfunction(deoptim.parameters, hist.temp)

# Plot check that model simulation fits the data.
#plot(year, slr/100, pch=20, ylab="Sea-level anomaly [m]", xlab="Year")
#lines(year, slr.est$sle/100, col="blue", lwd=2)

# Calculate residuals from the fit to the data. Equation (S5)
res = slr-slr.est$sle

# Apply the auto-correlation function to determine a starting value for rho,
# the correlation coefficient.
#pdf(file="acf.pdf", family="Helvetica", pointsize=11, height=4.5, width=4.5)
rho=rep(NA,3)
ac=acf(res, lag.max=5, plot = FALSE, main="")
rho[1]=ac$acf[1]
rho[2]=ac$acf[2]
rho[3]=ac$acf[3]
rho[4]=ac$acf[4]
rho[5]=ac$acf[5]
#dev.off()

#--------------------------------- Run MCMC Calibration ----------------------------------------
# Set up priors.
bound.lower = c(0, -3, err_neg[1], 0, -0.99)
bound.upper = c(2,  2, err_pos[1], 1,  0.99)

# Name the model parameters and specify the number of model parameters.
# Sigma and rho are statistical parameters and are not counted in the number.
parnames = c("alpha","base temp","initialvalue", "sigma.y", "rho.y")
model.p = 3

# Set the measurement errors for heteroskedastic assumption.
y.meas.err = err.obs

# Load the likelihood model assuming correlated residuals.
source("../Scripts/Robs_likelihood_AR.R")

# Optimize the likelihood function to estimate initial starting values.
p = c(outDEoptim$optim$bestmem[1], outDEoptim$optim$bestmem[2], 
      outDEoptim$optim$bestmem[3], sd(res), rho[2]) 
p0 = c(0.34, -0.5, slr[1], 0.6, 0.5) # Rahmstorf estimated best guess parameters
p0 = optim(p0, function(p) -log.post(p))$par
print(round(p0,4))

# Set the step size and number of iterations.
step = c(0.02, 0.02, 0.1, 0.01, 0.01)
NI = 2.5e7

# Run MCMC calibration.
mcmc.out1 = metrop(log.post, p0, nbatch=NI, scale=step)
prechain1 = mcmc.out1$batch

# Print the acceptance rate as a percent. Should be ~ 25%
acceptrate = mcmc.out1$accept * 100
cat("Accept rate =", acceptrate, "%\n")

#-------------------------- Estimate Parameter PDFs & Median Estimates ----------------------------
# Identify the burn-in period and subtract it from the chains.
burnin = seq(1, 0.01*NI, 1) # 1% burnin
hetChainBurnin <- prechain1[-burnin,]

# Find the probability density function for each of the estimated parameters.
heter.pdfa <- density(hetChainBurnin[ ,1])
heter.pdfTo <- density(hetChainBurnin[ ,2])
heter.pdfinitialval <- density(hetChainBurnin[ ,3])
heter.pdfsigma <- density(hetChainBurnin[ ,4])
heter.pdfrho <- density(hetChainBurnin[ ,5])

# Find median of the estimated parameters.
heter.med = c(median(hetChainBurnin[ ,1]), median(hetChainBurnin[ ,2]), median(hetChainBurnin[ ,3]))
print(heter.med)

# Estimate model hindcast from parameter medians.
to=hindcast_length #122
heter.med.hindcast = rahmfunction(heter.med, hist.temp)

# Estimate model projection from parameter medians.
to = projection_length #421
heter.med.projection = rahmfunction(heter.med, rcp85)

# Thin the chain to a subset; ~20,000 is sufficient.
heter_sub_chain = prechain1[seq(length(burnin), NI, 1237), ]
heter_subset_length = length(heter_sub_chain[ ,1])

## Run multiple times with different seeds to check for convergence &
## robustness of the results:
# save.image(file = "Workspace/heter1780.RData") # seed 3
# save.image(file = "Workspace/heter1234.RData") # seed 5
# save.image(file = "Workspace/heter1014.RData") # seed 4
# save.image(file = "Workspace/heter111.RData")  # seed 1
# save.image(file = "Workspace/heter1.RData")    # seed 2

# Check for simularities between full chain and the subset.
par(mfrow=c(3,2))
for(i in 1:5){
    plot(density(hetChainBurnin[ ,i]), type="l",
    xlab=paste('Parameter =',' ', parnames[i], sep=''), ylab="PDF", main="")
    lines(density(heter_sub_chain[ ,i]), col="red")
}

#------------- Hindcast Sea-level Rise with Uncertainty & Find Plausible Parameters ------------------
# Extract parameter vectors from the chain to enhance code readability.
alpha.heter.chain = heter_sub_chain[ ,1]
T_0.heter.chain = heter_sub_chain[ ,2]
H_0.heter.chain = heter_sub_chain[ ,3]
sigma.heter.chain = heter_sub_chain[ ,4]
rho.heter.chain = heter_sub_chain[ ,5]

# Set up empty matrices for sea level rate and sea level output.
nyears.obs=length(year) #number of years in observational time series
new.RATE = mat.or.vec(heter_subset_length, nyears.obs)
mcmc.fit = mat.or.vec(heter_subset_length, nyears.obs)

to=hindcast_length

# Loop over the sea level model to generate a distribution of sea level rates
# and sea level simulations.
for(i in 1:heter_subset_length) {
    # Estimate the sea level rate of change: equation (S17)
    new.RATE[i, ] = alpha.heter.chain[i]*(hist.temp - T_0.heter.chain[i])
    mcmc.fit[i,1] = H_0.heter.chain[i]  # Initial value
    
    # Use Forward Euler to estimate sea level over time.
    for (n in from:to){
        mcmc.fit[i,n] = mcmc.fit[i,n-1] + new.RATE[i,n-1]*timestep
    }
}

# Estimate the residuals with the AR(1) coefficient and sigma.
Resid_hindcast_heter = mat.or.vec(heter_subset_length, nyears.obs) #(nr,nc)
for(n in 1:heter_subset_length) {
    for(i in 2:nyears.obs) {
        # Equation (4-5 & S3-S4)
        Resid_hindcast_heter[n,i] = rho.heter.chain[n]*Resid_hindcast_heter[n,i-1] +
          rnorm(1, mean = 0, sd = sigma.heter.chain[n])
    }
}

# Estimate the hindcasts: add the residuals onto the model simulations. Equation (2) & (S1)
SLR.heter.hindcasts = mat.or.vec(heter_subset_length, nyears.obs) #(nr,nc)
for(i in 1:heter_subset_length) {
    SLR.heter.hindcasts[i,] = mcmc.fit[i,] + Resid_hindcast_heter[i,]
}

#----------------------------- Project Sea-level Rise with Uncertainty --------------------------------
years.mod=(alltime) # all time represent the years from 1880 to 2300
nyears.mod=length(years.mod)
to=projection_length #421

# Set up empty matrices for sea level rate and sea level output.
proj.heter.RATE = mat.or.vec(heter_subset_length, nyears.mod) #(nr,nc)
proj.heter.sim = mat.or.vec(heter_subset_length, nyears.mod) #(nr,nc)

# Loop over the sea level model to generate a distribution of sea level rates
# and sea level simulations.
for(n in 1:heter_subset_length) {
    # Estimate the sea level rate of change: equation (S17)
    proj.heter.RATE[n, ] = alpha.heter.chain[n]*(rcp85 - T_0.heter.chain[n])
    proj.heter.sim[n,1] = H_0.heter.chain[n] # initial value in 1880
    
    # Use Forward Euler to estimate sea level over time.
    for (i in from:to){
        proj.heter.sim[n,i] = proj.heter.sim[n,i-1] + proj.heter.RATE[n,i-1]*timestep
    }
}

# Estimate the residuals with the AR(1) coefficient and sigma.
Resid_projection_heter = mat.or.vec(heter_subset_length, nyears.mod) #(nr,nc)
for(n in 1:heter_subset_length) {
    for(i in 2:nyears.mod) {
        # Equation (4-5 & S3-S4)
        Resid_projection_heter[n,i] = rho.heter.chain[n]*Resid_projection_heter[n,i-1] +
          rnorm(1, mean = 0,sd = sigma.heter.chain[n])
    }
}

# Estimate the projections: add the residuals onto the model simulations. Equation (2) & (S1)
SLR.projections.heter = mat.or.vec(heter_subset_length, nyears.mod) #(nr,nc)
for(i in 1:heter_subset_length) {
    SLR.projections.heter[i,] = proj.heter.sim[i,] + Resid_projection_heter[i,]
}

#--------------------- Estimate PDF, CDF, and SF of SLR in 2100 & 2050 --------------------------
# Load survival function Function
source("../Scripts/plot_sf.r")

# Set up a vector for the sea-level anomaly distribution in 2050.
# The year 2050 is the 171 number in the sequence.
prob_proj2050_heter = mat.or.vec(heter_subset_length, 1)
prob_proj2050_heter = SLR.projections.heter[ ,171]

# Estimate the probability density function of sea-level anomalies in 2050.
pdf2050_heter <- density(prob_proj2050_heter/100)

# Estimate the cumulative density function.
cdf2050_heter = ecdf(prob_proj2050_heter/100)

# Estimate the survivial function.
survival2050_heter <- plot.sf(prob_proj2050_heter/100, make.plot=F)

#---

# Set up a vector for the sea-level anomaly distribution in 2100.
# The year 2100 is the 221 number in the sequence.
prob_proj2100_heter = mat.or.vec(heter_subset_length, 1)
prob_proj2100_heter = SLR.projections.heter[ ,221]

# Estimate the probability density function of sea-level anomalies in 2100.
pdf2100_heter <- density(prob_proj2100_heter/100)

# Estimate the cumulative density function.
cdf2100_heter = ecdf(prob_proj2100_heter/100)

# Estimate the survivial function.
survival2100_heter <- plot.sf(prob_proj2100_heter/100, make.plot=F)

#### SAVE THE CURRENT WORKSPACE ####
#save.image(file = "mega_R_methods_workspace_hetcon.RData")
############################################ END #############################################
