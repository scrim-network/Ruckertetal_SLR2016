#################################################################################
#
#  -file = "bootstrap_grinsted.R"   Code written March 2014, Updated Sept. 22 2016
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
#   -NOTE: This program applies the bootstrap method to the Grinsted et al. (2010)
#       semi-empirical sea-level model:
#
#       -preserving the AR(1) structure of the observations
#       -creating random realizations of "natural variability"
#       -projecting new realizations to 2300
#
#----------------------------- DATA DETAILS ------------------------------------
#  -Model can be found in Grinsted et al. (2010)
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
###################################################################################

rm(list =ls()) #Clear all previous variables
library(DEoptim)
library(compiler)
enableJIT(3)

# Set the seed
set.seed(111)

## Run multiple times with different seeds to check for convergence &
## robustness of the results:
# set.seed(1780)
# set.seed(1)
# set.seed(1014)
# set.seed(1234)

# Read in the sea-level rise observations, observation errors, years, and historic and emission temps.
source("../Data/temp_sea_2300_grinstead.R")
hindcast_length = 122 # there are 122 years from 1880 to 2002
projection_length = 421 # from 1880 to 2300

#------------------------------ Find Initial Parameter ------------------------
# Physical model parameters
# [1] alpha 
# [2] beta   
# [3] H_0   
# [4] 1/tau   

timestep=1 # timesteps are 1 year
from=2 # start from the second year since the first year is an uncertain parameter
to=hindcast_length #122

# Run differential evolution optimization to find initial starting values..
source("../Scripts/Deoptim_grinsted_model.R")      # physical model
source("../Scripts/minimize_grinsted_residuals.R") # function to minimize the residuals

lower=c(0,-1, err_neg[1], 0)
upper=c(4, 2.5, err_pos[1], 1)
iter=1000  # specify number of iterations
outDEoptim <- DEoptim(min_res, lower, upper, 
                      DEoptim.control(itermax=iter,
                                      trace=FALSE))
print(outDEoptim$optim$bestmem)# print best initial parameters
deoptim.parameters = c(outDEoptim$optim$bestmem[1], outDEoptim$optim$bestmem[2], outDEoptim$optim$bestmem[3], outDEoptim$optim$bestmem[4])

#------------------------ Calculate the Residuals  & AR(1) Coefficient --------------------------
# Load the physical sea-level model converted to R from the equations in Grinsted et al. (2010).
source("../Scripts/sealevel_grinsted_model.R")

# Use the optimized parameters to generate a fit to the data.
slr.est = grinsted_sealevel(deoptim.parameters, hist.temp)

# Plot check that model simulation fits the data.
#plot(year, slr, pch=20, ylab="Sea-level anomaly [m]", xlab="Year")
#lines(year, slr.est$sle, col="blue", lwd=2)

# Calculate residuals from the fit to the data. Equation (S5)
res = slr - slr.est$sle

# Apply the auto-correlation function to determine rho,
# the correlation coefficient, and retain the time-dependent structure of the residuals.
#pdf(file="acf.pdf", family="Helvetica", pointsize=11, height=4.5, width=4.5)
rho=rep(NA,3)
ac=acf(res, lag.max=5, plot=TRUE, main="")
rho[1]=ac$acf[1]
rho[2]=ac$acf[2]
rho[3]=ac$acf[3]
rho[4]=ac$acf[4]
rho[5]=ac$acf[5]
#dev.off()

#---------- Bootstrap the Residuals & Save AR(1) Coefficient and Sigma ---------
nyears.obs = length(year) #number of years in observational time series

# Set the number of bootstrap samples.
NI.boot = 20000

# Set up empty matrices/vectors for bootstrapped residuals and the stationary variance output.
resample.RESID = mat.or.vec(NI.boot, nyears.obs)
sigma.boot_sd = rep(NA, NI.boot)

# Draw replicates of the residuals by resampling the residuals with replacement to form bootstrap samples.
# This method assumes errors are homoskedastic
for(i in 1:NI.boot){
  resample.RESID[i,1:nyears.obs] = sample(res, size = nyears.obs, replace=TRUE)
  
  # Estimate the stationary process variance for each bootstrap sample
  sigma.boot_sd[i] = sd(resample.RESID[i,])
}

# IMPORTANT: estimate the uncorrelated white noise variance from the stationary process vaiance to avoid accounting
# for autocorrelation twice. Equation (S6)
boot.sigma = sqrt((sigma.boot_sd^2)*(1-(rho[2]^2)))

# Estimate the bootstrapped residuals with the AR(1) coefficient and sigma.
boostrapped_RES_hind = mat.or.vec(NI.boot, nyears.obs) #(nr,nc)
for(n in 1:NI.boot) {
  for(i in 2:nyears.obs) {
    # Equation (S4)
    boostrapped_RES_hind[n,i] = rho[2]*boostrapped_RES_hind[n,i-1] +
      rnorm(1, mean = 0, sd = boot.sigma[n])
  }
}

# Estimate the bootstrap simulations: : add the residuals onto the best fit model simulation. Equation (2) & (S1)
slr.boot.sim = mat.or.vec(NI.boot, nyears.obs) #(nr,nc)
for(i in 1:NI.boot) {
  slr.boot.sim[i,] = slr.est$sle + boostrapped_RES_hind[i,]
}

#------------- Estimate model parameters from each boostrap simulation ------------------
### calculate polynomial coefficients to obtain "N" number of simulated observation vectors###
from=2
to=hindcast_length #122
timestep=1

# Set up empty matrix for parameter outputs.
bootstrap_parameters = mat.or.vec(NI.boot, 5)
bootstrap_parameters[ ,5] = boot.sigma

# Loop over the DEoptim command to generate a distribution of model parameters.
for(i in 1:NI.boot) {
  source("../Scripts/min_res_grinsted_bootstrapped.R") # function to minimize the residuals of the bootstrap simulations
  
  lower=c(0,-1, err_neg[1], 0)
  upper=c(4, 2.5, err_pos[1], 1)
  iter=100  # specify number of iterations
  outDE <- DEoptim(min_res, lower, upper,
                   DEoptim.control(itermax=iter,trace=FALSE))
  
  bootstrap_parameters[i, 1]=outDE$optim$bestmem[1]
  bootstrap_parameters[i, 2]=outDE$optim$bestmem[2]
  bootstrap_parameters[i, 3]=outDE$optim$bestmem[3]
  bootstrap_parameters[i, 4]=outDE$optim$bestmem[4]
}

## Run multiple times with different seeds to check for robustness of the results:
# save.image(file = "Workspace/bootstrap_seed1780.RData") # seed 3
# save.image(file = "Workspace/bootstrap_seed1234.RData") # seed 5
# save.image(file = "Workspace/bootstrap_seed1014.RData") # seed 4
# save.image(file = "Workspace/bootstrap_seed111.RData")  # seed 1
# save.image(file = "Workspace/bootstrap_seed1.RData")    # seed 2

#------------- Hindcast Sea-level Rise with Uncertainty ------------------
# Extract parameter vectors from the distributions to enhance code readability.
alpha.boot = bootstrap_parameters[ ,1]
beta.boot = bootstrap_parameters[ ,2]
H_0.boot = bootstrap_parameters[ ,3]
tau.boot = bootstrap_parameters[ ,4]

# Set up empty matrices for sea level rate and sea level output.
boot.Seq = mat.or.vec(NI.boot, nyears.obs)
boot.RATE = mat.or.vec(NI.boot, nyears.obs)
boot.fit = mat.or.vec(NI.boot, nyears.obs)

to = hindcast_length

# Loop over the sea level model to generate a distribution of sea level rates
# and sea level simulations.
for(n in 1:NI.boot) {
  # Estimate the sea level rate of change: equation (1)
  boot.Seq[n, ] = alpha.boot[n]*hist.temp + beta.boot[n]
  boot.fit[n, 1] = H_0.boot[n]
  boot.RATE[n, 1] = (boot.Seq[n,1] - boot.fit[n,1])*tau.boot[n]

  # Use Forward Euler to estimate sea level over time.
  for (i in from:to){
    boot.RATE[n,i] = (boot.Seq[n, i-1] - boot.fit[n, i-1])*tau.boot[n]
    boot.fit[n,i] = boot.fit[n,i-1] + boot.RATE[n,i-1]*timestep
  }
}

# Estimate the hindcasts: add the residuals onto the model simulations. Equation (2) & (S1)
SLR.boot.hindcasts = mat.or.vec(NI.boot, nyears.obs) #(nr,nc)
for(i in 1:NI.boot) {
  SLR.boot.hindcasts[i,] = boot.fit[i,] + boostrapped_RES_hind[i,]
}

#----------------------------- Project Sea-level Rise with Uncertainty --------------------------------
years.mod=(alltime) # all time represent the years from 1880 to 2300
nyears.mod=length(years.mod)
to=projection_length #421

# Set up empty matrices for sea level rate and sea level output.
proj.boot.Seq = mat.or.vec(NI.boot, nyears.mod)
proj.boot.RATE = mat.or.vec(NI.boot, nyears.mod) #(nr,nc)
proj.boot.sim = mat.or.vec(NI.boot, nyears.mod) #(nr,nc)

# Loop over the sea level model to generate a distribution of sea level rates
# and sea level simulations.
for(n in 1:NI.boot) {
  # Estimate the sea level rate of change: equation (1)
  proj.boot.Seq[n, ] = alpha.boot[n]*rcp85 + beta.boot[n]
  proj.boot.sim[n, 1] = H_0.boot[n]
  proj.boot.RATE[n, 1] = (proj.boot.Seq[n,1] - proj.boot.sim[n,1])*tau.boot[n]
  
  # Use Forward Euler to estimate sea level over time.
  for (i in from:to){
    proj.boot.RATE[n,i] = (proj.boot.Seq[n,i-1] - proj.boot.sim[n,i-1])*tau.boot[n]
    proj.boot.sim[n,i] = proj.boot.sim[n,i-1] + proj.boot.RATE[n,i-1]*timestep
  }
}

# Estimate the residuals with the AR(1) coefficient and sigma.
Resid_projection_boot = mat.or.vec(NI.boot, nyears.mod) #(nr,nc)
for(n in 1:NI.boot) {
  for(i in 2:nyears.mod) {
    # Equation (S4)
    Resid_projection_boot[n,i] = rho[2]*Resid_projection_boot[n,i-1] +
      rnorm(1, mean = 0, sd = boot.sigma[n])
  }
}

# Estimate the hindcasts: add the residuals onto the model simulations. Equation (2) & (S1)
SLR.projections.boot = mat.or.vec(NI.boot, nyears.mod) #(nr,nc)
for(i in 1:NI.boot) {
  SLR.projections.boot[i,] = proj.boot.sim[i,] + Resid_projection_boot[i,]
}

#-------------------------- Estimate Parameter PDFs & Median Estimates ----------------------------
# Find the probability density function for each of the estimated parameters.
boot.pdfa <- density(bootstrap_parameters[ ,1])
boot.pdfb <- density(bootstrap_parameters[ ,2])
boot.pdfinitialval <- density(bootstrap_parameters[ ,3])
boot.pdfTau <- density(bootstrap_parameters[ ,4])
boot.pdfsigma <- density(bootstrap_parameters[ ,5])

# Find median of the estimated parameters
boot.med = c(median(bootstrap_parameters[,1]), median(bootstrap_parameters[,2]), median(bootstrap_parameters[,3]), median(bootstrap_parameters[,4]))
print(boot.med)

# Estimate model hindcast from parameter medians.
to=hindcast_length #122
boot.med.hindcast = grinsted_sealevel(boot.med, hist.temp) # median fit hindcast

# Estimate model projection from parameter medians.
to=projection_length #421
boot.med.projection = grinsted_sealevel(boot.med, rcp85)

#--------------------- Estimate PDF, CDF, and SF of SLR in 2100 & 2050 --------------------------
# Load survival function Function
source("../Scripts/plot_sf.r")

# Set up a vector for the sea-level anomaly distribution in 2050.
# The year 2050 is the 171 number in the sequence.
prob_proj2050_boot = mat.or.vec(NI.boot, 1)
prob_proj2050_boot = SLR.projections.boot[ ,171]

# Estimate the probability density function of sea-level anomalies in 2050.
pdf2050_boot <- density(prob_proj2050_boot)

# Estimate the cumulative density function.
cdf2050_boot = ecdf(prob_proj2050_boot)

# Estimate the survivial function.
survival2050_boot <- plot.sf(prob_proj2050_boot, make.plot=F)

#---

# Set up a vector for the sea-level anomaly distribution in 2100.
# The year 2100 is the 221 number in the sequence.
prob_proj2100_boot = mat.or.vec(NI.boot, 1)
prob_proj2100_boot = SLR.projections.boot[ ,221]

# Estimate the probability density function of sea-level anomalies in 2100.
pdf2100_boot <- density(prob_proj2100_boot)

# Estimate the cumulative density function.
cdf2100_boot = ecdf(prob_proj2100_boot)

# Estimate the survivial function.
survival2100_boot <- plot.sf(prob_proj2100_boot, make.plot=F)

#### SAVE THE CURRENT WORKSPACE ####
save.image(file = "Workspace/Bootstrap_Grinsted_SLR_model_workspace.RData")

##################################### END ####################################################

