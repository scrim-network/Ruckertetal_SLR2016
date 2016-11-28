#################################################################################
#
#  -file = "cali_homo_model_AR_VR.R"   Code written Sept. 21 2016
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
#       heteroskedastic errors to the Vermeer and Rahmstorf (2009) semi-empirical sea-level model:
#
#       -preserving the AR(1) structure of the observations
#       -creating random realizations of "natural variability"
#       -projecting new realizations to 2300
#
#----------------------------- DATA DETAILS ------------------------------------
#  -Model can be found in Vermeer and Rahmstorf (2009)
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
hindcast_length=122 # there are 122 years from 1880 to 2002
projection_length=421 # from 1880 to 2300

#------------------------------ Find Initial Parameter & Initial Hindcast ------------------------
# Physical model parameters
# [1] alpha 
# [2] T0
# [3] H_0   
# [4] b

timestep=1 # timesteps are 1 year
from=2 # start from the second year since the first year is an uncertain parameter
to=hindcast_length #122

timestep=1 # timesteps are 1 year
from=2 # start from the second year since the first year is an uncertain parameter
to=hindcast_length #122

# Run differential evolution optimization to find initial starting values..
source("../Scripts/DEoptim_VR_model.R")      # physical model
source("../Scripts/minimize_VR_residuals.R") # function to minimize the residuals

lower=c(0,-3, err_neg[1],-1)
upper=c(2, 2, err_pos[1], 1)
iter=1000  # specify number of iterations
outDEoptim <- DEoptim(min_res, lower, upper, 
                      DEoptim.control(itermax=iter,
                                      trace=FALSE))
print(outDEoptim$optim$bestmem)# print best initial parameters
deoptim.parameters = c(outDEoptim$optim$bestmem[1], outDEoptim$optim$bestmem[2], outDEoptim$optim$bestmem[3], outDEoptim$optim$bestmem[4])

#------------------------ Calculate the Residuals  & AR(1) Coefficient --------------------------
# Load the physical sea-level model converted to R from the equations in Grinsted et al. (2010).
source("../Scripts/sealevel_VR_2009.R")

# Use the optimized parameters to generate a fit to the data.
slr.est = VR_SL_model(deoptim.parameters, hist.temp)

# Plot check that model simulation fits the data.
#plot(year, slr, pch=20, ylab="Sea-level anomaly [m]", xlab="Year")
#lines(year, slr.est$sle, col="blue", lwd=2)

# Calculate residuals from the fit to the data. Equation (S6)
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

#--------------------------------- Run MCMC Calibration ----------------------------------------
# Set up priors.
bound.lower = c(0, -3, err_neg[1], -1, 0, -0.99)
bound.upper = c(2,  2, err_pos[1],  1, 1,  0.99)

# Name the model parameters and specify the number of model parameters.
# Sigma and rho are statistical parameters and are not counted in the number.
model.p=4
parnames=c("alpha","base temp","initialvalue","b", "sigma.y", "rho")

# Set the measurement errors for homoskedastic assumption.
y.meas.err = rep(0, length(slr))

# Load the likelihood model assuming correlated residuals.
source("../Scripts/VR_obs_likelihood_AR.R")

# Optimize the likelihood function to estimate initial starting values.
p = c(outDEoptim$optim$bestmem[1], outDEoptim$optim$bestmem[2], 
      outDEoptim$optim$bestmem[3],outDEoptim$optim$bestmem[4]/2, sd(res), rho[2]) 
p0 = c(0.34, 0.2, slr[1], 0.005, 0.6, 0.5)
p0 = optim(p0, function(p) -log.post(p))$par
print(round(p0,4))

# Set the step size and number of iterations.
step = c(0.01, 0.01, 0.01, 0.0001, 0.0001, 0.01)
NI = 5e6

library(mcmc)
library(adaptMCMC)

# # Run MCMC calibration.
# mcmc.out.homo = metrop(log.post, p0, nbatch=NI, scale=step)
# homoskchain = mcmc.out.homo$batch
# 
# # Print the acceptance rate as a percent. Should be ~ 25%
# acceptrate = mcmc.out.homo$accept * 100
# cat("Accept rate =", acceptrate, "%\n")

# Set optimal acceptance rate as # parameters->infinity (Gelman et al, 1996; Roberts et al, 1997)
# Set number of iterations and rate of adaptation (between 0.5 and 1, lower is faster adaptation)
accept.mcmc = 0.234
gamma.mcmc = 0.5											
step.mcmc = step

# Run MCMC calibration.
mcmc.out.homo = MCMC(log.post, NI, p0, scale=step.mcmc, adapt=TRUE, acc.rate=accept.mcmc,gamma=gamma.mcmc, list=TRUE,
                 n.start=round(0.01*NI))
homoskchain = mcmc.out.homo$samples

# Identify the burn-in period and subtract it from the chains.
burnin = seq(1, 0.01*NI, 1) # 1% burnin
homChainBurnin <- homoskchain[-burnin,]

#--------------------------------- Test for convergence ----------------------------------------
library(coda)
# heidel.diag(prechain1, eps=0.1, pvalue=0.05)
#source("Scripts/Trace_plots.R")
set.seed(1780)

mcmc.out1780= MCMC(log.post, NI, p0, scale=step.mcmc, adapt=TRUE, acc.rate=accept.mcmc,gamma=gamma.mcmc, list=TRUE,
                   n.start=round(0.01*NI))
prechain1780 = mcmc.out1780$samples
# mcmc.out1780 = metrop(log.post, p0, nbatch=NI, scale=step)
# prechain1780 = mcmc.out1780$batch

homo = as.mcmc(homoskchain)
homo2 = as.mcmc(prechain1780)
homolist = mcmc.list(list(homo, homo2))
gelman.diag(homolist)

# Remove the second MCMC run
rm(mcmc.out1780, prechain1780, homo, homo2, homolist)

# set seed back to original
set.seed(111)

#-------------------------- Estimate Parameter PDFs & Median Estimates ----------------------------
# Find the probability density function for each of the estimated parameters
homo.pdfa <- density(homChainBurnin[ ,1])
homo.pdfT_0 <- density(homChainBurnin[ ,2])
homo.pdfinitialval <- density(homChainBurnin[ ,3])
homo.pdfb <- density(homChainBurnin[ ,4])
homo.pdfsigma <- density(homChainBurnin[ ,5])
homo.pdfrho <- density(homChainBurnin[ ,6])

# Find median of the estimated parameters.
homo.med = c(median(homChainBurnin[ ,1]), median(homChainBurnin[ ,2]), median(homChainBurnin[ ,3]), median(homChainBurnin[ ,4]))
print(homo.med)

# Estimate model hindcast from parameter medians.
to=hindcast_length #122
homo.med.hindcast = VR_SL_model(homo.med, hist.temp)

# Estimate model projection from parameter medians.
to = projection_length
homo.med.projection = VR_SL_model(homo.med, rcp85) #421 #best fit projection

# Thin the chain to a subset; ~20,000 is sufficient.
homo_subset_length = 20000
homo_sub_chain = homChainBurnin[sample(nrow(homChainBurnin), size=homo_subset_length, replace=FALSE), ]

## Run multiple times with different seeds to check for convergence &
## robustness of the results:
# save.image(file = "Workspace/homo1780.RData") # seed 3
# save.image(file = "Workspace/homo1234.RData") # seed 5
# save.image(file = "Workspace/homo1014.RData") # seed 4
# save.image(file = "Workspace/homo111.RData")  # seed 1
# save.image(file = "Workspace/homo1.RData")    # seed 2

# Check for simularities between full chain and the subset.
pdf(file="simplecheck_homo_subset.pdf")
par(mfrow=c(3,2))
for(i in 1:6){
  plot(density(homChainBurnin[ ,i]), type="l",
       xlab=paste('Parameter =',' ', parnames[i], sep=''), ylab="PDF", main="")
  lines(density(homo_sub_chain[ ,i]), col="red")
}
dev.off()

#------------- Plausible Parameters ------------------
# Extract parameter vectors from the chain to enhance code readability.
alpha.homo.chain = homo_sub_chain[ ,1]
T_0.homo.chain = homo_sub_chain[ ,2]
H_0.homo.chain = homo_sub_chain[ ,3]
b.homo.chain = homo_sub_chain[ ,4]
sigma.homo.chain = homo_sub_chain[ ,5]
rho.homo.chain = homo_sub_chain[ ,6]

#----------------------------- Hindcast and Project Sea-level Rise with Uncertainty --------------------------------
years.mod=(alltime) # all time represent the years from 1880 to 2300
nyears.mod=length(years.mod)
to=projection_length #421

# Set up empty matrices for sea level output.
proj.homo.sim = mat.or.vec(homo_subset_length, nyears.mod) #(nr,nc)

# Loop over the sea level model to generate a distribution of sea level simulations.
for(n in 1:homo_subset_length) {
    homo.est = c(alpha.homo.chain[n], T_0.homo.chain[n], H_0.homo.chain[n], b.homo.chain[n])
    proj.homo.sim[n, ] = VR_SL_model(homo.est, rcp85)$sle
}

# Estimate the residuals with the AR(1) coefficient and sigma.
Resid_projection_homo = mat.or.vec(homo_subset_length, nyears.mod) #(nr,nc)
for(n in 1:homo_subset_length) {
  for(i in 2:nyears.mod) {
    # Equation (S4)
    Resid_projection_homo[n,i] = rho.homo.chain[n]*Resid_projection_homo[n,i-1] +
      rnorm(1, mean = 0,sd = sigma.homo.chain[n])
  }
}

# Estimate the projections: add the residuals onto the model simulations. Equation (2) & (S1)
SLR.projections.homo = mat.or.vec(homo_subset_length, nyears.mod) #(nr,nc)
for(i in 1:homo_subset_length) {
  SLR.projections.homo[i,] = proj.homo.sim[i,] + Resid_projection_homo[i,]
}

#--------------------- Estimate PDF, CDF, and SF of SLR in 2100 & 2050 --------------------------
# Load survival function Function
source("../Scripts/plot_sf.r")

# Set up a vector for the sea-level anomaly distribution in 2050.
# The year 2050 is the 171 number in the sequence.
prob_proj2050_homo = mat.or.vec(homo_subset_length, 1)
prob_proj2050_homo = SLR.projections.homo[ ,171]

# Estimate the probability density function of sea-level anomalies in 2050.
pdf2050_homo <- density(prob_proj2050_homo)

# Estimate the cumulative density function.
cdf2050_homo = ecdf(prob_proj2050_homo)

# Estimate the survivial function.
survival2050_homo <- plot.sf(prob_proj2050_homo, make.plot=F)

#---

# Set up a vector for the sea-level anomaly distribution in 2100.
# The year 2100 is the 221 number in the sequence.
prob_proj2100_homo = mat.or.vec(homo_subset_length, 1)
prob_proj2100_homo = SLR.projections.homo[ ,221]

# Estimate the probability density function of sea-level anomalies in 2100.
pdf2100_homo <- density(prob_proj2100_homo)

# Estimate the cumulative density function.
cdf2100_homo = ecdf(prob_proj2100_homo)

# Estimate the survivial function.
survival2100_homo <- plot.sf(prob_proj2100_homo, make.plot=F)

#### SAVE THE CURRENT WORKSPACE ####
save.image(file = "Workspace/Homoskedastic_vermeer_SLR_model_workspace.RData")
############################################ END #############################################
