#################################################################################
#
#  -file = "Mega_Rahmstorf.R"   Code written August 2014
#  - Author: Kelsey Ruckert (klr324@psu.edu)
#
#  -This program loads and analyses global sea-level and temperature
#       data as described in Ruckert et al. (2016). It also produces
#       the graphs shown in that paper. The graphs will be saved as pdf files
#       in the current working directory. For further description and references,
#       please read the paper.
#
# THIS CODE IS PROVIDED AS-IS WITH NO WARRANTY (NEITHER EXPLICIT
# NOT IMPLICIT).  I SHARE THIS CODE IN HOPES THAT IT IS USEFUL,
# BUT I AM NOT LIABLE FOR THE BEHAVIOR OF THIS CODE IN YOUR OWN
# APPLICATION.  YOU ARE FREE TO SHARE THIS CODE SO LONG AS THE
# AUTHOR(S) AND VERSION HISTORY REMAIN INTACT.
#
#   -NOTE: This program recreates the estimates found in Rahmstorf (2007) and
#       applies several methods and assumptions to the Rahmstorf (2007)
#       semi-empirical sea-level model including:
#
#       -preserves the AR(1) structure of the observations
#       -creates random realizations of "natural variability"
#       -bootstraping correlated data
#       -Markov Chain Monte Carlo AR(1) Heteroskedastic errors
#       -Markov Chain Monte Carlo AR(1) Homoskedastic errors
#       -projects new realizations to 2300
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
#  -The output from this file once ran will be saved as:
#       -Workspace file: "mega_R_methods_workspace.RData"
#
###################################################################################

rm(list =ls()) #Clear global environment

# Load in R packages
#install.packages("mcmc")
library(mcmc)
library(ncdf4)
library(coda)
library(mvtnorm)
library(DEoptim)
library(compiler)
enableJIT(3)
enableJIT(3)

#--------------------- Step 1: Read in all data. ----------------------------------
# This includes historical temperatures, emission temperatures,
#       and sea-level data.
source("Data/temp_sea_2300.R")

#--------------------- Step 2: Reproduce the Rahmstorf projections: ----------------
# Load the physical sea-level model converted to R from the equaitons in Rahmstorf (2007).
source("Scripts/sealevel_rahm_model.R")

# Physical model parameters
# [1] alpha =  .34      sensitivity of SLR to temperature change (cm/year/C)
# [2] T_0   = -0.5      baseline temp at which the sea level anomaly is zero (C)
# [3] H_0   = -15      initial sea-level anomaly (cm)

# Set vector equal to parameters in Rahmstorf (2007).
original = c(0.34, -0.5, slr[1])

# Run the hindcast in one year incriments from 1880 to 2002.
from = 2 # Start
to = 122 # Number of observations from 1880 to 2002
timestep = 1
hindcast = rahmfunction(original, hist.temp)

# Run the projections using the temperature derived from the IPCC emission scenarios
to = 45     # There are 45 points in this projection from 1880 to 2100
timestep = 5

max_p = rahmfunction(original, max)
min_p = rahmfunction(original, min)
a1fi_p = rahmfunction(original, a1fi)
a1b_p = rahmfunction(original, a1b)
a1t_p = rahmfunction(original, a1t)
a2_p = rahmfunction(original, a2)
b1_p = rahmfunction(original, b1)
b2_p = rahmfunction(original, b2)

#-------------------- Step 3: Fit the model in three different ways to the observational data ------------------------------
## NOTE: Running all the methods takes ~5hr on a single high-performance core

# Method: Bootstrap (assumes AR(1) and homoskedastic residuals)
source("bootstrap_Rahm_new.R")

# Method: Bayesian homoskedastic (Markov Chain Monte Carlo assuming AR(1) and homoskedastic residuals)
source("Rcali_homo_model_AR.R")

# Method: Bayesian heteroskedastic (Markov Chain Monte Carlo assuming AR(1) and heteroskedastic residuals)
source("Rcali_heter_model_AR.R")

#------------------- Step 4: Compare the median, mode, and 99% parameter estimates ------------------------
# Calculate mode http://www.tutorialspoint.com/r/r_mean_median_mode.htm
getmode <- function(v) {
    uniqv <- unique(v)
    uniqv[which.max(tabulate(match(v, uniqv)))]
}

# Original: Rahmstorf [2007]:
print(original)

# Method: Bootstrap: --------------------
# Print the median estimates (a, T0, H0, sigma).
median(bootstrap_parameters[,1]); median(bootstrap_parameters[,2])
median(bootstrap_parameters[,3]); median(bootstrap_parameters[,4])

# Print the mode estimates (a, T0, H0, sigma, rho).
getmode(bootstrap_parameters[,1]); getmode(bootstrap_parameters[,2]);
getmode(bootstrap_parameters[,3]); getmode(bootstrap_parameters[,4]); rho[2]

# Print the 99% estimates (a, T0, H0, sigma).
for(i in 1:4){
    print(quantile(bootstrap_parameters[,i], 0.99))
}

# Method: Bayesian homoskedastic: --------------------
# Print the median estimates (a, T0, H0, sigma, rho).
median(homChainBurnin[,1]); median(homChainBurnin[,2]); median(homChainBurnin[,3])
median(homChainBurnin[,4]); median(homChainBurnin[,5])

# Print the mode estimates (a, T0, H0, sigma, rho).
hom.rounded = round(homChainBurnin,5)
getmode(hom.rounded[,1]); getmode(hom.rounded[,2]); getmode(hom.rounded[,3])
getmode(hom.rounded[,4]); getmode(hom.rounded[,5])

# Print the 99% estimates (a, T0, H0, sigma, rho).
for(i in 1:5){
    print(quantile(homChainBurnin[,i], 0.99))
}

# Method: Bayesian heteroskedastic: --------------------
# Print the median estimates (a, T0, H0, sigma, rho).
median(hetChainBurnin[,1]); median(hetChainBurnin[,2]); median(hetChainBurnin[,3])
median(hetChainBurnin[,4]); median(hetChainBurnin[,5])

# Print the mode estimates (a, T0, H0, sigma, rho).
het.rounded = round(hetChainBurnin,5)
getmode(het.rounded[,1]); getmode(het.rounded[,2]); getmode(het.rounded[,3])
getmode(het.rounded[,4]); getmode(het.rounded[,5])

# Print the 99% estimates (a, T0, H0, sigma, rho).
for(i in 1:5){
    print(quantile(hetChainBurnin[,i], 0.99))
}

#------------------- Step 5: Save the Current Workspace --------------------------
save.image(file = "Workspace/mega_R_methods_workspace.RData")

#------------------- Step 6: Compare the results in plots ----------------------
## NOTE: All of the plots will be saved to the current working directory as .eps
#source("PlotRuckert_etal.R")
###################################### END #######################################
