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
library(ncdf)
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

#--------------------- Step 2: Recreate the Rahmstorf predictions: ----------------
# Source in the function specifying the model (found in Rahmstorf 2007):
source("Scripts/sealevel_rahm_model.R")

# Set up the parameters specified in Rahmstorf 2007:
# The parameters in the model include (in order of appreciance in the vector)
# a = sensitivity of sea-level to changes in temperature. Given as meters/Celsius/year
# T0 = base temperature or temperature when sea-level is zero. Given as Celius
# Initial = the value of sea-level in the year 1880. 1880 is the starting year for this
#       analysis. Given as meters.

original = c(0.34, -0.5, slr[1]) #slr[1] is equal to -0.1464 meters

# Run the hindcast in one year incriments from 1880 to 2002
from = 2 # Start from the 2nd year since the first value is defined as a uncertain parameter
to = 122 # Number of observations from 1880 to 2002
timestep = 1 # Data is in criments of one year
hindcast = rahmfunction(original, hist.temp)  # hist.temp are the historical temperatures

# Run the projections using the emission scenarios used in Rahmstorf 2007
to = 45     # There are 45 points in this projection from 1880 to 2100
timestep = 5 # Projections are in intervals of 5 years

# The names of the projections correspond to the names of the emission scenario
max_p = rahmfunction(original, max)
min_p = rahmfunction(original, min)
a1fi_p = rahmfunction(original, a1fi)
a1b_p = rahmfunction(original, a1b)
a1t_p = rahmfunction(original, a1t)
a2_p = rahmfunction(original, a2)
b1_p = rahmfunction(original, b1)
b2_p = rahmfunction(original, b2)

#-------------------- Step 3: Calculate Uncertainty ------------------------------
# Use three statistical methods with differening assumptions to calculate the uncertainty
#       of sea-level rise out to 2300 and parameter uncertainties:
## NOTE: Running all the methods takes ~2hr on a single high-performance core

# Method A: Bootstrap:
source("bootstrap_Rahm_new.R")

# Method B: Markov Chain Monte Carlo AR(1) Homoskedastic:
source("Rcali_homo_model_AR.R")

# Method C: Markov Chain Monte Carlo AR(1) Heteroskedastic:
source("Rcali_heter_model_AR.R")

#------------------- Step 4: Compare Estimates Parameters ------------------------
# Compare the max, median, min parameters estimated from each method:
# Original: Rahmstorf [2007]:
print(original)

# Method A: Bootstrap:
max.boot.all.par <- c(max(parm[,1]), max(parm[,2]), max(initialvalue), max(sigma))
min.boot.all.par <- c(min(parm[,1]), min(parm[,2]), min(initialvalue), min(sigma))

print(b.new) # Print the median model parameters (a, T0, H0)
print(median(sigma)) # Print the median sigma parameter (calculation NEGLECTS heteroskedasticity)
print(rho[2]) # Print the autocorrelation coefficient 
print(max.boot.all.par) # Print the maximum estimate of parameters (a, T0, H0, sigma)
print(min.boot.all.par) # Print the minimum estimate of parameters (a, T0, H0, sigma)

# Method B: Markov Chain Monte Carlo AR(1) Homoskedastic:
summary(homChainBurnin) # Print the median model parameters (a, T0, H0, sigma, rho)

# Method C: Markov Chain Monte Carlo AR(1) Heteroskedastic:
summary(hetChainBurnin) # Print the median model parameters (a, T0, H0, sigma, rho)

# Find mode and 99%
getmode <- function(v) { #Calculate mode http://www.tutorialspoint.com/r/r_mean_median_mode.htm
    uniqv <- unique(v)
    uniqv[which.max(tabulate(match(v, uniqv)))]
}
# Bootstrap
getmode(parm[,1]); getmode(parm[,2]); getmode(initialvalue); getmode(sigma); rho[2]

# homoskedastic
hom.rounded = round(homChainBurnin,5)
getmode(hom.rounded[,1]); getmode(hom.rounded[,2]); getmode(hom.rounded[,3])
getmode(hom.rounded[,4]); getmode(hom.rounded[,5])

# Heteroskedastic
het.rounded = round(hetChainBurnin,5)
getmode(het.rounded[,1]); getmode(het.rounded[,2]); getmode(het.rounded[,3])
getmode(het.rounded[,4]); getmode(het.rounded[,5])

# Print the 99% model parameters (a, T0, H0, sigma, rho)
for(i in 1:5){
  print(quantile(homChainBurnin[,i], 0.99))
}
  for(i in 1:5){
  print(quantile(hetChainBurnin[,i], 0.99))
  }
# "Bootstrap: (a, T0)
for(i in 1:2){
  print(quantile(parm[,i], 0.99))
}
# H0, sigma
quantile(initialvalue, 0.99); quantile(sigma, 0.99)

#------------------- Step 5: Save the Current Workspace --------------------------
save.image(file = "Workspace/mega_R_methods_workspace.RData")

#------------------- Step 6: Compare the results with plots ----------------------
# Compare all the different uncertainty methods against each other and
#       the Rahmstorf(2007) estimates:
## NOTE: All of the plots will be saved to the current working directory as .pdf

#source("PlotRuckert_etalNEW_v2.R")

###################################### END #######################################
