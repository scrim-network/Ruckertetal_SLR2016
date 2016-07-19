##################################################################################################
#
#  -file = "RCP_temp_scenarios.R"   Code written April 2016
#  - Author: Kelsey Ruckert (klr324@psu.edu)
#
#  -This program loads in temperatures from various RCP scenarios. Using the different temperature
#       scenarios each method (Bootstrap, MCMC homoskedastic, and MCMC heteroskedastic) is projected.
#       This script produces supplementary figure 4 in Ruckert et al. (2016). For further description and
#       references, please read the paper.
#
#   -NOTE: The graph will be saved as a eps file in the working directory.
#
# THIS CODE IS PROVIDED AS-IS WITH NO WARRANTY (NEITHER EXPLICIT
# NOT IMPLICIT).  I SHARE THIS CODE IN HOPES THAT IT IS USEFUL,
# BUT I AM NOT LIABLE FOR THE BEHAVIOR OF THIS CODE IN YOUR OWN
# APPLICATION.  YOU ARE FREE TO SHARE THIS CODE SO LONG AS THE
# AUTHOR(S) AND VERSION HISTORY REMAIN INTACT.
#
#
###################################################################################################
# Various RCP temperature scenarios:
# install.packages('ncdf')
library(ncdf)
library(ncdf4)
source("Data/temp_sea_2300.R")

# RCP 2.6 -----------------------------------------------------------------
# fid1<-open.ncdf("Data/global.tas.aann.HadGEM2-ES.historical+rcp26.r1i1p1.18600101-22991230.nc")

# temp_K<-get.var.ncdf(fid1,"tas")
# year_nc<-get.var.ncdf(fid1,"time")

#Close file
# close.ncdf(fid1)

# If ncdf is unavailable un comment the following code and use ncdf4 library
fid1 <- nc_open("Data/global.tas.aann.HadGEM2-ES.historical+rcp26.r1i1p1.18600101-22991230.nc")
temp_K <- ncvar_get(fid1,"tas")
year_nc <- ncvar_get(fid1,"time")
nc_close(fid1)

# Convert Kelvin to Celsius
# [°C] = [K] − 273.15
temp_C = temp_K - 273.15

# Calculate anomaly with respect to 20th century mean
respect_20th = mean(temp_C[42:141]) #1901 - 2000
tempRCP26_anom = temp_C - respect_20th
tempRCP26_anom = c(data[1:134,2], tempRCP26_anom[155:241])

# RCP 4.5 -----------------------------------------------------------------
# fid1<-open.ncdf("Data/global.tas.aann.HadGEM2-ES.historical+rcp45.r1i1p1.18600101-22991230.nc")

# temp_K<-get.var.ncdf(fid1,"tas")
# year_nc<-get.var.ncdf(fid1,"time")

#Close file
# close.ncdf(fid1)

# If ncdf is unavailable un comment the following code and use ncdf4 library
fid1 <- nc_open("Data/global.tas.aann.HadGEM2-ES.historical+rcp45.r1i1p1.18600101-22991230.nc")
temp_K <- ncvar_get(fid1,"tas")
year_nc <- ncvar_get(fid1,"time")
nc_close(fid1)

temp_C = temp_K - 273.15
# Calculate anomaly with respect to 20th century mean
respect_20th = mean(temp_C[42:141]) #1901 - 2000
tempRCP45_anom = temp_C - respect_20th
tempRCP45_anom = c(data[1:134,2], tempRCP45_anom[155:241])

# RCP 6.0 -----------------------------------------------------------------
# fid1<-open.ncdf("Data/global.tas.aann.HadGEM2-ES.historical+rcp60.r1i1p1.18600101-20991230.nc")

# temp_K<-get.var.ncdf(fid1,"tas")
# year_nc<-get.var.ncdf(fid1,"time")

#Close file
# close.ncdf(fid1)

# If ncdf is unavailable un comment the following code and use ncdf4 library
fid1 <- nc_open("Data/global.tas.aann.HadGEM2-ES.historical+rcp60.r1i1p1.18600101-20991230.nc")
temp_K <- ncvar_get(fid1,"tas")
year_nc <- ncvar_get(fid1,"time")
nc_close(fid1)

temp_C = temp_K - 273.15
# Calculate anomaly with respect to 20th century mean
respect_20th = mean(temp_C[42:141])
tempRCP60_anom = temp_C - respect_20th
tempRCP60_anom = c(data[1:134,2], tempRCP60_anom[155:241])

# RCP 8.5 -----------------------------------------------------------------
# fid1<-open.ncdf("Data/global.tas.aann.HadGEM2-ES.historical+rcp85.r1i1p1.18600101-22991230.nc")

# temp_K<-get.var.ncdf(fid1,"tas")
# year_nc<-get.var.ncdf(fid1,"time")

#Close file
# close.ncdf(fid1)

# If ncdf is unavailable un comment the following code and use ncdf4 library
fid1 <- nc_open("Data/global.tas.aann.HadGEM2-ES.historical+rcp85.r1i1p1.18600101-22991230.nc")
temp_K <- ncvar_get(fid1,"tas")
year_nc <- ncvar_get(fid1,"time")
nc_close(fid1)

temp_C = temp_K - 273.15
# Calculate anomaly with respect to 20th century mean
respect_20th = mean(temp_C[42:141])
tempRCP85_anom = temp_C - respect_20th
tempRCP85_anom = c(data[1:134,2], tempRCP85_anom[155:241])

remove(year_nc, temp_C, temp_K, respect_20th, fid1)

# Load saved workspace -----------------------------------------------------------------
load("Workspace/mega_R_methods_workspace.RData") # Load in the saved workspace from file Mega_Rahmstorf.R

#----------------------------- Project Sea-level Rise with Uncertainty --------------------------------
timestep=1
years.rcp=1880:2100 # all time represent the years from 1880 to 2100
to=length(years.rcp)

# Bootstrap RCP26 ---------------------------------------------------------
boot.rcp26.rate=mat.or.vec(NI.boot, length(years.rcp)) #(nr,nc)
boot.rcp26.proj=mat.or.vec(NI.boot, length(years.rcp)) #(nr,nc) # fit.boot_proj is the projected SLR simulation without noise

for(n in 1:NI.boot) {
  # Estimate the sea level rate of change: equation (1)
  boot.rcp26.rate[n,] = bootstrap_parameters[n,1]*(tempRCP26_anom - bootstrap_parameters[n,2])
  boot.rcp26.proj[n,1] = bootstrap_parameters[n,3]
  
  # Use Forward Euler to estimate sea level over time.
  for (i in from:to){
    boot.rcp26.proj[n,i]=boot.rcp26.proj[n,i-1]+boot.rcp26.rate[n,i-1]*timestep
  }
}

# Estimate the residuals with the AR(1) coefficient and sigma.
res.boot.26=mat.or.vec(NI.boot, length(years.rcp)) #(nr,nc)
boot.slr.rcp26=res.boot.26
for(n in 1:NI.boot) {
  for(i in 2:length(years.rcp)) {
    res.boot.26[n,i] = rho[2]*res.boot.26[n,i-1] + rnorm(1, mean = 0, sd = bootstrap_parameters[n,4])
  }
}

# Estimate the projections: add the residuals onto the model simulations. Equation (2) & (S1)
for(i in 1:NI.boot) {
  boot.slr.rcp26[i,]=boot.rcp26.proj[i,] + res.boot.26[i,]
}

remove(res.boot.26, boot.rcp26.proj, boot.rcp26.rate)
# Bootstrap RCP45 ---------------------------------------------------------
boot.rcp45.rate=mat.or.vec(NI.boot, length(years.rcp)) #(nr,nc)
boot.rcp45.proj=mat.or.vec(NI.boot, length(years.rcp)) #(nr,nc) # fit.boot_proj is the projected SLR simulation without noise

for(n in 1:NI.boot) {
  # Estimate the sea level rate of change: equation (1)
  boot.rcp45.rate[n,] = bootstrap_parameters[n,1]*(tempRCP45_anom - bootstrap_parameters[n,2])
  boot.rcp45.proj[n,1] = bootstrap_parameters[n,3]
  
  # Use Forward Euler to estimate sea level over time.
  for (i in from:to){
    boot.rcp45.proj[n,i]=boot.rcp45.proj[n,i-1]+boot.rcp45.rate[n,i-1]*timestep
  }
}

# Estimate the residuals with the AR(1) coefficient and sigma.
res.boot.45=mat.or.vec(NI.boot, length(years.rcp)) #(nr,nc)
boot.slr.rcp45=res.boot.45
for(n in 1:NI.boot) {
  for(i in 2:length(years.rcp)) {
    res.boot.45[n,i] = rho[2]*res.boot.45[n,i-1] + rnorm(1, mean = 0, sd = bootstrap_parameters[n,4])
  }
}

# Estimate the projections: add the residuals onto the model simulations. Equation (2) & (S1)
for(i in 1:NI.boot) {
  boot.slr.rcp45[i,] = boot.rcp45.proj[i,] + res.boot.45[i,]
}

remove(res.boot.45, boot.rcp45.proj, boot.rcp45.rate)
# Bootstrap RCP60 ---------------------------------------------------------
boot.rcp60.rate=mat.or.vec(NI.boot, length(years.rcp)) #(nr,nc)
boot.rcp60.proj=mat.or.vec(NI.boot, length(years.rcp)) #(nr,nc) # fit.boot_proj is the projected SLR simulation without noise

for(n in 1:NI.boot) {
  # Estimate the sea level rate of change: equation (1)
  boot.rcp60.rate[n,] = bootstrap_parameters[n,1]*(tempRCP60_anom - bootstrap_parameters[n,2])
  boot.rcp60.proj[n,1] = bootstrap_parameters[n,3]
  
  # Use Forward Euler to estimate sea level over time.
  for (i in from:to){
    boot.rcp60.proj[n,i]=boot.rcp60.proj[n,i-1]+boot.rcp60.rate[n,i-1]*timestep
  }
}

# Estimate the residuals with the AR(1) coefficient and sigma.
res.boot.60=mat.or.vec(NI.boot, length(years.rcp)) #(nr,nc)
boot.slr.rcp60=res.boot.60
for(n in 1:NI.boot) {
  for(i in 2:length(years.rcp)) {
    res.boot.60[n,i] = rho[2]*res.boot.60[n,i-1] + rnorm(1, mean = 0, sd = bootstrap_parameters[n,4])
  }
}

# Estimate the projections: add the residuals onto the model simulations. Equation (2) & (S1)
for(i in 1:NI.boot) {
  boot.slr.rcp60[i,]=boot.rcp60.proj[i,] + res.boot.60[i,]
}

remove(res.boot.60, boot.rcp60.proj, boot.rcp60.rate)
# Bootstrap RCP85 ---------------------------------------------------------
boot.rcp85.rate=mat.or.vec(NI.boot, length(years.rcp)) #(nr,nc)
boot.rcp85.proj=mat.or.vec(NI.boot, length(years.rcp)) #(nr,nc) # fit.boot_proj is the projected SLR simulation without noise

for(n in 1:NI.boot) {
  # Estimate the sea level rate of change: equation (1)
  boot.rcp85.rate[n,] =  bootstrap_parameters[n,1]*(tempRCP85_anom - bootstrap_parameters[n,2])
  boot.rcp85.proj[n,1] = bootstrap_parameters[n,3]
  
  # Use Forward Euler to estimate sea level over time.
  for (i in from:to){
    boot.rcp85.proj[n,i]=boot.rcp85.proj[n,i-1] + boot.rcp85.rate[n,i-1]*timestep
  }
}

# Estimate the residuals with the AR(1) coefficient and sigma.
res.boot.85=mat.or.vec(NI.boot, length(years.rcp)) #(nr,nc)
boot.slr.rcp85=res.boot.85
for(n in 1:NI.boot) {
  for(i in 2:length(years.rcp)) {
    res.boot.85[n,i] = rho[2]*res.boot.85[n,i-1] + rnorm(1, mean = 0, sd = bootstrap_parameters[n,4])
  }
}

# Estimate the projections: add the residuals onto the model simulations. Equation (2) & (S1)
for(i in 1:NI.boot) {
  boot.slr.rcp85[i,]=boot.rcp85.proj[i,] + res.boot.85[i,]
}

remove(res.boot.85, boot.rcp85.proj, boot.rcp85.rate)
###############################################################################
# Homoskedastic RCP26 ---------------------------------------------------------
homo.rcp26.rate=mat.or.vec(homo_subset_length, length(years.rcp)) #(nr,nc)
homo.rcp26.proj=mat.or.vec(homo_subset_length, length(years.rcp)) #(nr,nc)

for(n in 1:homo_subset_length) {
  # Estimate the sea level rate of change: equation (1)
  homo.rcp26.rate[n,] = homo_sub_chain[n,1]*(tempRCP26_anom - homo_sub_chain[n,2])
  homo.rcp26.proj[n,1] = homo_sub_chain[n,3]
  
  # Use Forward Euler to estimate sea level over time.
  for (i in from:to){
    homo.rcp26.proj[n,i]=homo.rcp26.proj[n,i-1] + homo.rcp26.rate[n,i-1]*timestep
  }
}

# Estimate the residuals with the AR(1) coefficient and sigma.
res.homo_26=mat.or.vec(homo_subset_length, length(years.rcp)) #(nr,nc)
homo.slr.rcp26=res.homo_26
for(n in 1:homo_subset_length) {
  for(i in 2:length(years.rcp)) {
    res.homo_26[n,i] = homo_sub_chain[n,5]*res.homo_26[n,i-1] +
      rnorm(1, mean = 0, sd = homo_sub_chain[n,4])
  }
}

# Estimate the projections: add the residuals onto the model simulations. Equation (2) & (S1)
for(i in 1:homo_subset_length) {
  homo.slr.rcp26[i,]=homo.rcp26.proj[i,] + res.homo_26[i,]
}

remove(res.homo_26, homo.rcp26.proj, homo.rcp26.rate)
# Homoskedastic RCP45 ---------------------------------------------------------
homo.rcp45.rate=mat.or.vec(homo_subset_length, length(years.rcp)) #(nr,nc)
homo.rcp45.proj=mat.or.vec(homo_subset_length, length(years.rcp)) #(nr,nc)

for(n in 1:homo_subset_length) {
  # Estimate the sea level rate of change: equation (1)
  homo.rcp45.rate[n,] = homo_sub_chain[n,1]*(tempRCP45_anom - homo_sub_chain[n,2])
  homo.rcp45.proj[n,1] = homo_sub_chain[n,3]
  
  # Use Forward Euler to estimate sea level over time.
  for (i in from:to){
    homo.rcp45.proj[n,i]=homo.rcp45.proj[n,i-1]+homo.rcp45.rate[n,i-1]*timestep
  }
}

# Estimate the residuals with the AR(1) coefficient and sigma.
res.homo_45=mat.or.vec(homo_subset_length, length(years.rcp)) #(nr,nc)
homo.slr.rcp45=res.homo_45
for(n in 1:homo_subset_length) {
  for(i in 2:length(years.rcp)) {
    res.homo_45[n,i] = homo_sub_chain[n,5]*res.homo_45[n,i-1] +
      rnorm(1, mean = 0, sd = homo_sub_chain[n,4])
  }
}

# Estimate the projections: add the residuals onto the model simulations. Equation (2) & (S1)
for(i in 1:homo_subset_length) {
  homo.slr.rcp45[i,]=homo.rcp45.proj[i,] + res.homo_45[i,]
}

remove(res.homo_45, homo.rcp45.proj, homo.rcp45.rate)
# Homoskedastic RCP60 ---------------------------------------------------------
homo.rcp60.rate=mat.or.vec(homo_subset_length, length(years.rcp)) #(nr,nc)
homo.rcp60.proj=mat.or.vec(homo_subset_length, length(years.rcp)) #(nr,nc)

for(n in 1:homo_subset_length) {
  # Estimate the sea level rate of change: equation (1)
  homo.rcp60.rate[n,] = homo_sub_chain[n,1]*(tempRCP60_anom - homo_sub_chain[n,2])
  homo.rcp60.proj[n,1] = homo_sub_chain[n,3]
  
  # Use Forward Euler to estimate sea level over time.
  for (i in from:to){
    homo.rcp60.proj[n,i]=homo.rcp60.proj[n,i-1] + homo.rcp60.rate[n,i-1]*timestep
  }
}

# Estimate the residuals with the AR(1) coefficient and sigma.
res.homo_60=mat.or.vec(homo_subset_length, length(years.rcp)) #(nr,nc)
homo.slr.rcp60=res.homo_60
for(n in 1:homo_subset_length) {
  for(i in 2:length(years.rcp)) {
    res.homo_60[n,i] = homo_sub_chain[n,5]*res.homo_60[n,i-1] +
      rnorm(1, mean = 0, sd = homo_sub_chain[n,4])
  }
}

# Estimate the projections: add the residuals onto the model simulations. Equation (2) & (S1)
for(i in 1:homo_subset_length) {
  homo.slr.rcp60[i,]=homo.rcp60.proj[i,] + res.homo_60[i,]
}

remove(res.homo_60, homo.rcp60.proj, homo.rcp60.rate)
# Homoskedastic RCP85 ---------------------------------------------------------
homo.rcp85.rate=mat.or.vec(homo_subset_length, length(years.rcp)) #(nr,nc)
homo.rcp85.proj=mat.or.vec(homo_subset_length, length(years.rcp)) #(nr,nc)

for(n in 1:homo_subset_length) {
  # Estimate the sea level rate of change: equation (1)
  homo.rcp85.rate[n,] = homo_sub_chain[n,1]*(tempRCP85_anom - homo_sub_chain[n,2])
  homo.rcp85.proj[n,1] = homo_sub_chain[n,3]
  
  # Use Forward Euler to estimate sea level over time.
  for (i in from:to){
    homo.rcp85.proj[n,i]=homo.rcp85.proj[n,i-1] + homo.rcp85.rate[n,i-1]*timestep
  }
}

# Estimate the residuals with the AR(1) coefficient and sigma.
res.homo_85=mat.or.vec(k, length(years.rcp)) #(nr,nc)
homo.slr.rcp85=res.homo_85
for(n in 1:homo_subset_length) {
  for(i in 2:length(years.rcp)) {
    res.homo_85[n,i] = homo_sub_chain[n,5]*res.homo_85[n,i-1] +
      rnorm(1, mean = 0, sd = homo_sub_chain[n,4])
  }
}

# Estimate the projections: add the residuals onto the model simulations. Equation (2) & (S1)
for(i in 1:homo_subset_length) {
  homo.slr.rcp85[i,] = homo.rcp85.proj[i,] + res.homo_85[i,]
}

remove(res.homo_85, homo.rcp85.proj, homo.rcp85.rate)
###############################################################################
# Heteroskedastic RCP26 ---------------------------------------------------------
heter.rcp26.rate=mat.or.vec(heter_subset_length, length(years.rcp)) #(nr,nc)
heter.rcp26.proj=mat.or.vec(heter_subset_length, length(years.rcp)) #(nr,nc)

for(n in 1:heter_subset_length) {
  # Estimate the sea level rate of change: equation (1)
  heter.rcp26.rate[n,] = heter_sub_chain[n,1]*(tempRCP26_anom - heter_sub_chain[n,2])
  heter.rcp26.proj[n,1] = heter_sub_chain[n,3]
  
  # Use Forward Euler to estimate sea level over time.
  for (i in from:to){
    heter.rcp26.proj[n,i]=heter.rcp26.proj[n,i-1] + heter.rcp26.rate[n,i-1]*timestep
  }
}

# Estimate the residuals with the AR(1) coefficient and sigma.
res.heter_26=mat.or.vec(heter_subset_length, length(years.rcp)) #(nr,nc)
heter.slr.rcp26=res.heter_26
for(n in 1:heter_subset_length) {
  for(i in 2:length(years.rcp)) {
    res.heter_26[n,i] = heter_sub_chain[n,5]*res.heter_26[n,i-1] +
      rnorm(1, mean = 0, sd = heter_sub_chain[n,4])
  }
}

# Estimate the projections: add the residuals onto the model simulations. Equation (2) & (S1)
for(i in 1:heter_subset_length) {
  heter.slr.rcp26[i,] = heter.rcp26.proj[i,] + res.heter_26[i,]
}

remove(res.heter_26, heter.rcp26.proj, heter.rcp26.rate)
# Heteroskedastic RCP45 ---------------------------------------------------------
heter.rcp45.rate=mat.or.vec(heter_subset_length, length(years.rcp)) #(nr,nc)
heter.rcp45.proj=mat.or.vec(heter_subset_length, length(years.rcp)) #(nr,nc)

for(n in 1:heter_subset_length) {
  # Estimate the sea level rate of change: equation (1)
  heter.rcp45.rate[n,] = heter_sub_chain[n,1]*(tempRCP45_anom - heter_sub_chain[n,2])
  heter.rcp45.proj[n,1] = heter_sub_chain[n,3]
  
  # Use Forward Euler to estimate sea level over time.
  for (i in from:to){
    heter.rcp45.proj[n,i]=heter.rcp45.proj[n,i-1] + heter.rcp45.rate[n,i-1]*timestep
  }
}

# Estimate the residuals with the AR(1) coefficient and sigma.
res.heter_45=mat.or.vec(heter_subset_length, length(years.rcp)) #(nr,nc)
heter.slr.rcp45=res.heter_45
for(n in 1:heter_subset_length) {
  for(i in 2:length(years.rcp)) {
    res.heter_45[n,i] = heter_sub_chain[n,5]*res.heter_45[n,i-1] +
      rnorm(1, mean = 0, sd = heter_sub_chain[n,4])
  }
}

# Estimate the projections: add the residuals onto the model simulations. Equation (2) & (S1)
for(i in 1:heter_subset_length) {
  heter.slr.rcp45[i,] = heter.rcp45.proj[i,] + res.heter_45[i,]
}

remove(res.heter_45, heter.rcp45.proj, heter.rcp45.rate)
# Heteroskedastic RCP60 ---------------------------------------------------------
heter.rcp60.rate=mat.or.vec(heter_subset_length, length(years.rcp)) #(nr,nc)
heter.rcp60.proj=mat.or.vec(heter_subset_length, length(years.rcp)) #(nr,nc)

for(n in 1:heter_subset_length) {
  # Estimate the sea level rate of change: equation (1)
  heter.rcp60.rate[n,] = heter_sub_chain[n,1]*(tempRCP60_anom - heter_sub_chain[n,2])
  heter.rcp60.proj[n,1] = heter_sub_chain[n,3]
  
  # Use Forward Euler to estimate sea level over time.
  for (i in from:to){
    heter.rcp60.proj[n,i] = heter.rcp60.proj[n,i-1] + heter.rcp60.rate[n,i-1]*timestep
  }
}

# Estimate the residuals with the AR(1) coefficient and sigma.
res.heter_60=mat.or.vec(heter_subset_length, length(years.rcp)) #(nr,nc)
for(n in 1:heter_subset_length) {
  for(i in 2:length(years.rcp)) {
    res.heter_60[n,i] = heter_sub_chain[n,5]*res.heter_60[n,i-1] +
      rnorm(1, mean = 0, sd = heter_sub_chain[n,4])
  }
}
heter.slr.rcp60 = res.heter_60

# Estimate the projections: add the residuals onto the model simulations. Equation (2) & (S1)
for(i in 1:heter_subset_length) {
  heter.slr.rcp60[i,] = heter.rcp60.proj[i,] + res.heter_60[i,]
}

remove(res.heter_60, heter.rcp60.proj, heter.rcp60.rate)
# Heteroskedastic RCP85 ---------------------------------------------------------
heter.rcp85.rate=mat.or.vec(heter_subset_length, length(years.rcp)) #(nr,nc)
heter.rcp85.proj=mat.or.vec(heter_subset_length, length(years.rcp)) #(nr,nc)

for(n in 1:heter_subset_length) {
  # Estimate the sea level rate of change: equation (1)
  heter.rcp85.rate[n,] = heter_sub_chain[n,1]*(tempRCP85_anom - heter_sub_chain[n,2])
  heter.rcp85.proj[n,1] = heter_sub_chain[n,3]
  
  # Use Forward Euler to estimate sea level over time.
  for (i in from:to){
    heter.rcp85.proj[n,i] = heter.rcp85.proj[n,i-1] + heter.rcp85.rate[n,i-1]*timestep
  }
}

# Estimate the residuals with the AR(1) coefficient and sigma.
res.heter_85=mat.or.vec(heter_subset_length, length(years.rcp)) #(nr,nc)
for(n in 1:heter_subset_length) {
  for(i in 2:length(years.rcp)) {
    res.heter_85[n,i] = heter_sub_chain[n,5]*res.heter_85[n,i-1] +
      rnorm(1, mean = 0, sd = heter_sub_chain[n,4])
  }
}

# Estimate the projections: add the residuals onto the model simulations. Equation (2) & (S1)
heter.slr.rcp85 = res.heter_85
for(i in 1:heter_subset_length) {
  heter.slr.rcp85[i,]=heter.rcp85.proj[i,] + res.heter_85[i,]
}

remove(res.heter_85, heter.rcp85.proj, heter.rcp85.rate)

###############################################################################
# Plot ------------------------------------------------------------------------
source("Scripts/put_fig_letter.R")
######################################## MAIN FIGURES #############################################

# Calculate the 99% confidence interval for each method
boot_rcp26.005 <-
  boot_rcp26.995 <-
  boot_rcp45.005 <-
  boot_rcp45.995 <-
  boot_rcp60.005 <-
  boot_rcp60.995 <-
  boot_rcp85.005 <-
  boot_rcp85.995 <-
  hom_rcp26.005 <-
  hom_rcp26.995 <-
  hom_rcp45.005 <-
  hom_rcp45.995 <-
  hom_rcp60.005 <-
  hom_rcp60.995 <-
  hom_rcp85.005 <-
  hom_rcp85.995 <-
  het_rcp26.005 <-
  het_rcp26.995 <-
  het_rcp45.005 <-
  het_rcp45.995 <-
  het_rcp60.005 <-
  het_rcp60.995 <-
  het_rcp85.005 <-
  het_rcp85.995 <-rep(NA,length(years.rcp)) # nyears.mod is the total number of data (221 years)
for(i in 1:length(years.rcp)){
  boot_rcp26.005[i] <-quantile(boot.slr.rcp26[,i]/100,0.005)
    boot_rcp26.995[i] <-quantile(boot.slr.rcp26[,i]/100,0.995)
    boot_rcp45.005[i] <-quantile(boot.slr.rcp45[,i]/100,0.005)
    boot_rcp45.995[i] <-quantile(boot.slr.rcp45[,i]/100,0.995)
    boot_rcp60.005[i] <-quantile(boot.slr.rcp60[,i]/100,0.005)
    boot_rcp60.995[i] <-quantile(boot.slr.rcp60[,i]/100,0.995)
    boot_rcp85.005[i] <-quantile(boot.slr.rcp85[,i]/100,0.005)
    boot_rcp85.995[i] <-quantile(boot.slr.rcp85[,i]/100,0.995)
    hom_rcp26.005[i] <-quantile(homo.slr.rcp26[,i]/100,0.005)
    hom_rcp26.995[i] <-quantile(homo.slr.rcp26[,i]/100,0.995)
    hom_rcp45.005[i] <-quantile(homo.slr.rcp45[,i]/100,0.005)
    hom_rcp45.995[i] <-quantile(homo.slr.rcp45[,i]/100,0.995)
    hom_rcp60.005[i] <-quantile(homo.slr.rcp60[,i]/100,0.005)
    hom_rcp60.995[i] <-quantile(homo.slr.rcp60[,i]/100,0.995)
    hom_rcp85.005[i] <-quantile(homo.slr.rcp85[,i]/100,0.005)
    hom_rcp85.995[i] <-quantile(homo.slr.rcp85[,i]/100,0.995)
    het_rcp26.005[i] <-quantile(heter.slr.rcp26[,i]/100,0.005)
    het_rcp26.995[i] <-quantile(heter.slr.rcp26[,i]/100,0.995)
    het_rcp45.005[i] <-quantile(heter.slr.rcp45[,i]/100,0.005)
    het_rcp45.995[i] <-quantile(heter.slr.rcp45[,i]/100,0.995)
    het_rcp60.005[i] <-quantile(heter.slr.rcp60[,i]/100,0.005)
    het_rcp60.995[i] <-quantile(heter.slr.rcp60[,i]/100,0.995)
    het_rcp85.005[i] <-quantile(heter.slr.rcp85[,i]/100,0.005)
    het_rcp85.995[i] <-quantile(heter.slr.rcp85[,i]/100,0.995)
}

#--------------- Estimate survival function
source("Scripts/plot_sf.r")
survival.boot.26 <- plot.sf(boot.slr.rcp26[,221]/100, make.plot=F)
survival.boot.45 <- plot.sf(boot.slr.rcp45[,221]/100, make.plot=F)
survival.boot.60 <- plot.sf(boot.slr.rcp60[,221]/100, make.plot=F)
survival.boot.85 <- plot.sf(boot.slr.rcp85[,221]/100, make.plot=F)

survival.homo.26 <- plot.sf(homo.slr.rcp26[,221]/100, make.plot=F)
survival.homo.45 <- plot.sf(homo.slr.rcp45[,221]/100, make.plot=F)
survival.homo.60 <- plot.sf(homo.slr.rcp60[,221]/100, make.plot=F)
survival.homo.85 <- plot.sf(homo.slr.rcp85[,221]/100, make.plot=F)

survival.het.26 <- plot.sf(heter.slr.rcp26[,221]/100, make.plot=F)
survival.het.45 <- plot.sf(heter.slr.rcp45[,221]/100, make.plot=F)
survival.het.60 <- plot.sf(heter.slr.rcp60[,221]/100, make.plot=F)
survival.het.85 <- plot.sf(heter.slr.rcp85[,221]/100, make.plot=F)

#Save the workspace:
save.image(file = "Workspace/rcp_temp_scenarios.RData")
# load("Workspace/rcp_temp_scenarios.RData")
#------------------------------------- Transparent Color Function -------------------------------------------
makeTransparent<- function(somecolor, alpha=100){
  someColor = someColor
  newColor<-col2rgb(someColor)
  apply(newColor,2 ,
        function(curcoldata)
        {rgb(red=curcoldata[1],
             green=curcoldata[2],
             blue=curcoldata[3], alpha=alpha,
             maxColorValue=255)})
}
library(RColorBrewer)
test.colors = brewer.pal(9, "YlGnBu")

# FOR PDF
# someColor = c("slateblue", "paleturquoise2", "gray")
# trans_con_colors = makeTransparent(someColor,150)
# someColor = "gray"
# trans_gray = makeTransparent(someColor,200)

# FOR EPS
trans_con_colors = c(test.colors[4], test.colors[6], test.colors[8])
trans_gray = "gray"

mm_TO_inches = function(mm){
  mm * 0.039370
}

single_column = mm_TO_inches(84)
double_column = mm_TO_inches(174)
maximum_width = mm_TO_inches(234)
column_height=2.7

#------------------------------------- Calculate 90% and 99% credible intervals -------------------------------------------
# pdf(file="nRuckertetal_supp5_tempscen.pdf", family="Helvetica", width=6.7, height=6.7, pointsize=11)
setEPS()
# postscript(file="SuppFigures/sfigure4a_4d.eps", horizontal = FALSE, onefile = FALSE, paper = "special", family="Helvetica", width=6.7, height=5.4, pointsize=11)
postscript(file="SuppFigures/SFig6.eps", horizontal = FALSE, onefile = FALSE, paper = "special", family="Helvetica",
width=double_column, height=column_height*2, pointsize=11)

par(mfrow=c(2,2), mgp=c(1.5,.5,0),  mar=c(3.5,4,1,1)) # set figure dimensions
y_99=c(years.rcp, rev(years.rcp))

# a. RCP8.5
plot(year, slr/100, type="l", xlab="",ylab="Sea-level anomalies [m]",
col=trans_gray,lwd=1,ylim=c(0,1.8), xlim=c(2000,2090))
#Method C
het_x_99=c(het_rcp85.005, rev(het_rcp85.995))
polygon(y_99, het_x_99, col=trans_con_colors[3], border=NA)
#Method B
hom_x_99=c(hom_rcp85.005, rev(hom_rcp85.995))
polygon(y_99, hom_x_99, col=trans_con_colors[2], border=NA)
#Method A
boot_x_99=c(boot_rcp85.005, rev(boot_rcp85.995))
polygon(y_99, boot_x_99, col=trans_con_colors[1], border=NA) # plot 90% polygon

# lines(scenario_time, a1fi_p$sle/100, lwd=1.5, col="salmon") # Add A2 SLR estimated
# text(2080,0.8, "A1FI", col="salmon3", cex=0.85)
# lines(scenario_time, a2_p$sle/100, lwd=1.5, col="salmon") # Add A2 SLR estimated
# text(2080,0.45, "A2", col="salmon3", cex=0.85)
axis(side=4, labels=FALSE)
put.fig.letter("a. RCP8.5",font=2,  x=0.15, y=0.98)

legend("topleft", c('99% CI Bootstrap','99% CI Bayesian (homoskedastic)','99% CI Bayesian (heteroskedastic)'), cex=0.85,
       pch=15, bty="n", col=trans_con_colors[1:3])

# b. RCP6.0
plot(year, slr/100, type="l", xlab="",ylab="",
col=trans_gray,lwd=1,ylim=c(0,1.8), xlim=c(2000,2090))
#Method C
het_x_99=c(het_rcp60.005, rev(het_rcp60.995))
polygon(y_99, het_x_99, col=trans_con_colors[3], border=NA)
#Method B
hom_x_99=c(hom_rcp60.005, rev(hom_rcp60.995))
polygon(y_99, hom_x_99, col=trans_con_colors[2], border=NA)
#Method A
boot_x_99=c(boot_rcp60.005, rev(boot_rcp60.995))
polygon(y_99, boot_x_99, col=trans_con_colors[1], border=NA) # plot 90% polygon

# lines(scenario_time, a1b_p$sle/100, lwd=1.5, col="salmon") # Add A1b SLR estimated
# text(2080,0.5, "A1B", col="salmon3", cex=0.85)
axis(side=4, labels=FALSE)
put.fig.letter("b. RCP6.0",font=2,  x=0.15, y=0.98)

# c. RCP4.5
plot(year, slr/100, type="l", xlab="Year",ylab="Sea-level anomalies [m]",
col=trans_gray,lwd=1,ylim=c(0,1.8), xlim=c(2000,2090))
#Method C
het_x_99=c(het_rcp45.005, rev(het_rcp45.995))
polygon(y_99, het_x_99, col=trans_con_colors[3], border=NA)
#Method B
hom_x_99=c(hom_rcp45.005, rev(hom_rcp45.995))
polygon(y_99, hom_x_99, col=trans_con_colors[2], border=NA)
#Method A
boot_x_99=c(boot_rcp45.005, rev(boot_rcp45.995))
polygon(y_99, boot_x_99, col=trans_con_colors[1], border=NA) # plot 90% polygon

# lines(scenario_time, b2_p$sle/100, lwd=1.5, col="salmon") # Add b1 SLR estimated
# text(2080,0.5, "B2", col="salmon3", cex=0.85)
axis(side=4, labels=FALSE)
put.fig.letter("c. RCP4.5",font=2,  x=0.15, y=0.98)

# d. RCP2.6
plot(year, slr/100, type="l", xlab="Year",ylab="",
col=trans_gray,lwd=1,ylim=c(0,1.8), xlim=c(2000,2090))
#Method C
het_x_99=c(het_rcp26.005, rev(het_rcp26.995))
polygon(y_99, het_x_99, col=trans_con_colors[3], border=NA)
#Method B
hom_x_99=c(hom_rcp26.005, rev(hom_rcp26.995))
polygon(y_99, hom_x_99, col=trans_con_colors[2], border=NA)
#Method A
boot_x_99=c(boot_rcp26.005, rev(boot_rcp26.995))
polygon(y_99, boot_x_99, col=trans_con_colors[1], border=NA) # plot 90% polygon

# lines(scenario_time, b1_p$sle/100, lwd=1.5, col="salmon") # Add b1 SLR estimated
# text(2080,0.45, "B1", col="salmon3", cex=0.85)
# lines(scenario_time, min_p$sle/100, lwd=2, col="salmon") # Add minimum SLR estimated
axis(side=4, labels=FALSE)
put.fig.letter("d. RCP2.6",font=2, x=0.15, y=0.98)

dev.off()

################################ survival function ###################################################
setEPS()
postscript(file="SuppFigures/SFig7.eps", horizontal = FALSE, onefile = FALSE, paper = "special", family="Helvetica",
width=double_column, height=column_height*2, pointsize=11)
par(mfrow=c(2,2), mgp=c(1.5,.5,0),  mar=c(3.5,4,1,1)) # set figure dimensions
# survival function for SLR values in 2100
plot(survival.het.85$sf.num, survival.het.85$sf, log="y", type="l",
col=trans_con_colors[3], ylab="Survival function
[1-cumulative probability]", lwd=3,
xlab="",main="",ylim=c(1e-04,1), yaxt="n")
lines(survival.boot.85$sf.num, survival.boot.85$sf, lwd=3, col=trans_con_colors[1])
lines(survival.homo.85$sf.num, survival.homo.85$sf, lwd=3, col=trans_con_colors[2])
abline(h=1/100, lty=2)
axis(2, at=10^(-4:0), label=parse(text=paste("10^", -4:0, sep="")))
text(0.5, 0.015, expression(paste("10"^"-2", sep="")))
legend("topright", c('Bootstrap','Bayesian (homoskedastic)','Bayesian (heteroskedastic)'), cex=0.85,
pch=15, bty="n",col=trans_con_colors[1:3])
put.fig.letter("a. RCP8.5",font=2,  x=0.15, y=0.98)

plot(survival.het.60$sf.num, survival.het.60$sf, log="y", type="l",
col=trans_con_colors[3], ylab="", lwd=3,
xlab="",main="",ylim=c(1e-04,1), yaxt="n")
lines(survival.boot.60$sf.num, survival.boot.60$sf, lwd=3, col=trans_con_colors[1])
lines(survival.homo.60$sf.num, survival.homo.60$sf, lwd=3, col=trans_con_colors[2])
abline(h=1/100, lty=2)
axis(2, at=10^(-4:0), label=parse(text=paste("10^", -4:0, sep="")))
put.fig.letter("b. RCP6.0",font=2,  x=0.15, y=0.98)

plot(survival.het.45$sf.num, survival.het.45$sf, log="y", type="l",
col=trans_con_colors[3], ylab="Survival function
[1-cumulative probability]", lwd=3,
xlab="Projected sea-level 2100 [m]",main="",ylim=c(1e-04,1), yaxt="n")
lines(survival.boot.45$sf.num, survival.boot.45$sf, lwd=3, col=trans_con_colors[1])
lines(survival.homo.45$sf.num, survival.homo.45$sf, lwd=3, col=trans_con_colors[2])
abline(h=1/100, lty=2)
axis(2, at=10^(-4:0), label=parse(text=paste("10^", -4:0, sep="")))
put.fig.letter("c. RCP4.5",font=2,  x=0.15, y=0.98)

plot(survival.het.26$sf.num, survival.het.26$sf, log="y", type="l",
col=trans_con_colors[3], ylab="", lwd=3,
xlab="Projected sea-level 2100 [m]",main="",ylim=c(1e-04,1), yaxt="n")
lines(survival.boot.26$sf.num, survival.boot.26$sf, lwd=3, col=trans_con_colors[1])
lines(survival.homo.26$sf.num, survival.homo.26$sf, lwd=3, col=trans_con_colors[2])
abline(h=1/100, lty=2)
axis(2, at=10^(-4:0), label=parse(text=paste("10^", -4:0, sep="")))
put.fig.letter("d. RCP2.6",font=2, x=0.15, y=0.98)
dev.off()

################################ END ###################################################

