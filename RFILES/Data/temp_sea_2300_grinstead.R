########################################################################################
#
#  -file = "temp_sea_2300_grinsted.R"   Code written April 2014, updated Sept. 2016
#  - Author: Kelsey Ruckert (klr324@psu.edu)
#
#  -This program loads in temperature and global sea-level data for use in model
#       and uncertainty calculations described in Ruckert et al. (accepted). For
#       further description and references, please read the paper.
#
# THIS CODE IS PROVIDED AS-IS WITH NO WARRANTY (NEITHER EXPLICIT
# NOR IMPLICIT).  I SHARE THIS CODE IN HOPES THAT IT IS USEFUL,
# BUT I AM NOT LIABLE FOR THE BEHAVIOR OF THIS CODE IN YOUR OWN
# APPLICATION.  YOU ARE FREE TO SHARE THIS CODE SO LONG AS THE
# AUTHOR(S) AND VERSION HISTORY REMAIN INTACT.
#
#   -NOTE: This file contains data that is sourced into the other programs. Information
#       regarding this data can be found below:
#
#       -RCP8.5 is used to create temperature simulations to 2300
#       -RCP8.5 simulates "Business as usual" and is similar to the A2 scenario
#       -The RCP8.5 temperatures are from the CNRM-CM5 model
#           (Centre National de Recherches Meteorologiques)
#       -These RCP8.5 temperatures closly resemble historical temperatures
#           from the National Oceanic and Atmospheric Administration (NOAA)
#
#       -Annual global land and ocean temperature anomalies (C)
#       -Anomalies with respect to the 20th century average
#       -http://www.ncdc.noaa.gov/cag/time-series/global
#
#       -Tide guage data from Church & White_GRL_2006
#       -http://www.psmsl.org/products/reconstructions/church.php
#       -SLR values are in relation to the 1990 mean value
#
########################################################################################

# Read in the information from Smith (2008) historical temperatures and the RCP8.5 scenario.
# alltime, time in years from 1880 to 2300 (5yr increments)
# hist.temp, historical global mean temperatures from 1880 to 2002 in C
# rcp85, merged historical + rcp 8.5 temperatures from 1880 to 2300 in C
# scenario_time, , time in years from 1880 to 2100 (5yr increments)
# max to b2, merged historical + IPCC temperature scenarios from 1880-2100 in C (5yr increments)
temp.data = read.csv("../Data/NOAA_IPCC_RCPtempsscenarios.csv")
alltime = temp.data[, 1]
hist.temp = temp.data[1:122, 2]
rcp85 = temp.data[, 4]
scenario_time = temp.data[1:45, 5]
max = temp.data[1:45, 6]
min = temp.data[1:45, 7]
a1fi = temp.data[1:45, 8]
a1b = temp.data[1:45, 9]
a1t = temp.data[1:45, 10]
a2 = temp.data[1:45, 11]
b1 = temp.data[1:45, 12]
b2 = temp.data[1:45, 13]

# Read in and extract three key vectors from the Church and White (2006) sea-level file.
# year, time in years from 1880 to 2001, 1yr increments
# slr, global mean sea level in mm
# err.obs, global mean sea level measurement error in mm
# Divide vectors by 1000 to convert to m
church = read.table("../Data/church_13221.txt")
year = church[11:132, 1]
slr = church[11:132, 2]/1000
err.obs = church[11:132,3]/1000

## Estimate sea-level anomalies with respect to the 1980-2000 period mean sea-level value (matching Grinsted (2010))
SLR1980_2000 = mean(slr[100:120])
slr = slr - SLR1980_2000

# Calculate sea level -/+ observation errors.
err_pos=slr+err.obs
err_neg=slr-err.obs

######################################## END ###############################################
