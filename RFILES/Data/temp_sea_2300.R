########################################################################################
#
#  -file = "temp_sea_2300.R"   Code written April 2014, updated July 2014
#  - Author: Kelsey Ruckert (klr324@psu.edu)
#
#  -This program loads in temperature and global sea-level data for use in model
#       and uncertainty calculations described in Ruckert et al. (2016). For
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

#Read in temperature and sea level data
data = read.csv("Data/NOAA_IPCC_RCPtempsscenarios.csv")

#Historical time frame and temperatures from NOAA
hist.temp = data[1:122,2] #temperature data
alltime = data[,1] #1880-2300 (1 yr incriments)

#IPCC temperature scenarios added to NOAA temperatures from 1880-2100
scenario_time = data[1:45,5] #1880-2100 (5yr incriments)
max = data[1:45,6]  #6.195 C in 2100
min = data[1:45,7]  #1.773 C in 2100
a1fi = data[1:45,8] #4.892 C in 2100
a1b = data[1:45,9]  #3.347 C in 2100
a1t = data[1:45,10] #2.938 C in 2100
a2 = data[1:45,11]  #4.194 C in 2100
b1 = data[1:45,12]  #2.382 C in 2100
b2 = data[1:45,13]  #3.087 C in 2100

#RCP 8.5 temperatures from 2014-2300 and NOAA historical temps from 1880-2013
rcp85 = data[,4]    #4.284 C in 2100 and 9.849284058 C in 2300

# Historical global mean sea-levels from tide gauges & estimated errors
church = read.table("Data/church_13221.txt")
year = church[11:132, 1] #timeframe 1880 to 2002 in 1yre incriments
slr = church[11:132, 2]/10 #mm to cm
err.obs = church[11:132,3]/10 #mm to cm

## To match Rahmstorf (2007) Sea-leve values are set in reference to the 1990 mean SLR value
SLR1990 = slr[111]  #111 equals the year 1990
slr = slr - SLR1990

# Set the observational errors by adding and subtracting the errors to the sea-level values
err_pos=slr+err.obs
err_neg=slr-err.obs

######################################## END ###############################################
