
#######################################################################
#
# Estimate_surprise_indicies.R    15 Sept. 2015
#
# Author: Kelsey Ruckert (klr324@psu.edu)
#
# Script with function that helps to create the reliability diagram
#
# Note: This script/function is used in part to help create the
# MODELNAME_surprise.csv files in Ruckert et al. (accepted). See
# section 2 and equations S13-14 for more details on testing
# hindcast reliability.
#
# THIS CODE IS PROVIDED AS-IS WITH NO WARRANTY (NEITHER EXPLICIT
# NOR IMPLICIT).  I SHARE THIS CODE IN HOPES THAT IT IS USEFUL,
# BUT I AM NOT LIABLE FOR THE BEHAVIOR OF THIS CODE IN YOUR OWN
# APPLICATION.  YOU ARE FREE TO SHARE THIS CODE SO LONG AS THE
# AUTHOR(S) AND VERSION HISTORY REMAIN INTACT.
#
# Function Name: outliers
# Parameters:
#   slr - vector of data points
#   lower.quant - vector of the lower values in the credible interval
#   upper.quant - vector of the upper values in the credible interval
#
# The function returns the number of outliers (the number of data points
# that fall outside of a specified credible/confidence interval).
#
#######################################################################
# Choose method -----------------------------------------------------------

# If working with Bootstrap:
#load("Workspace/Bootstrap_Grinsted_SLR_model_workspace.RData")
#load("Workspace/Bootstrap_VR_SLR_model_workspace.RData")
hindcast = SLR.projections.boot


# If working with MCMC homoskedastic:
#load("Workspace/Homoskedastic_Grinsted_SLR_model_workspace.RData")
#load("Workspace/Homoskedastic_vermeer_SLR_model_workspace.RData")
hindcast = SLR.projections.homo


# If working with MCMC heteroskedastic:
#load("Workspace/Heteroskedastic_Grinsted_SLR_model_workspace.RData")
load("Workspace/Heteroskedastic_vermeer_SLR_model_workspace.RData")
hindcast = SLR.projections.heter


# Estimate quantiles ------------------------------------------------------
datalength = nyears.mod

ar.005 <-
  ar.995 <- rep(NA,datalength) # 99%
ar.45 <- ar.995; ar.55 <- ar.995 # 10%
ar.40 <- ar.995; ar.60 <- ar.995 # 20%
ar.35 <- ar.995; ar.65 <- ar.995 # 30%
ar.30 <- ar.995; ar.70 <- ar.995 # 40%
ar.25 <- ar.995; ar.75 <- ar.995 # 50%
ar.20 <- ar.995; ar.80 <- ar.995 # 60%
ar.15 <- ar.995; ar.85 <- ar.995 # 70%
ar.10 <- ar.995; ar.90 <- ar.995 # 80%
ar.5 <- ar.995; ar.95 <- ar.995  # 90%
ar.4 <- ar.995; ar.96 <- ar.995  # 92%
ar.025 <- ar.995; ar.975 <- ar.995 # 95%
ar.02 <- ar.995; ar.98 <- ar.995   # 96%
ar.015 <- ar.995; ar.985 <- ar.995 # 97%
ar.01 <- ar.995; ar.99 <- ar.995   # 98%
ar.0 <- ar.995; ar.100 <- ar.995   # 98%

for(i in 1:datalength){
  ar.45[i] = quantile(hindcast[,i],0.45); ar.55[i] = quantile(hindcast[,i],0.55)
  ar.40[i] = quantile(hindcast[,i],0.40); ar.60[i] = quantile(hindcast[,i],0.60)
  ar.35[i] = quantile(hindcast[,i],0.35); ar.65[i] = quantile(hindcast[,i],0.65)
  ar.30[i] = quantile(hindcast[,i],0.30); ar.70[i] = quantile(hindcast[,i],0.70)
  ar.25[i] = quantile(hindcast[,i],0.25); ar.75[i] = quantile(hindcast[,i],0.75)
  ar.20[i] = quantile(hindcast[,i],0.20); ar.80[i] = quantile(hindcast[,i],0.80)
  ar.15[i] = quantile(hindcast[,i],0.15); ar.85[i] = quantile(hindcast[,i],0.85)
  ar.10[i] = quantile(hindcast[,i],0.10); ar.90[i] = quantile(hindcast[,i],0.90)
  ar.5[i] = quantile(hindcast[,i],0.05); ar.95[i] = quantile(hindcast[,i],0.95)
  ar.4[i] = quantile(hindcast[,i],0.04); ar.96[i] = quantile(hindcast[,i],0.96)
  ar.025[i] = quantile(hindcast[,i],0.025); ar.975[i] = quantile(hindcast[,i],0.975)
  ar.02[i] = quantile(hindcast[,i],0.02); ar.98[i] = quantile(hindcast[,i],0.98)
  ar.015[i] = quantile(hindcast[,i],0.015); ar.985[i] = quantile(hindcast[,i],0.985)
  ar.01[i] = quantile(hindcast[,i],0.01); ar.99[i] = quantile(hindcast[,i],0.99)
  ar.005[i] = quantile(hindcast[,i],0.005); ar.995[i] = quantile(hindcast[,i],0.995)
  ar.0[i] = quantile(hindcast[,i],0); ar.100[i] = quantile(hindcast[,i],1)
}

range_y=c(alltime, rev(alltime))
range_ar_10=c(ar.45, rev(ar.55)); range_ar_20=c(ar.40, rev(ar.60))
range_ar_30=c(ar.35, rev(ar.65)); range_ar_40=c(ar.30, rev(ar.70))
range_ar_50=c(ar.25, rev(ar.75)); range_ar_60=c(ar.20, rev(ar.80))
range_ar_70=c(ar.15, rev(ar.85)); range_ar_80=c(ar.10, rev(ar.90))
range_ar_90=c(ar.5, rev(ar.95)); range_ar_92=c(ar.4, rev(ar.96))
range_ar_95=c(ar.025, rev(ar.975)); range_ar_96=c(ar.02, rev(ar.98))
range_ar_97=c(ar.015, rev(ar.985)); range_ar_98=c(ar.01, rev(ar.99))
range_ar_99=c(ar.005, rev(ar.995)); range_ar_100=c(ar.0, rev(ar.100))

# Estimate number of outliers
outliers = function(slr, lower.quant, upper.quant){
  in.range = rep(0,length(slr))
  
  for(i in 1:length(slr)){
    in.range[i] = (slr[i] > lower.quant[i]) & (slr[i] < upper.quant[i])
  }
  length(slr) - sum(in.range)
}

outliers(slr, ar.45, ar.55)
outliers(slr, ar.40, ar.60)
outliers(slr, ar.35, ar.65)
outliers(slr, ar.30, ar.70)
outliers(slr, ar.25, ar.75)
outliers(slr, ar.20, ar.80)
outliers(slr, ar.15, ar.85)
outliers(slr, ar.10, ar.90)
outliers(slr, ar.5, ar.95)
outliers(slr, ar.4, ar.96)
outliers(slr, ar.025, ar.975)
outliers(slr, ar.02, ar.98)
outliers(slr, ar.015, ar.985)
outliers(slr, ar.01, ar.99)
outliers(slr, ar.005, ar.995)
outliers(slr, ar.0, ar.100)

# Plot --------------------------------------------------------------------
plot(year, slr, type="l",xlab="x",ylab="Observations", ylim=c(-0.16,0.07))
points(year, slr, pch=20, cex=0.85)
polygon(range_y, range_ar_10, col="red")
polygon(range_y, range_ar_20, col="red")
polygon(range_y, range_ar_30, col="red")
polygon(range_y, range_ar_40, col="red")
polygon(range_y, range_ar_50, col="red")
polygon(range_y, range_ar_60, col="red")
polygon(range_y, range_ar_70, col="red")
polygon(range_y, range_ar_80, col="red")
polygon(range_y, range_ar_90, col="red")
polygon(range_y, range_ar_92, col="red")
polygon(range_y, range_ar_95, col="red")
polygon(range_y, range_ar_96, col="red")
polygon(range_y, range_ar_97, col="red")
polygon(range_y, range_ar_98, col="red")
polygon(range_y, range_ar_99, col="red")
polygon(range_y, range_ar_100, col="red")

# End ---------------------------------------------------------------------













