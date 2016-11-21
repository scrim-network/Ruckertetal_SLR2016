#--------- Markov Chain Monte Carlo DAIS Simulations
# file ~ DAIS_convergence_plots.R
########### Plot the results for further analysis on convergence #################################

################################## CONVERGENCE ####################################

# Remove "burn-in"
results <- prechain1[length(burnin):NI,]

################################## TRACE PLOTS ####################################
## 1) - 8) Plot the Trace and Histograms for each parameter in the chain to insure convergence:
jpeg(file="Converge_alpha.jpeg", family="Helvetica", width=1200, height=700, units="px", pointsize=20)
par(mfrow=c(2,1), mar=c(3, 7, 1, 7), mgp=c(1.5,.5,0))
plot(results[,1], type="l",
     ylab="alpha",
     xlab="Number of Runs", main="")
par(mar=c(3,7,1,7))
hist(results[,1], freq=FALSE, col="gray",
     xlab="alpha",
     ylab="Density [Dimensionless]", main="")
dev.off()

# ------------------------------------------------------------------------
jpeg(file="Converge_beta.jpeg", family="Helvetica", width=1200, height=700, units="px", pointsize=20)
par(mfrow=c(2,1), mar=c(3, 7, 1, 7), mgp=c(1.5,.5,0))
plot(results[,2], type="l",
     ylab="beta",
     xlab="Number of Runs", main="")
par(mar=c(3,7,1,7))
hist(results[,2], freq=FALSE, col="gray",
     xlab="beta",
     ylab="Density [Dimensionless]", main="")
dev.off()

# ------------------------------------------------------------------------
jpeg(file="Converge_S_1.jpeg", family="Helvetica", width=1200, height=700, units="px", pointsize=20)
par(mfrow=c(2,1), mar=c(3, 7, 1, 7), mgp=c(1.5,.5,0))
plot(results[,3], type="l",
     ylab="initial sea-level value",
     xlab="Number of Runs", main="")
par(mar=c(3,7,1,7))
hist(results[,3], freq=FALSE, col="gray",
     xlab="initial sea-level value",
     ylab="Density [Dimensionless]", main="")
dev.off()

# ------------------------------------------------------------------------
jpeg(file="Converge_tau.jpeg", family="Helvetica", width=1200, height=700, units="px", pointsize=20)
par(mfrow=c(2,1), mar=c(3, 7, 1, 7), mgp=c(1.5,.5,0))
plot(results[,4], type="l",
     ylab="tau",
     xlab="Number of Runs", main="")
par(mar=c(3,7,1,7))
hist(results[,4], freq=FALSE, col="gray",
     xlab="tau",
     ylab="Density [Dimensionless]", main="")
dev.off()

# ------------------------------------------------------------------------
jpeg(file="Converge_sigma.jpeg", family="Helvetica", width=1200, height=700, units="px", pointsize=20)
par(mfrow=c(2,1), mar=c(3, 7, 1, 7), mgp=c(1.5,.5,0))
plot(results[,5], type="l",
     ylab="sigma",
     xlab="Number of Runs", main="")
par(mar=c(3,7,1,7))
hist(results[,5], freq=FALSE, col="gray",
     xlab="sigma",
     ylab="Density [Dimensionless]", main="")
dev.off()

# ------------------------------------------------------------------------
jpeg(file="Converge_rho.jpeg", family="Helvetica", width=1200, height=700, units="px", pointsize=20)
par(mfrow=c(2,1), mar=c(3, 7, 1, 7), mgp=c(1.5,.5,0))
plot(results[,6], type="l",
     ylab="rho",
     xlab="Number of Runs", main="")
par(mar=c(3,7,1,7))
hist(results[,6], freq=FALSE, col="gray",
     xlab="rho",
     ylab="Density [Dimensionless]", main="")
dev.off()

#------------------------------ END ---------------------------------------------------------------------


