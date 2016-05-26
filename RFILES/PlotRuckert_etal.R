##################################################################################################
#
#  -file = "PlotRuckert_etal.R"   Code written September 2014
#  - Author: Kelsey Ruckert (klr324@psu.edu)
#
#  -This program loads in the workspace generated from the file Mega_Rahmstorf.R
#       and produces the graphs shown in Ruckert et al. (2016). It also produces
#       the graphs shown in that paper. For further description and references,
#       please read the paper.
#
#   -NOTE: The graphs will be saved as pdf files in the current working directory.
#
# THIS CODE IS PROVIDED AS-IS WITH NO WARRANTY (NEITHER EXPLICIT
# NOT IMPLICIT).  I SHARE THIS CODE IN HOPES THAT IT IS USEFUL,
# BUT I AM NOT LIABLE FOR THE BEHAVIOR OF THIS CODE IN YOUR OWN
# APPLICATION.  YOU ARE FREE TO SHARE THIS CODE SO LONG AS THE
# AUTHOR(S) AND VERSION HISTORY REMAIN INTACT.
#
#----------------------------- GRAPH DETAILS -----------------------------------------------------
#  -Figure 1: 90% Confidence interval of sea-level rise both as a hindcast and
#        projection out to 2050
#
#  -Figure 2: The surprise index associated with each method producing uncertainty and
#        a zoomed in plot of the surprise index above 90%
#
#  -Figure 3: Probability density functions of each of the parameters estimated
#
#  -Figure 4: A probability density function of sea-level rise in the year 2050. Followed by
#       the cumulative desity function and the surival function of sea-level rise estimated
#       in each method for the year 2050
#
#  -Supplementary Figure 1: A hindcast plot showing the best fit estimated from each method
#       including the best fit generated in Rahmstorf (2007)
#
#  -Supplementary Figure 2: Trend of observation errors (heteroskedastic errors) from the sea-level record.   
#         The right panel plots the autocorrelation coefficient showing its significance for sea-level.
#
#  -Supplementary Figure 3: 90% Confidence interval of sea-level rise projected out to 2100.
#        The bottom panel plots the 99% credible interval of sea-level estimates projected to 2100
#
#  -Supplementary Figure 7: A probability density function of sea-level rise in the year 2100.
#       Followed by the cumulative desity function and the surival function of sea-level rise
#       estimated in each method for the year 2100
#
###################################################################################################

load("Workspace/mega_R_methods_workspace.RData") # Load in the saved workspace from file Mega_Rahmstorf.R

source("Scripts/put_fig_letter.R")
######################################## MAIN FIGURES #############################################

# Calculate the 90% confidence interval for each method
het_5 <-
het_95 <-
hom_5 <-
hom_95 <-
boot_5 <-
boot_95 <-
het_point5 <-
het_995 <-
hom_point5 <-
hom_995 <-
boot_point5 <-
boot_995 <- rep(NA,nyears.mod) # nyears.mod is the total number of data (421 years)
for(i in 1:nyears.mod){
    het_5[i] = quantile(slr.mcmc_proj[,i]/100,0.05) #heteroskedastic SLR values 90%
    het_95[i] = quantile(slr.mcmc_proj[,i]/100,0.95)
    hom_5[i] = quantile(h.slr.mcmc_proj[,i]/100,0.05) #homostedastic SLR values 90%
    hom_95[i] = quantile(h.slr.mcmc_proj[,i]/100,0.95)
    boot_5[i] = quantile(slr.boot_proj[,i]/100,0.05) # Bootstrap SLR values 90%
    boot_95[i] = quantile(slr.boot_proj[,i]/100,0.95)
    het_point5[i] = quantile(slr.mcmc_proj[,i]/100,0.005) #heteroskedastic SLR values 99%
    het_995[i] = quantile(slr.mcmc_proj[,i]/100,0.995)
    hom_point5[i] = quantile(h.slr.mcmc_proj[,i]/100,0.005) #homostedastic SLR values 99%
    hom_995[i] = quantile(h.slr.mcmc_proj[,i]/100,0.995)
    boot_point5[i] = quantile(slr.boot_proj[,i]/100,0.005) # Bootstrap SLR values 99%
    boot_995[i] = quantile(slr.boot_proj[,i]/100,0.995)
}

#---------------------------------- Figure 1 ----------------------------------------------
# 1) 90% Confidence interval of sea-level rise both as a hindcast + projection out to 2050:
#pdf(file="nRuckertetal1_revise.pdf", family="Helvetica", width=6.7, height=2.7, pointsize=9)
setEPS()
postscript(file="Figures/figure1a_1b.eps", horizontal = FALSE, onefile = FALSE, paper = "special", family="Helvetica", width=6.7, height=2.7, pointsize=9)

# Hindcast
par(mfrow=c(1,2), mgp=c(1.5,.5,0), mar=c(3.5,4,1,1))
# par(mfrow=c(1,2), mgp=c(1.5,.5,0), mar=c(4, 4, 3, 1)) # set figure dimensions
plot(year, slr/100, type="l",xlab="Year",ylab="Sea-level anomalies [m]", 
    ylim=c(-0.16,0.07), xlim=c(1884,1998))

# MCMC heteroskedastic#
het_x_90=c(het_5, rev(het_95)); het_y_90=c(years.mod, rev(years.mod))
polygon(het_y_90, het_x_90, col="slateblue", border=NA)
# MCMC Homoskedastic#
hom_x_90=c(hom_5, rev(hom_95)); hom_y_90=c(years.mod, rev(years.mod))
polygon(hom_y_90, hom_x_90, col="paleturquoise2", border=NA)
# Method A#
boot_x_90=c(boot_5, rev(boot_95)); boot_y_90=c(years.mod, rev(years.mod)) # set up the polygon
polygon(boot_y_90, boot_x_90, col="gray", border=NA) # plot 90% polygon
lines(scenario_time, a2_p$sle/100, lwd=2, col="salmon") # Add A2 SLR estimated
points(year, slr/100, pch=20, cex=0.5, col="black") # add in the observations
arrows(year, err_pos/100, year, err_neg/100, length=0)  # add in the observational errors
axis(side=4, labels=FALSE)
put.fig.letter("a.",font=2)
legend("topleft", c("Rahmstorf [2007] A2 estimate",'90% CI Method A',
                    '90% CI Method B','90% CI Method C', "Observations"), cex=0.85,
       pch=c(NA,15,15,15,20), lty=c(1,NA,NA,NA,NA),bty="n", lwd=2, 
       col=c("salmon","gray","paleturquoise2","slateblue","black"))
legend(1880, -0.005, "Observation error",pch="l", bty="n", col="black", cex=0.85)

#Projection to 2050
#par(mgp=c(1.5,.5,0), mar=c(4, 3, 3, 2)) # set figure dimensions
plot(year, slr/100, type="l", xlab="Year",ylab="Sea-level anomalies [m]",
     col="gray",lwd=1,ylim=c(0,0.4), xlim=c(2000,2050))

# MCMC Heteroskedastic:
polygon(het_y_90, het_x_90, col="slateblue", border=NA)
# MCMC Homoskedastic:
polygon(hom_y_90, hom_x_90, col="paleturquoise2", border=NA)
# Bootstrap:
polygon(boot_y_90, boot_x_90, col="gray", border=NA)
lines(scenario_time, a2_p$sle/100, lwd=2, col="salmon") # Add A2 SLR estimated
axis(side=4, labels=FALSE)
put.fig.letter("b.", font=2)
dev.off()

#---------------------------------- Figure 2 ----------------------------------------------
# 2) The surprise index associated with each method producing uncertainty and
#        a zoomed in plot of the surprise index above 90%
#
# Load in previously calculated surprise index information
surprise = read.csv("Data/rahm_surprise.csv")
percent = surprise[1:16,1]
boots = surprise[1:16,3]
homos = surprise[1:16,5]
heters = surprise[1:16,7]
boot_deviation = surprise[1:16,9]
homo_deviation = surprise[1:16,10]
heter_deviation = surprise[1:16,11]

# Full surprise index
# pdf(file="nRuckertetal2b_revise.pdf", family="Helvetica", width=3.5, height=8.1, pointsize=13)
# par(mfrow=c(3,1), mgp=c(1.5,.5,0),mar=c(2.5,4,4,2))

setEPS()
postscript(file="Figures/figure2a_2c.eps", horizontal = FALSE, onefile = FALSE, paper = "special", family="Helvetica", width=3.5, height=8.1, pointsize=13)
par(mfrow=c(3,1), mgp=c(1.5,.5,0), mar=c(3.5,4,1,1))

plot(percent, percent, typ="l", lty=2, ylab="Percent covered [%]", 
     xlab="Credible interval [%]", lwd=2)
points(percent, boots*100, pch=10, col="gray", lwd=2)
points(percent, homos*100, pch=8, col="paleturquoise2", lwd=2)
points(percent, heters*100, pch=4, col="slateblue", lwd=2)
text(50,65, "underconfident", srt=35)
text(55,45, "overconfident", srt=35)

put.fig.letter("a.",font=2)

# set up box in the upper corner to show where panel b is located
lines(c(88,102),c(102,102))
lines(c(88,88),c(88,102))
lines(c(102,102),c(88,102))
lines(c(88,102),c(88,88))
text(100,92, font=2, "b")

# zoomed plot of the section above 90%
#par(mgp=c(1.5,.5,0), mar=c(3.75,4,2.75,2)) 
plot(percent, percent, type="l", lty=2, ylab="Percent covered [%]", 
     xlab="Credible interval [%]", lwd=2, xlim=c(90,100), ylim=c(90,100))
points(percent, boots*100, pch=10, col="gray", lwd=2)
points(percent, homos*100, pch=8, col="paleturquoise2", lwd=2)
points(percent, heters*100, pch=4, col="slateblue", lwd=2)
text(94,95, "underconfident", srt=35)
text(95,94, "overconfident", srt=35)

put.fig.letter("b.",font=2)

legend("bottomright", c("Method A","Method B","Method C","Perfect 1:1 line"), pch=c(10,8,4,NA),
       lty=c(NA,NA,NA,2), lwd=2,col=c("gray","paleturquoise2","slateblue","black"))

homo_100=which(homos == 1)
heter_100=which(heters == 1)
boot_100=which(boots == 1)

#par(mgp=c(1.5,.5,0), mar=c(5,4,2.5,2))
plot(percent[1:homo_100[1]-1], homo_deviation[1:homo_100[1]-1], typ="l", col="paleturquoise2", lwd=2, ylim=c(-0.04,0.08),
     xlim=c(10,100), ylab="Surprise index [Dimensionless] ", xlab="Credible interval [%]")
points(percent[1:homo_100[1]-1], homo_deviation[1:homo_100[1]-1], pch=8, col="paleturquoise2", lwd=2)
lines(percent[1:boot_100[1]-1], boot_deviation[1:boot_100[1]-1], col="gray", lwd=2)
points(percent[1:boot_100[1]-1], boot_deviation[1:boot_100[1]-1], pch=10, col="gray", lwd=2)
lines(percent[1:heter_100[1]-1], heter_deviation[1:heter_100[1]-1], col="slateblue", lwd=2)
points(percent[1:heter_100[1]-1], heter_deviation[1:heter_100[1]-1], pch=4, col="slateblue", lwd=2)

abline(h=0, lty=2, lwd=2)
text(30,0.02, "underconfident")
text(30,-0.017, "overconfident")

put.fig.letter("c.",font=2)
dev.off()
#---------------------------------- Figure 3 ----------------------------------------------
# 3) Probability density functions of each of the parameters estimated

# pdf(file="nRuckertetal3_revise.pdf", family="Helvetica", width=3.5, height=8.1, pointsize=13)
# par(mfrow=c(3,1), mgp=c(1.5,.5,0),mar=c(2.5,4,4,2))
setEPS()
postscript(file="Figures/figure3a_3c.eps", horizontal = FALSE, onefile = FALSE, paper = "special", family="Helvetica", width=3.5, height=8.1, pointsize=13)
par(mfrow=c(3,1), mgp=c(1.5,.5,0), mar=c(3.5,4,1,1))

# Sea-level sensitivity (alpha)
plot(b.pdfalpha, main="", lwd=3, col="gray",xlab=expression(paste(alpha, " [cm/yr/", degree,"C]")),# "Sensitivity of SLR (a) [cm/yr/C]",
ylab="Probability Density", xlim=c(-0.01,2), yaxt="n") #homoskedastic
lines(h.pdfa, lwd=3, col="paleturquoise3") #bootstrap
lines(pdfa, lwd=3, col="slateblue")  #heteroskedastic
lines(c(0,0),c(0,1), lwd=3, lty=2)  # add in a lines to create the prior distribution
lines(c(2,2),c(0,1), lwd=3, lty=2)
lines(c(0,2),c(1,1), lwd=3, lty=2)
abline(v=original[1], lwd=3, col="salmon") # add in the original rahmstorf (2007) estimate
put.fig.letter("a.",font=2)
legend("topright", c("Rahmstorf [2007]
parameter estimate","Method A", "Method B","Method C",
                     "Prior"), lty=c(1,1,1,1,2), lwd=2,
       col=c("salmon","gray","paleturquoise3","slateblue","black"), bty="n")

#Base Temperature (T0)
#par(mgp=c(1.5,.5,0),mar=c(3.75,4,2.75,2))
plot(b.pdfbasetemp, lwd=3, col="gray", main="",xlab=expression(paste(T[0], " [",degree,"C]")),
     ylab="Probability Density", xlim=c(-3.3,2.1), yaxt="n") #homoskedastic
lines(h.pdfTo, lwd=3, col="paleturquoise3") #bootstrap
lines(pdfTo, lwd=3, col="slateblue") #heteroskedastic
lines(c(-3,-3),c(0,0.5), lwd=3, lty=2) # add in a lines to create the prior distribution
lines(c(2,2),c(0,0.5), lwd=3, lty=2)
lines(c(-3,2),c(0.5,0.5), lwd=3, lty=2)
abline(v=original[2], lwd=3, col="salmon") # add in the original rahmstorf (2007) estimate
put.fig.letter("b.",font=2)

#Initial value of sea-level in 1880
#par(mgp=c(1.5,.5,0), mar=c(5,4,2.5,2))
plot(b.pdfinitialval, lwd=3, col="gray", main="",xlab=expression(paste(H[0], " [cm]")),
     ylab="Probability Density", xlim=c(-17.3,-12), yaxt="n") #homoskedastic
lines(h.pdfinitialval, lwd=3, col="paleturquoise3") #bootstrap
lines(pdfinitialval, lwd=3, col="slateblue") #heteroskedastic
lines(c(-16.9,-16.9),c(0,.25), lwd=3, lty=2)  # add in a lines to create the prior distribution
lines(c(-12.4,-12.4),c(0,.25), lwd=3, lty=2)
lines(c(-16.9,-12.4),c(.25,.25), lwd=3, lty=2)
abline(v=original[3], lwd=3, col="salmon") # add in the original rahmstorf (2007) estimate
put.fig.letter("c.",font=2)
dev.off()

#---------------------------------- Figure 4 ----------------------------------------------
# 4) A probability density function of sea-level rise in the year 2050. Followed by
#       the cumulative desity function and the surival function of sea-level rise estimated
#       in each method for the year 2050

# pdf(file="nRuckertetal4.pdf", family="Helvetica", width=3.5, height=8.1, pointsize=13)
# par(mfrow=c(3,1), mgp=c(1.5,.5,0),mar=c(2.5,4,4,2))
setEPS()
postscript(file="Figures/figure4a_4c.eps", horizontal = FALSE, onefile = FALSE, paper = "special", family="Helvetica", width=3.5, height=8.1, pointsize=13)
par(mfrow=c(3,1), mgp=c(1.5,.5,0), mar=c(3.5,4,1,1))

# probability density function for SLR values in 2050
plot(b.pdf2050, main="", col="gray", lwd=3, ylab="Probabilty density", 
xlab="Projected sea-level 2050 [m]", xlim=c(-0.2,1), yaxt="n") #homoskedastic SLR values
lines(h.pdf2050, lwd=3, col="paleturquoise3") #bootstrap SLR values
lines(pdf2050, lwd=3, col="slateblue") #heteroskedastic SLR values
rahm_50 = c(min_p$sle[35]/100,max_p$sle[35]/100) # set up the range generated in rahmstorf (2007) for 2050
lines(rahm_50, c(3,3), col="salmon", lwd=3)
points(rahm_50, c(3,3), pch="|", col="salmon", cex=2)
points(a2_p$sle[35]/100, 3, pch="|", col="salmon3", cex=2) # NOTE 35 equals the year 2050
put.fig.letter("a.",font=2)

legend("topright", c("Rahmstorf [2007] 
A2 estimate","Rahmstorf [2007]
range","Method A","Method B","Method C"),
       bty="n", lty=c(1,1,1,1,1), lwd=2, y.intersp=1.5, cex=0.9, pt.cex=0.9,
       col=c("salmon3","salmon","gray","paleturquoise3","slateblue"))

# cumulative density function for SLR values in 2050
#par(mgp=c(1.5,.5,0), mar=c(3.75,4,2.75,2)) 
plot(cdf2050, lwd=3, ylab="Cumulative density", col="slateblue",
     xlab="Projected sea-level 2050 [m]", main="") #bootstrap SLR values
lines(b.cdf2050, lwd=3, col="gray") #heteroskedastic SLR values
lines(h.cdf2050, lwd=3, col="paleturquoise3") #homoskedastic SLR values
put.fig.letter("b.",font=2)
# legend("bottomright", c("Rahmstorf [2007] 
# A2 estimate","Rahmstorf [2007]
# range","Method A","Method B","Method C"),
#        bty="n", lty=c(1,1,1,1,1), lwd=2, y.intersp=1.5, cex=0.9, pt.cex=0.9,
#        col=c("salmon3","salmon","gray","paleturquoise3","slateblue"))

# survival function for SLR values in 2050
#par(mgp=c(1.5,.5,0), mar=c(5,4,2.5,2))
plot(survival2050$sf.num, survival2050$sf, log="y", type="l", xlim=c(0,0.8),
     col="slateblue", ylab="Survival function [1-cumulative frequency]", lwd=3,
     xlab="Projected sea-level 2050 [m]",main="",ylim=c(1e-04,1), yaxt="n") 
lines(b.survival2050$sf.num, b.survival2050$sf, lwd=3, col="gray") #heteroskedastic SLR values
lines(h.survival2050$sf.num, h.survival2050$sf, lwd=3, col="paleturquoise2") #homoskedastic SLR values
abline(h=c(0.01, 0.0001), lty=2) # add in the 1 in 100 and 1 in 10,000 lines
abline(v=a2_p$sle[35]/100, col="salmon3", lwd=3) # add in the Rahmstorf (2007) best A2 estimate
axis(2, at=10^(-4:0), label=parse(text=paste("10^", -4:0, sep="")))
text(0.10, 0.02, "1:100 level")
text(0.12, 0.0003, "1:10,000 
level")
put.fig.letter("c.",font=2)
dev.off()

########################################### SUPPLEMENTARY FIGURES #######################################
#---------------------------------- Supplementary Figure 1 ----------------------------------------------
# sup1: Tide gauge observations with method best fits and errors

# pdf(file="nRuckertetal_sup1_revise.pdf", family="Helvetica", width=6.7, height=2.7, pointsize=9)
# par(mfrow=c(1,2), mgp=c(1.5,.5,0), mar=c(4, 4, 3, 1)) # set figure dimensions
setEPS()
postscript(file="SuppFigures/sfigure1a_1b.eps", horizontal = FALSE, onefile = FALSE, paper = "special", family="Helvetica", width=6.7, height=2.7, pointsize=9)
par(mfrow=c(1,2), mgp=c(1.5,.5,0), mar=c(3.5,4,1,1))

plot(year,slr/100, pch=20, cex=0.5, xlab = "Year", ylab="Sea-level Anomaly [m]") # SLR observations
lines(year, err_neg/100, col="seashell4", lwd=2) # Negative observational errors
lines(year, err_pos/100, col="seashell4", lwd=2) # Positive observational errors
lines(year, hindcast$sle/100, col="salmon", lwd=2) # Best fit found in Rahmstorf (2007)
lines(year, b.new.slr$sle/100, col="gray", lwd=2)  # Bootstrap best fit
lines(year, new.est$sle/100, col="slateblue", lwd=2) # Heteroskedastic best fit
lines(year, h.new.est$sle/100, col="paleturquoise2", lwd=2) # Homoskedastic best fit
axis(side=4, labels=FALSE)
put.fig.letter("a.",font=2)

legend("topleft", c("Rahmstorf [2007] A2 estimate","Method A median fit","Method B median fit","Method C median fit",
                    "Observation error", "Observations"),bty="n", lwd=2, pch=c(NA,NA,NA,NA,NA,20), cex=0.85, 
       lty=c(1,1,1,1,1,NA), col=c("salmon","gray","paleturquoise2","slateblue","seashell4","black")) 

#Projection to 2050
#par(mgp=c(1.5,.5,0), mar=c(4, 3, 3, 2)) # set figure dimensions
plot(year, slr/100, type="l", xlab="Year",ylab="Sea-level anomalies [m]",
     col="gray",lwd=1,ylim=c(0,0.3), xlim=c(2000,2050))

lines(scenario_time, a2_p$sle/100, col="salmon", lwd=2) # Best fit found in Rahmstorf (2007)
lines(years.mod, b.new.rcp85$sle/100, col="gray", lwd=2)  # Bootstrap best fit
lines(years.mod, proj.est$sle/100, col="slateblue", lwd=2) # Heteroskedastic best fit
lines(years.mod, h.proj.est$sle/100, col="paleturquoise2", lwd=2) # Homoskedastic best fit

axis(side=4, labels=FALSE)
put.fig.letter("b.",font=2)

dev.off()

#---------------------------------- Supplementary Figure 2 ----------------------------------------------
# sup2: Trend of observation errors (heteroskedastic errors) from the sea-level record. The right panel plots  
#         the autocorrelation coefficient showing its significance in the sea-level record.
# pdf(file="nRuckertetal_sup2.pdf", family="Helvetica", width=6.7, height=2.7, pointsize=9)
# par(mfrow=c(1,2), mgp=c(1.5,.5,0), mar=c(4, 4, 3, 1)) # set figure dimensions

setEPS()
postscript(file="SuppFigures/sfigure2a_2b.eps", horizontal = FALSE, onefile = FALSE, paper = "special", family="Helvetica", width=6.7, height=2.7, pointsize=9)
par(mfrow=c(1,2), mgp=c(1.5,.5,0), mar=c(3.5,4,1,1))
plot(year, err.obs*10, typ="l", lwd=3,xlab="Years", ylab="Observation error [mm]", 
     col="black", ylim=c(0,max(err.obs*10)))
put.fig.letter("a.",font=2)

#par(mgp=c(1.5,.5,0), mar=c(4, 3, 3, 2)) # set figure dimensions
ac=acf(res, lag.max=5, plot=TRUE, main="",
       ylab="Residual autocorrelation coeffient",
       xlab = "Time lag")
put.fig.letter("b.",font=2)
dev.off()

#---------------------------------- Supplementary Figure 3 ----------------------------------------------
# sup3: 90% Credible interval of sea-level rise projected out to 2100. The bottom panel plots the
#         99% Credible intervalf sea-level estimates projected to 2100

#pdf(file="nRuckertetal_sup3.pdf", family="Helvetica",  width=3.5, height=5.4, pointsize=9)
#par(mfrow=c(2,1), mgp=c(1.5,.5,0), mar=c(3.5,4,4,2))
setEPS()
postscript(file="SuppFigures/sfigure3a_3b.eps", horizontal = FALSE, onefile = FALSE, paper = "special", family="Helvetica", width=3.5, height=5.4, pointsize=9)
par(mfrow=c(2,1), mgp=c(1.5,.5,0), mar=c(3.5,4,1,1))

# 90% Confidence intervals out to 2100.
plot(year, slr/100, type="l", xlab="Year",ylab="Sea-level anomalies [m]",
     col="gray",lwd=1,ylim=c(-.15,1.6), xlim=c(1950,2095))
# MCMC Heteroskedastic:
polygon(het_y_90, het_x_90, col="slateblue", border=NA)
# MCMC Homoskedastic:
polygon(hom_y_90, hom_x_90, col="paleturquoise2", border=NA)
# Bootstrap:
polygon(boot_y_90, boot_x_90, col="gray", border=NA) #90% polygon
lines(scenario_time, max_p$sle/100, lty=2, lwd=2) # Add in the max and min SLR estimated
lines(scenario_time, min_p$sle/100, lty=2, lwd=2) #  in Rahmstorf (2007)
axis(side=4, labels=FALSE)
legend("topleft", c("Rahmstorf [2007] range","90% CI Method A",
                    "90% CI Method B","90% CI Method C"),
       lty=c(2,NA,NA,NA),bty="n", lwd=2, pch=c(NA,15,15,15),
       col=c("black","gray","paleturquoise2","slateblue"))
put.fig.letter("a.",font=2)

# The entire range projected to 2100
#par(mgp=c(1.5,.5,0), mar=c(5,4,2.5,2))
plot(year, slr/100, type="l", xlab="Year",
     ylab="Sea-level anomalies [m]",col="gray",lwd=1,
     ylim=c(-.15,1.6), xlim=c(1950,2095))
# MCMC heteroskedastic#
het_x_99=c(het_point5, rev(het_995)); het_y_99=c(years.mod, rev(years.mod))
polygon(het_y_99, het_x_99, col="slateblue", border=NA)
# MCMC Homoskedastic#
hom_x_99=c(hom_point5, rev(hom_995)); hom_y_99=c(years.mod, rev(years.mod))
polygon(hom_y_99, hom_x_99, col="paleturquoise2", border=NA)
# Method A#
boot_x_99=c(boot_point5, rev(boot_995)); boot_y_99=c(years.mod, rev(years.mod)) # set up the polygon
polygon(boot_y_99, boot_x_99, col="gray", border=NA) # plot 99% polygon

lines(scenario_time, max_p$sle/100, lty=2, lwd=2) # Add in the max and min SLR estimated
lines(scenario_time, min_p$sle/100, lty=2, lwd=2) #  in Rahmstorf (2007)
axis(side=4, labels=FALSE)
put.fig.letter("b.",font=2)
legend("topleft", c("99% CI Method A","99% CI Method B","99% CI Method C"), 
       pch=c(15,15,15),bty="n",col=c("gray","paleturquoise2","slateblue"))
dev.off()

#---------------------------------- Supplementary Figure 5 ----------------------------------------------
# sup5: Comparison of SLR projections

studies.num = 1:9
width = 0.25  # Width of bars

# Method A, Method B, Method C, Kopp 2014, IPCC AR4 (A2),IPCC AR5 (RCP85), Rahmstorf 2012 (default with CW06), Grinsted 2009, VR09 with JE08 #0.98 to 1.55 was the +/- 1 standard deviation for VR
low_5_dist = c(0.68, 0.50, 0.32, 0.5, 0.23, 0.52,1.125, 0.6, 0.9)
high_95_dist = c(0.98, 0.96, 1.07, 1.2, 0.51, 0.98,1.375, 0.85, 1.65)

studies.colors <- c("gray", "paleturquoise2", "slateblue", "black", "black", "black", "black", "black", "black")

setEPS()
postscript(file="SuppFigures/sfigure5.eps", horizontal = FALSE, onefile = FALSE, paper = "special", family="Helvetica", width=3.5, height=2.7, pointsize=9)
par(mfrow=c(1,1), mgp=c(1.5,.5,0),  mar=c(3.5,4,1,1)) # set figure dimensions
#par(mfrow=c(1,1), mgp=c(1.5,.5,0),mar=c(4, 4, 3, 2))
plot(c(0.5,9), c(0,2), pch = 20, col = "white", 
     xaxt="n", ylab="Sea-level anomaly in 2100 [m]", xlab="")

for(i in 1:length(studies.num)) {
  polygon(x = c((studies.num[i]), (studies.num[i]),  (studies.num[i] - width), (studies.num[i]-width)),
          y = c(low_5_dist[i], high_95_dist[i], high_95_dist[i], low_5_dist[i]), 
          border = NA,
          col = studies.colors[i])
}

ticks = studies.num
axis(side=1, at=ticks, labels=expression("A","B","C","K14","AR4","AR5", "R12", "G10", "VR09"))

dev.off() 

#---------------------------------- Supplementary Figure 7 ----------------------------------------------
# sup7: A probability density function of sea-level rise in the year 2100.
#       Followed by the cumulative desity function and the surival function of sea-level rise
#       estimated in each method for the year 2100

#pdf(file="nRuckertetal_sup6.pdf", family="Helvetica", width=3.5, height=8.1, pointsize=13)
#par(mfrow=c(3,1), mgp=c(1.5,.5,0),mar=c(2.5,4,4,2))

setEPS()
postscript(file="SuppFigures/sfigure9a_9c.eps", horizontal = FALSE, onefile = FALSE, paper = "special", family="Helvetica", width=3.5, height=8.1, pointsize=13)
par(mfrow=c(3,1), mgp=c(1.5,.5,0), mar=c(3.5,4,1,1))

# Probability desity function for 2100
plot(b.pdf2100, main="",lwd=3, col="gray", xlab="Projected sea-level 2100 [m]", 
ylab="Probability Density", xlim=c(-0.5,3), yaxt="n") # homoskedastic SLR values
lines(h.pdf2100, col="paleturquoise2", lwd=3) # bootstrap SLR values
lines(pdf2100, col="slateblue", lwd=3) # heteroskedastic SLR values
rahm_50 = c(min_p$sle[45]/100,max_p$sle[45]/100)
lines(rahm_50, c(1,1), col="salmon", lwd=3) # set up the range generated in rahmstorf (2007) for 2050
points(a2_p$sle[45]/100, 1, pch="|", col="salmon3", cex=2) # NOTE: 45 equals the year 2100
points(rahm_50, c(1,1), pch="|", col="salmon", cex=2)
put.fig.letter("a.",font=2)
legend("topright", c("Rahmstorf [2007] 
A2 estimate","Rahmstorf [2007]
range","Method A","Method B","Method C"),
       bty="n", lty=c(1,1,1,1,1), lwd=2, y.intersp=1.5, cex=0.9, pt.cex=0.9,
       col=c("salmon3","salmon","gray","paleturquoise3","slateblue"))

# Cumulative desity function for 2100
par(mgp=c(1.5,.5,0), mar=c(3.75,4,2.75,2)) 
plot(cdf2100, col="slateblue", ylab="Cumulative density", main="",
xlab="Projected sea-level 2100 [m]", lwd=3) # bootstrap SLR values
lines(b.cdf2100, col="gray", lwd=3) # heteroskedastic SLR values
lines(h.cdf2100, col="paleturquoise2", lwd=3) # homoskedastic SLR values
put.fig.letter("b.",font=2)

# Survival function for 2100
par(mgp=c(1.5,.5,0), mar=c(5,4,2.5,2))
plot(survival2100$sf.num, survival2100$sf, log="y", type="l", xlim=c(0,3),
     col="slateblue", ylab="Survival function [1-cumulative frequency]", lwd=3,
     xlab="Projected sea-level 2100 [m]",main="",ylim=c(1e-04,1), yaxt="n") 
lines(b.survival2100$sf.num, b.survival2100$sf, lwd=3, col="gray") #heteroskedastic SLR values
lines(h.survival2100$sf.num, h.survival2100$sf, lwd=3, col="paleturquoise2") #homoskedastic SLR values
abline(h=c(0.01, 0.0001), lty=2) # add in the 1 in 100 and 1 in 10,000 lines
abline(v=a2_p$sle[45]/100, col="salmon3", lwd=3) # add in the Rahmstorf (2007) best A2 estimate
axis(2, at=10^(-4:0), label=parse(text=paste("10^", -4:0, sep="")))
text(0.35, 0.015, "1:100 level")
text(0.30, 0.0003, "1:10,000 
level")
put.fig.letter("c.",font=2)
dev.off()

####################################### END ################################################################
