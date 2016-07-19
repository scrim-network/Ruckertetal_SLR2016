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
#
###################################################################################################

load("Workspace/mega_R_methods_workspace.RData") # Load in the saved workspace from file Mega_Rahmstorf.R

source("Scripts/put_fig_letter.R")
######################################## MAIN FIGURES #############################################

# Calculate the 90% confidence interval for each method and upper tail estimates
het_1in30 <-
het_1in100 <-
het_1in1000 <-
hom_1in30 <-
hom_1in100 <-
hom_1in1000 <-
boot_1in30 <-
boot_1in100 <-
boot_1in1000 <-

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
    het_1in30[i] <- quantile(SLR.projections.heter[,i]/100, (1-(1/30)))
    het_1in100[i] <- quantile(SLR.projections.heter[,i]/100, (1-(1/100)))
    het_1in1000[i] <- quantile(SLR.projections.heter[,i]/100, (1-(1/1000)))
    hom_1in30[i] <- quantile(SLR.projections.homo[,i]/100, (1-(1/30)))
    hom_1in100[i] <- quantile(SLR.projections.homo[,i]/100, (1-(1/100)))
    hom_1in1000[i] <- quantile(SLR.projections.homo[,i]/100, (1-(1/1000)))
    boot_1in30[i] <- quantile(SLR.projections.boot[,i]/100, (1-(1/30)))
    boot_1in100[i] <- quantile(SLR.projections.boot[,i]/100, (1-(1/100)))
    boot_1in1000[i] <- quantile(SLR.projections.boot[,i]/100, (1-(1/1000)))
    
    het_5[i] = quantile(SLR.projections.heter[,i]/100,0.05) #heteroskedastic SLR values 90%
    het_95[i] = quantile(SLR.projections.heter[,i]/100,0.95)
    hom_5[i] = quantile(SLR.projections.homo[,i]/100,0.05) #homostedastic SLR values 90%
    hom_95[i] = quantile(SLR.projections.homo[,i]/100,0.95)
    boot_5[i] = quantile(SLR.projections.boot[,i]/100,0.05) # Bootstrap SLR values 90%
    boot_95[i] = quantile(SLR.projections.boot[,i]/100,0.95)
    het_point5[i] = quantile(SLR.projections.heter[,i]/100,0.005) #heteroskedastic SLR values 99%
    het_995[i] = quantile(SLR.projections.heter[,i]/100,0.995)
    hom_point5[i] = quantile(SLR.projections.homo[,i]/100,0.005) #homostedastic SLR values 99%
    hom_995[i] = quantile(SLR.projections.homo[,i]/100,0.995)
    boot_point5[i] = quantile(SLR.projections.boot[,i]/100,0.005) # Bootstrap SLR values 99%
    boot_995[i] = quantile(SLR.projections.boot[,i]/100,0.995)
}

# Print the sea-level anomaly with 3.3, 1, and 0.1% probability in the year 2050 and 2100:
# 2050         ; 2100
het_1in30[171] ; het_1in30[221]
het_1in100[171] ; het_1in100[221]
het_1in1000[171] ; het_1in1000[221]

hom_1in30[171] ; hom_1in30[221]
hom_1in100[171] ; hom_1in100[221]
hom_1in1000[171] ; hom_1in1000[221]

boot_1in30[171] ; boot_1in30[221]
boot_1in100[171] ; boot_1in100[221]
boot_1in1000[171] ; boot_1in1000[221]

# Print the percent increase: in 2050
# 1/30 (3.3%)
(het_1in30[171] - boot_1in30[171])/ boot_1in30[171])
# 1/100 (1%)
(het_1in100[171] - boot_1in100[171])/ boot_1in100[171])
# 1/1,000 (0.1%)
(het_1in1000[171] - boot_1in1000[171])/ boot_1in1000[171])

# Print the percent increase: in 2100
# 1/30 (3.3%)
(het_1in30[221] - boot_1in30[221])/ boot_1in30[221])
# 1/100 (1%)
(het_1in100[221] - boot_1in100[221])/ boot_1in100[221])
# 1/1,000 (0.1%)
(het_1in1000[221] - boot_1in1000[221])/ boot_1in1000[221])

# Set plotting dimensions
mm_TO_inches = function(mm){
    mm * 0.039370
}

single_column = mm_TO_inches(84)
med_column = mm_TO_inches(129)
double_column = mm_TO_inches(174)
maximum_width = mm_TO_inches(234)
column_height=2.7

library(RColorBrewer)
test.colors = brewer.pal(9, "YlGnBu")
test.colors = c(test.colors[4], test.colors[6], test.colors[8])

# set up the polygons
het_x_90=c(het_5, rev(het_95)); het_y_90=c(years.mod, rev(years.mod))
hom_x_90=c(hom_5, rev(hom_95)); hom_y_90=c(years.mod, rev(years.mod))
boot_x_90=c(boot_5, rev(boot_95)); boot_y_90=c(years.mod, rev(years.mod))

#---------------------------------- Figure 1 ----------------------------------------------
# Residuals (a); trend of observation errors (b); residual auto-correlation coefficient (c)

setEPS()
# postscript(file="SuppFigures/sfigure2a_2b.eps", horizontal = FALSE, onefile = FALSE, paper = "special", family="Helvetica", width=6.7, height=2.7, pointsize=9)
postscript(file="Figures/Fig1.eps", horizontal = FALSE, onefile = FALSE, paper = "special", family="Helvetica",
width=double_column, height=column_height*2, pointsize=11)

layout(matrix(c(1,1,2,3), 2, 2, byrow = TRUE))
par(mgp=c(1.5,.5,0), mar=c(3.5,4,1,1))

# Residuals (a)
plot(year, res*10, type="l", xlab="Year", ylab="Sea-level residuals [mm]
(observations - best fit)")
points(year, res*10, pch=20)
abline(h=0, lty=3, lwd=2)
put.fig.letter("a.",font=2)

# Trend of observation errors
plot(year, err.obs*10, typ="l", xlab="Years", ylab="Observation error [mm]",
col="black", ylim=c(0, max(err.obs*10)))
points(year, err.obs*10, pch=20)
put.fig.letter("b.",font=2)

# Residual auto-correlation coefficient
ac=acf(res, lag.max=20, plot=TRUE, main="", lwd=4,
ylab="Residual autocorrelation coeffient",
xlab = "Time lag")
put.fig.letter("c.",font=2)
dev.off()

#---------------------------------- Figure 2 ----------------------------------------------
# 2) 90% credible interval hindcast (a); relibability diagram and subset (b); 90% credible interval projection to 2050 (c);
#    90% credible interval projections to 2100 (c)
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

setEPS()
# postscript(file="Figures/figure1a_1b.eps", horizontal = FALSE, onefile = FALSE, paper = "special", family="Helvetica", width=6.7, height=2.7, pointsize=9)
postscript(file="Figures/Fig2.eps", horizontal = FALSE, onefile = FALSE, paper = "special", family="Helvetica",
width=double_column, height=column_height*2, pointsize=11)
par(mfrow=c(2,2), mgp=c(1.5,.5,0), mar=c(3.5,4,1,1))

# 90% credible interval hindcast
plot(year, slr/100, type="l",xlab="Year",ylab="Sea-level anomalies [m]",
ylim=c(-0.16,0.07), xlim=c(1884,1998))
polygon(het_y_90, het_x_90, col=test.colors[3], border=NA)
polygon(hom_y_90, hom_x_90, col=test.colors[2], border=NA)
polygon(boot_y_90, boot_x_90, col=test.colors[1], border=NA) # plot 90% polygon
points(year, slr/100, pch=20, cex=0.5, col="black") # add in the observations
arrows(year, err_pos/100, year, err_neg/100, length=0)  # add in the observational errors
axis(side=4, labels=FALSE)
put.fig.letter("a.",font=2)

legend("topleft", c('Bootstrap','Bayesian (homoskedastic)','Bayesian (heteroskedastic)', "Perfect reliability", "Observations"), cex=0.85,
pch=c(15,15,15,NA,20), lty=c(NA,NA,NA,2,NA), bty="n", col=c(test.colors[1:3],"black","black"))
legend("topleft", c('','', "","","", "Observation error"),
pch=c(NA,NA,NA,NA,NA,"l"), bty="n", col="black", cex=0.85)

# Relibability diagram
plot(percent, percent, typ="l", lty=2, ylab="Observed relative frequency [%]",
xlab="Forecast probability [%] (Credible interval)", lwd=2)
points(percent, boots*100, pch=10, col=test.colors[1], lwd=2)
points(percent, homos*100, pch=8, col=test.colors[2], lwd=2)
points(percent, heters*100, pch=4, col=test.colors[3], lwd=2)
text(25,40, "underconfident", srt=35)
text(35,25, "overconfident", srt=35)
axis(side=4, labels=FALSE)
put.fig.letter("b.",font=2)

# 90% credible interval projections to 2050
plot(year, slr/100, type="l", xlab="Year",ylab="Sea-level anomalies [m]",
col="gray",lwd=1,ylim=c(0,0.4), xlim=c(2000,2050))
polygon(het_y_90, het_x_90, col=test.colors[3], border=NA)
polygon(hom_y_90, hom_x_90, col=test.colors[2], border=NA)
polygon(boot_y_90, boot_x_90, col=test.colors[1], border=NA)
axis(side=4, labels=FALSE)
put.fig.letter("c.",font=2)

# 90% credible interval projections to 2100
plot(year, slr/100, type="l", xlab="Year",ylab="Sea-level anomalies [m]",
col="gray",lwd=1,ylim=c(0,1.1), xlim=c(2050,2100))
polygon(het_y_90, het_x_90, col=test.colors[3], border=NA)
polygon(hom_y_90, hom_x_90, col=test.colors[2], border=NA)
polygon(boot_y_90, boot_x_90, col=test.colors[1], border=NA)
axis(side=4, labels=FALSE)
put.fig.letter("d.", font=2)

# Relibability diagram subplot (zoom in on >90%)
par(fig = c(.63,.77,.82,.96), mar=c(0,0,0,0), new=TRUE)
plot(percent, percent, type="l", lty=2, ylab="Observed relative frequency [%]",
xlab="Forecast probability [%] (credible interval)", lwd=2, xlim=c(90,100), ylim=c(90,100))
points(percent, boots*100, pch=10, col=test.colors[1], lwd=2)
points(percent, homos*100, pch=8, col=test.colors[2], lwd=2)
points(percent, heters*100, pch=4, col=test.colors[3], lwd=2)
dev.off()

#---------------------------------- Figure 3 ----------------------------------------------
# 3) Probability density functions of each of the parameters estimated

setEPS()
# postscript(file="Figures/figure3a_3c.eps", horizontal = FALSE, onefile = FALSE, paper = "special", family="Helvetica", width=3.5, height=8.1, pointsize=13)
postscript(file="Figures/Fig3.eps", horizontal = FALSE, onefile = FALSE, paper = "special", family="Helvetica",
width=double_column, height=column_height*2, pointsize=11)
par(mfrow=c(2,2), mgp=c(1.5,.5,0), mar=c(3.5,4,1,1))

# Sea-level sensitivity (alpha)
plot(boot.pdfa, main="", lwd=3, col=test.colors[1],xlab=expression(paste(alpha, " [cm/yr/", degree,"C] (sensitivity of SLR)")),# "Sensitivity of SLR (a) [cm/yr/C]",
ylab="Probability Density", xlim=c(-0.01,2), yaxt="n")
lines(homo.pdfa, lwd=3, col=test.colors[2])
lines(heter.pdfa, lwd=3, col=test.colors[3])
# Add in a lines to create the prior distribution
lines(c(0,0),c(0,1), lwd=3, lty=2)
lines(c(2,2),c(0,1), lwd=3, lty=2)
lines(c(0,2),c(1,1), lwd=3, lty=2)
put.fig.letter("a.",font=2)

legend("topright", c("Bootstap","Bayesian
(homoskedastic)","Bayesian
(heteroskedastic)", "Prior"), bty="n", lty=c(1,1,1,2), lwd=2, y.intersp=1.5, cex=0.9,
pt.cex=0.9, col=c(test.colors[1:3],"black"))

#Base Temperature (T0)
plot(boot.pdfTo, lwd=3, col=test.colors[1], main="",xlab=expression(paste(T[0], " [",degree,"C] (equilibrium temperature)")),
ylab="Probability Density", xlim=c(-3.3,2.1), yaxt="n")
lines(homo.pdfTo, lwd=3, col=test.colors[2])
lines(heter.pdfTo, lwd=3, col=test.colors[3])
lines(c(-3,-3),c(0,0.5), lwd=3, lty=2)
lines(c(2,2),c(0,0.5), lwd=3, lty=2)
lines(c(-3,2),c(0.5,0.5), lwd=3, lty=2)
put.fig.letter("b.",font=2)

#Initial value of sea-level in 1880
plot(boot.pdfinitialval, lwd=3, col=test.colors[1], main="",xlab=expression(paste(H[0], " [cm] (initial value of sea-level)")),
ylab="Probability Density", xlim=c(-17.3,-12), yaxt="n")
lines(homo.pdfinitialval, lwd=3, col=test.colors[2])
lines(heter.pdfinitialval, lwd=3, col=test.colors[3])
lines(c(-16.9,-16.9),c(0,.25), lwd=3, lty=2)
lines(c(-12.4,-12.4),c(0,.25), lwd=3, lty=2)
lines(c(-16.9,-12.4),c(.25,.25), lwd=3, lty=2)
put.fig.letter("c.",font=2)

plot(homo.pdfsigma, lwd=3, col=test.colors[2], main="",xlab=expression(paste(rho, " (lag-1 autoregression coefficient)")),
ylab="Probability Density", xlim=c(-0.99, 0.99), yaxt="n")
abline(v=rho[2], lwd=3, col=test.colors[1])
lines(heter.pdfrho, lwd=3, col=test.colors[3])
lines(c(-0.99,-0.99),c(0,.5), lwd=3, lty=2)
lines(c(0.99,0.99),c(0,.5), lwd=3, lty=2)
lines(c(-0.99,0.99),c(.5,.5), lwd=3, lty=2)
put.fig.letter("d.",font=2)

dev.off()

#---------------------------------- Figure 4 ----------------------------------------------
# 4) A probability density function (a), CDF (b), and SF (c) of sea-level rise in the year 2050 & 2100.

setEPS()
# postscript(file="Figures/figure4a_4c.eps", horizontal = FALSE, onefile = FALSE, paper = "special", family="Helvetica", width=3.5, height=8.1, pointsize=13)
postscript(file="Figures/Fig4.eps", horizontal = FALSE, onefile = FALSE, paper = "special", family="Helvetica",
width=double_column, height=column_height*3, pointsize=13)
par(mfrow=c(3,2), mgp=c(1.5,.5,0), mar=c(3.5,4,1,1))

# probability density function for SLR values in 2050
plot(pdf2050_boot, main="", col=test.colors[1], lwd=3, ylab="Probabilty density",
xlab="", xlim=c(-0.1,0.8), yaxt="n")
lines(pdf2050_homo, lwd=3, col=test.colors[2])
lines(pdf2050_heter, lwd=3, col=test.colors[3])
put.fig.letter("a.",font=2)
legend("topright", c("Bootstap","Bayesian
(homoskedastic)","Bayesian
(heteroskedastic)"), bty="n", lty=c(1,1,1), lwd=2, y.intersp=1.5, cex=0.9,
pt.cex=0.9, col=test.colors[1:3])

# Probability desity function for 2100
plot(pdf2100_boot, main="",lwd=3, col=test.colors[1], xlab="",
ylab="", xlim=c(0,2.5), yaxt="n")
lines(pdf2100_homo, col=test.colors[2], lwd=3)
lines(pdf2100_heter, col=test.colors[3], lwd=3)
put.fig.letter("b.",font=2)

# cumulative density function for SLR values in 2050
plot(cdf2050_heter, lwd=3, ylab="Cumulative density", col=test.colors[3],
xlab="", main="", xlim=c(-0.1,0.8))
lines(cdf2050_boot, lwd=3, col=test.colors[1])
lines(cdf2050_homo, lwd=3, col=test.colors[2])
put.fig.letter("c.",font=2)

# Cumulative desity function for 2100
plot(cdf2100_heter, col=test.colors[3], ylab="", main="",
xlab="", lwd=3, xlim=c(0,2.5))
lines(cdf2100_boot, col=test.colors[1], lwd=3)
lines(cdf2100_homo, col=test.colors[2], lwd=3)
put.fig.letter("d.",font=2)

# survival function for SLR values in 2050
plot(survival2050_heter$sf.num, survival2050_heter$sf, log="y", type="l", xlim=c(-0.1,0.8),
col=test.colors[3], ylab="Survival function
[1-cumulative probability]", lwd=3,
xlab="Projected sea-level 2050 [m]",main="",ylim=c(1e-04,1), yaxt="n")
lines(survival2050_boot$sf.num, survival2050_boot$sf, lwd=3, col=test.colors[1])
lines(survival2050_homo$sf.num, survival2050_homo$sf, lwd=3, col=test.colors[2])
abline(h=1/100, lty=2)
axis(2, at=10^(-4:0), label=parse(text=paste("10^", -4:0, sep="")))
text(0.10, 0.015, expression(paste("10"^"-2", sep="")))
put.fig.letter("e.",font=2)

# Survival function for 2100
plot(survival2100_heter$sf.num, survival2100_heter$sf, log="y", type="l", xlim=c(0,2.5),
col=test.colors[3], ylab="", lwd=3,
xlab="Projected sea-level 2100 [m]",main="",ylim=c(1e-04,1), yaxt="n")
lines(survival2100_boot$sf.num, survival2100_boot$sf, lwd=3, col=test.colors[1])
lines(survival2050_homo$sf.num, survival2050_homo$sf, lwd=3, col=test.colors[2])
abline(h=1/100, lty=2)
axis(2, at=10^(-4:0), label=parse(text=paste("10^", -4:0, sep="")))
put.fig.letter("f.",font=2)

dev.off()

########################################### SUPPLEMENTARY FIGURES #######################################
#---------------------------------- Supplementary Figure 1 ----------------------------------------------
# sup1: Tide gauge observations with method best fits and errors

setEPS()
# postscript(file="SuppFigures/sfigure1a_1b.eps", horizontal = FALSE, onefile = FALSE, paper = "special", family="Helvetica", width=6.7, height=2.7, pointsize=9)
postscript(file="SuppFigures/SFig1.eps", horizontal = FALSE, onefile = FALSE, paper = "special", family="Helvetica",
width=single_column, height=column_height, pointsize=9)
par(mfrow=c(1,1), mgp=c(1.5,.5,0), mar=c(3.5,4,1,1))

plot(year,slr/100, pch=20, cex=0.5, xlab = "Year", ylab="Sea-level Anomaly [m]")
lines(year, err_neg/100, col="seashell4", lwd=2)
lines(year, err_pos/100, col="seashell4", lwd=2)
lines(year, boot.med.hindcast$sle/100, col=test.colors[1], lwd=2)
lines(year, heter.med.hindcast$sle/100, col=test.colors[3], lwd=2)
lines(year, homo.med.hindcast$sle/100, col=test.colors[2], lwd=2)
axis(side=4, labels=FALSE)

legend("topleft", c("Bootstrap median fit","Bayesian (homoskedastic) median fit","Bayesian (heteroskedastic) median fit",
"Observation error", "Observations"),bty="n", lwd=2, pch=c(NA,NA,NA,NA,20), cex=0.8,
lty=c(1,1,1,1,NA), col=c(test.colors[1:3],"seashell4","black"))

#---------------------------------- Supplementary Figure 2 ----------------------------------------------
# sup2: Comparison of SLR projections
studies.num = 1:7
width = 0.25  # Width of bars

# Method A, Method B, Method C,  Kopp 2016, Kopp et al 2014, # Rahmstorf 2012 (default with CW06),
# VR09 with JE08 #0.98 to 1.55 was the +/- 1 standard deviation for VR

#low_5_dist = c(0.68, 0.50, 0.32, 0.52,0.55, 1.125, 0.9)
#high_95_dist = c(0.98, 0.96, 1.07,1.31,1.21,1.375, 1.65)

low_5_dist =   c(round(boot_5,2),  round(hom_5,2),  round(het_5,2),  0.52, 0.55, 1.125, 0.9)
high_95_dist = c(round(boot_95,2), round(hom_95,2), round(het_95,2), 1.31, 1.21, 1.375, 1.65)

studies.colors <- c(test.colors[1:3], "black", "black", "black", "black", "black", "black")

setEPS()
# postscript(file="SuppFigures/sfigure5.eps", horizontal = FALSE, onefile = FALSE, paper = "special", family="Helvetica", width=3.5, height=2.7, pointsize=9)
postscript(file="SuppFigures/SFig2.eps", horizontal = FALSE, onefile = FALSE, paper = "special", family="Helvetica",
width=single_column, height=column_height, pointsize=9)
par(mfrow=c(1,1), mgp=c(1.5,.5,0),  mar=c(3.5,4,1,1)) # set figure dimensions

plot(c(0.5,7), c(0,2), pch = 20, col = "white",
xaxt="n", ylab="Sea-level anomaly in 2100 [m]", xlab="")

for(i in 1:length(studies.num)) {
    polygon(x = c((studies.num[i]), (studies.num[i]),  (studies.num[i] - width), (studies.num[i]-width)),
    y = c(low_5_dist[i], high_95_dist[i], high_95_dist[i], low_5_dist[i]),
    border = NA,
    col = studies.colors[i])
}

ticks = studies.num
axis(side=1, at=ticks, labels=expression("Boot.","Hom.","Het.","K16","K14", "R12", "VR09"))

dev.off()

####################################### END ################################################################
