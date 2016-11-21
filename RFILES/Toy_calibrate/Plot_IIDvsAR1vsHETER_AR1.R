#################################################################################
#
#  -file = "Plot_IIDvsAR1vsHETER_AR1.R"   Code written February 2015
#  - Author: Kelsey Ruckert (klr324@psu.edu)
#
#  -The toy scripts analyze three hypotheses for failure in the surprise index
#       as described in Ruckert et al. (2016). The three hypotheses include number
#       of observations, AR(1) observations, and observations assuming heteroskedastic
#       errors. It also produces the last three supplementary figures shown in that
#       paper. The graphs will be saved as pdf files in the current working directory.
#       For further description and references, please read the paper and appendix.
#
# THIS CODE IS PROVIDED AS-IS WITH NO WARRANTY (NEITHER EXPLICIT
# NOT IMPLICIT).  I SHARE THIS CODE IN HOPES THAT IT IS USEFUL,
# BUT I AM NOT LIABLE FOR THE BEHAVIOR OF THIS CODE IN YOUR OWN
# APPLICATION.  YOU ARE FREE TO SHARE THIS CODE SO LONG AS THE
# AUTHOR(S) AND VERSION HISTORY REMAIN INTACT.
#
#   -NOTE: This program applies a Markov Chain Monte Carlo method with varying
#       assumptions about the observations including:
#
#       -varies number of observations
#       -preserves the iid structure of the observations
#       -creates random realizations of "natural variability"
#       -preserves the AR(1) structure of the observations
#       -Markov Chain Monte Carlo AR(1) Heteroskedastic errors
#       -Markov Chain Monte Carlo AR(1) Homoskedastic errors
#
#  -The spreadsheet: "testhypoth_surprise.csv" is a previously calculated surprise index for
#       all the following tests.
#
###################################################################################
rm(list =ls()) #Clear global environment
library(compiler)
enableJIT(3)

source("Scripts/put_fig_letter.R")
library(RColorBrewer)
test.colors = brewer.pal(3, "RdBu")
test.colors[2] = "gray"
circlesize = seq(from=1.2, to=2.5, length.out=3)

################################################################################################
#------------------------------- Plot the Results -----------------------------------#
# NOTE: For Ease the surprise index for each test has been previously calculated for ease:
# load("Workspace/testiidall.RData") #Load saved workspace
# load("Workspace/testiidall_hettest.RData") #Load saved workspace

load("Workspace/test_IID_data.RData")
load("Workspace/test_AR1_data.RData")
load("Workspace/test_het_data.RData")

# i.i.d observations were deleted in the process. Regenarte the i.i.d. observations:
IIDobs = true.obs$mod.obs + residuals

mm_TO_inches = function(mm){
  mm * 0.039370
}

single_column = mm_TO_inches(84)
double_column = mm_TO_inches(174)
maximum_width = mm_TO_inches(234)
column_height=2.7
#------------------------------- Sup. Figure 8 -----------------------------------#
# Compare the surprise index and the pdfs for iid, AR(1), and heteroskedastic erros:
#pdf(file="SuppFigures/nRuckertetal_sup9.pdf", family="Helvetica", height=5.4, width=6.7, pointsize=11)
setEPS()
# postscript(file="SuppFigures/sfigure8a_8d.eps", horizontal = FALSE, onefile = FALSE, paper = "special", family="Helvetica", width=6.7, height=5.4, pointsize=11)
postscript(file="SuppFigures/sFig8a_8c_fixed.eps", horizontal = FALSE, onefile = FALSE, paper = "special", family="Helvetica",
           width=double_column, height=column_height*2, pointsize=11)
par(mfrow=c(2,2), mgp=c(1.5,.5,0),  mar=c(3.5,4,1,1)) # set figure dimensions
# 1) PDF
#par(mfrow=c(2,2),mgp=c(1.5,.5,0), mar=c(3.5,4,4,1))
plot(twoh_a, xlab="a", ylab="Probability density", main="",
     lwd=3, col=test.colors[1], yaxt="n", xlim=c(0,2)) #, xlim=c(-10,10))
lines(homtwoh_a, col=test.colors[2],lwd=3)
lines(hettwoh_a, col=test.colors[3],lwd=3)
abline(v=true_par[1], lty=3, lwd=2)
# lines(c(-10,-10),c(0,0.18), lwd=2, lty=3)
# lines(c(10,10),c(0,0.18), lwd=2, lty=3)
# lines(c(-10,10),c(0.18,0.18), lwd=2, lty=3)
put.fig.letter("a.",font=2)
legend("topright", c("True parameter", "1:1 line", "IID","Homoskedastic","Heteroskedastic"), cex=0.75, pt.cex=c(NA,NA,circlesize[1:3]), y.intersp = 1.5,
       lty=c(3,2,1,1,1),pch=c(NA,NA,20,20,20), lwd=2,col=c("black", "black", test.colors[1:3]))

# 2) Surprise index
# surprise = read.csv("ToyScripts/testhypoth_surprise.csv")
surprise = read.csv("ToyScripts/test_surprise_bias.csv")
percent = surprise[1:16,1]
homtwohund = surprise[1:16,3]
hettwohund = surprise[1:16,5]
twohund = surprise[1:16,7]
#par(mgp=c(1.5,.5,0), mar=c(3.5, 3, 4, 2))
plot(percent, percent, typ="l", lty=2, ylab="Observed relative frequency [%]", 
     xlab="Forecast probability [%] (credible interval)", lwd=2)
# points(percent, twohund*100,  pch=20, col=test.colors[1], cex=circlesize[1])
# points(percent, homtwohund*100,  pch=20, col=test.colors[2], cex=circlesize[2])
# points(percent, hettwohund*100,  pch=19, col=test.colors[3], cex=circlesize[3])
#----
points(percent, hettwohund*100,  pch=19, col=test.colors[3], cex=circlesize[3])
points(percent, homtwohund*100,  pch=20, col=test.colors[2], cex=circlesize[2])
points(percent, twohund*100,  pch=20, col=test.colors[1], cex=circlesize[1])
put.fig.letter("c.",font=2)
lines(c(88,102),c(102,102))
lines(c(88,88),c(88,102))
lines(c(102,102),c(88,102))
lines(c(88,102),c(88,88))
text(100,92, font=2, "d.")
text(50,65, "underconfident", srt=35)
text(55,35, "overconfident", srt=35)

#par(mgp=c(1.5,.5,0), mar=c(5,4,2.5,1))
plot(twoh_b, xlab="b", ylab="Probability density", main="",lwd=3,
     col=test.colors[1], yaxt="n", xlim=c(0.19,0.21)) #xlim=c(-10,10))
lines(homtwoh_b, col=test.colors[2],lwd=3)
lines(hettwoh_b, col=test.colors[3],lwd=3)
abline(v=true_par[2], lty=3, lwd=2)
# lines(c(-10,-10),c(0,50), lwd=2, lty=3)
# lines(c(10,10),c(0,50), lwd=2, lty=3)
# lines(c(-10,10),c(50,50), lwd=2, lty=3)
put.fig.letter("b.",font=2)

#par(mgp=c(1.5,.5,0), mar=c(5, 3, 2.5, 2))
plot(percent, percent, type="l", lty=2, ylab="Observed relative frequency [%]", 
     xlab="Forecast probability [%] (credible interval)", lwd=2, xlim=c(90,100), ylim=c(90,100))
# points(percent, twohund*100,  pch=20, col=test.colors[1], cex=circlesize[1])
# points(percent, homtwohund*100,  pch=20, col=test.colors[2], cex=circlesize[2])
# points(percent, hettwohund*100,  pch=19, col=test.colors[3], cex=circlesize[3])
#----
points(percent, hettwohund*100,  pch=19, col=test.colors[3], cex=circlesize[3])
points(percent, homtwohund*100,  pch=20, col=test.colors[2], cex=circlesize[2])
points(percent, twohund*100,  pch=20, col=test.colors[1], cex=circlesize[1])
put.fig.letter("d.",font=2)
text(92,92.75, "underconfident", srt=35)
text(93,92, "overconfident", srt=35)
# legend("bottomright", c("IID","Homoskedastic","Heteroskedastic","1:1 line"),
#        pch=c(10,8,4,NA), lty=c(NA,NA,NA,2), lwd=2,cex=0.75, 
#        col=c("gold","gray","cornflowerblue","black"))
dev.off()

######################################## END ###################################################
