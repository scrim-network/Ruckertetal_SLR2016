#################################################################################
#
#  -file = "Toy_TestMCMC_plots.R"   Code written February 2015
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
load("Workspace/testiidall.RData") #Load saved workspace
load("Workspace/testiidall_hettest.RData") #Load saved workspace

# i.i.d observations were deleted in the process. Regenarte the i.i.d. observations:
IIDobs = true.obs$mod.obs + residuals

#------------------------------- Sup. Figure 6 -----------------------------------#
# Compare the 90% Confidence intervals for iid, AR(1), and heteroskedastic error tests
# pdf(file="SuppFigures/nRuckertetal_sup7.pdf", family="Helvetica", width=3.5, height=2.7, pointsize=9)
# par(mfrow=c(1,1), mgp=c(1.5,.5,0),mar=c(4, 4, 3, 2))
setEPS()
postscript(file="SuppFigures/sfigure6a_6c.eps", horizontal = FALSE, onefile = FALSE, paper = "special", family="Helvetica", width=6.7, height=5.4, pointsize=11)
par(mfrow=c(2,2), mgp=c(1.5,.5,0),  mar=c(3.5,4,1,1)) # set figure dimensions
#IID
plot(x[1:200], IIDobs[1:200], type="l",xlab="x",ylab="Observations")
polygon(twoh_y, twoh_x, col=test.colors[1], border=NA)
points(x[1:200], IIDobs[1:200], pch=20, cex=0.5, col="black")
legend("topleft", c("90% CI Heteroskedastic","90% CI Homoskedastic","90% CI IID",
                    "Observations","Errors"), pch=c(15,15,15,20,NA),
       lty=c(NA,NA,NA,NA,1),bty="n", col=c(test.colors[3:1],
                                           "black","red"))
axis(side=4, labels=FALSE)
put.fig.letter("a.",font=2)

plot(x[1:200], H.obs[1:200], type="l",xlab="x",ylab="Observations")
polygon(hettwoh_y, hettwoh_x, col=test.colors[3], border=NA)
axis(side=4, labels=FALSE)
put.fig.letter("c.",font=2)

lines(x[1:200],err_pos, col="red",lty=1)
lines(x[1:200],err_neg, col="red",lty=1)
points(x[1:200], H.obs[1:200], pch=20, cex=0.5, col="black")
axis(side=4, labels=FALSE)

plot(x[1:200], obs[1:200], type="l",xlab="x",ylab="Observations")
polygon(homtwoh_y, homtwoh_x, col=test.colors[2], border=NA)
points(x[1:200], obs[1:200], pch=20, cex=0.5, col="black")
axis(side=4, labels=FALSE)
put.fig.letter("b.",font=2)

# plot(x[1:200], H.obs[1:200], type="l", col="white", xlab="x",ylab="Observations")
# polygon(twoh_y, twoh_x, col="gold", border=NA)
# polygon(hettwoh_y, hettwoh_x, col="cornflowerblue", border=NA)
# polygon(homtwoh_y, homtwoh_x, col="gray", border=NA)

dev.off()

#------------------------------- Sup. Figure 7 (IID version) -----------------------------------#
# Compare the surprise index and the pdfs for when the number of observations changes:
# pdf(file="SuppFigures/nRuckertetal_sup8.pdf", family="Helvetica", height=5.4, width=6.7,pointsize=11)
# par(mfrow=c(2,2),mgp=c(1.5,.5,0), mar=c(3.5,4,4,1))
setEPS()
postscript(file="SuppFigures/sfigure7a_7dIID_version.eps", horizontal = FALSE, onefile = FALSE, paper = "special", family="Helvetica", width=6.7, height=5.4, pointsize=11)
par(mfrow=c(2,2), mgp=c(1.5,.5,0),  mar=c(3.5,4,1,1)) # set figure dimensions

# 1) PDF
plot(one_a, xlab="a", ylab="Probability density", main="",
     ylim=c(0,1), lwd=3, yaxt="n") #iid
lines(five_a, col="gray",lwd=3)
lines(twoh_a, col="gold",lwd=3)
lines(hund_a, col="lightseagreen",lwd=3)
lines(fity_a, col="red",lwd=3)
lines(twent_a, col="blue",lwd=3)
lines(c(-10,-10),c(0,0.18), lwd=2, lty=2)
lines(c(10,10),c(0,0.18), lwd=2, lty=2)
lines(c(-10,10),c(0.18,0.18), lwd=2, lty=2)
put.fig.letter("a.",font=2)
legend("topright", c("Prior","1","5","20","50","100","200"), cex=0.75,
       lty=c(2,1,1,1,1,1,1),lwd=2, col=c("black","black","gray","blue","red",
                                         "lightseagreen","gold"))
# 2) Surprise index
surprise = read.csv("ToyScripts/testhypoth_surprise.csv")
percent = surprise[1:16,1]
twohund = surprise[1:16,7]
onehund = surprise[1:16,9]
fifty = surprise[1:16,11]
twenty = surprise[1:16,13]
five = surprise[1:16,15]
one = surprise[1:16,17]
#par(mgp=c(1.5,.5,0), mar=c(3.5, 3, 4, 2))
plot(percent, percent, typ="l", lty=2, ylab="Percent covered [%]", 
     xlab="Credible interval [%]", lwd=2, ylim=c(0,100))
points(percent, one*100, pch=8, col="black", lwd=2)
points(percent, five*100, pch=4, col="gray", lwd=2)
points(percent, twenty*100, pch=10, col="blue", lwd=2)
points(percent, fifty*100, pch=20, col="red", lwd=2)
points(percent, onehund*100, pch=15, col="lightseagreen", lwd=2)
points(percent, twohund*100, pch=17, col="gold", lwd=2)
put.fig.letter("c.",font=2)
lines(c(88,102),c(102,102))
lines(c(88,88),c(88,102))
lines(c(102,102),c(88,102))
lines(c(88,102),c(88,88))
text(100,92, font=2, "d.")
text(40,60, "underconfident", srt=30)
text(45,25, "overconfident", srt=30)

#par(mgp=c(1.5,.5,0), mar=c(5,4,2.5,1))
plot(one_b, xlab="b", ylab="Probability density", main="",lwd=3,
     ylim=c(0,350), yaxt="n") #iid
lines(five_b, col="gray",lwd=3)
lines(twoh_b, col="gold",lwd=3)
lines(hund_b, col="lightseagreen",lwd=3)
lines(fity_b, col="red",lwd=3)
lines(twent_b, col="blue",lwd=3)
lines(c(-10,-10),c(0,50), lwd=2, lty=2)
lines(c(10,10),c(0,50), lwd=2, lty=2)
lines(c(-10,10),c(50,50), lwd=2, lty=2)
put.fig.letter("b.",font=2)

#par(mgp=c(1.5,.5,0), mar=c(5, 3, 2.5, 2))
plot(percent, percent, type="l", lty=2, ylab="Percent covered [%]", 
     xlab="Credible interval [%]", lwd=2, xlim=c(90,100), ylim=c(90,100))
points(percent, one*100, pch=8, col="black", lwd=2)
points(percent, five*100, pch=4, col="gray", lwd=2)
points(percent, twenty*100, pch=10, col="blue", lwd=2)
points(percent, fifty*100, pch=20, col="red", lwd=2)
points(percent, onehund*100, pch=15, col="lightseagreen", lwd=2)
points(percent, twohund*100, pch=17, col="gold", lwd=2)
put.fig.letter("d.",font=2)
text(94,95, "underconfident", srt=35)
text(95,94, "overconfident", srt=35)
legend("bottomright", c("1","5","20","50","100","200","1:1 line"),
       pch=c(8,4,10,20,15,17,NA), lty=c(NA,NA,NA,NA,NA,NA,2), ncol=2, lwd=2, cex=0.75, 
       col=c("black","gray","blue","red","lightseagreen","gold","black"))
dev.off()

#------------------------------- Sup. Figure 8 -----------------------------------#
# Compare the surprise index and the pdfs for iid, AR(1), and heteroskedastic erros:
#pdf(file="SuppFigures/nRuckertetal_sup9.pdf", family="Helvetica", height=5.4, width=6.7, pointsize=11)
setEPS()
postscript(file="SuppFigures/sfigure8a_8d.eps", horizontal = FALSE, onefile = FALSE, paper = "special", family="Helvetica", width=6.7, height=5.4, pointsize=11)
par(mfrow=c(2,2), mgp=c(1.5,.5,0),  mar=c(3.5,4,1,1)) # set figure dimensions
# 1) PDF
#par(mfrow=c(2,2),mgp=c(1.5,.5,0), mar=c(3.5,4,4,1))
plot(hettwoh_a, xlab="a", ylab="Probability density", main="",
     ylim=c(0,4), lwd=3, col=test.colors[3], yaxt="n", xlim=c(-10,10))
lines(twoh_a, col=test.colors[1],lwd=3)
lines(homtwoh_a, col=test.colors[2],lwd=3)
lines(c(-10,-10),c(0,0.18), lwd=2, lty=3)
lines(c(10,10),c(0,0.18), lwd=2, lty=3)
lines(c(-10,10),c(0.18,0.18), lwd=2, lty=3)
put.fig.letter("a.",font=2)
legend("topright", c("Prior", "1:1 line", "IID","Homoskedastic","Heteroskedastic"), cex=0.75, pt.cex=c(NA,NA,circlesize[1:3]), y.intersp = 1.5,
       lty=c(3,2,1,1,1),pch=c(NA,NA,20,20,20), lwd=2,col=c("black", "black", test.colors[1:3]))

# 2) Surprise index
surprise = read.csv("ToyScripts/testhypoth_surprise.csv")
percent = surprise[1:16,1]
homtwohund = surprise[1:16,3]
hettwohund = surprise[1:16,5]
twohund = surprise[1:16,7]
#par(mgp=c(1.5,.5,0), mar=c(3.5, 3, 4, 2))
plot(percent, percent, typ="l", lty=2, ylab="Percent covered [%]", 
     xlab="Credible interval [%]", lwd=2)
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
plot(homtwoh_b, xlab="b", ylab="Probability density", main="",lwd=3,
     ylim=c(0,600), col=test.colors[2], yaxt="n", xlim=c(-10,10))
lines(hettwoh_b, col=test.colors[3],lwd=3)
lines(twoh_b, col=test.colors[1],lwd=3)
lines(c(-10,-10),c(0,50), lwd=2, lty=3)
lines(c(10,10),c(0,50), lwd=2, lty=3)
lines(c(-10,10),c(50,50), lwd=2, lty=3)
put.fig.letter("b.",font=2)

#par(mgp=c(1.5,.5,0), mar=c(5, 3, 2.5, 2))
plot(percent, percent, type="l", lty=2, ylab="Percent covered [%]", 
     xlab="Credible interval [%]", lwd=2, xlim=c(90,100), ylim=c(90,100))
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
