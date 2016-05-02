##################################################################################################
#
#  -file = "convergence_test.R"   Code written November 2015 edited April 2016
#  - Author: Kelsey Ruckert (klr324@psu.edu)
#
#  -This program loads in the workspace generated from the running the methods with different seeds.
#       In this script the potental scale reduction factor is used to check convergence of the MCMC methods.
#       For further description and references, please read the paper.
#
# THIS CODE IS PROVIDED AS-IS WITH NO WARRANTY (NEITHER EXPLICIT
# NOT IMPLICIT).  I SHARE THIS CODE IN HOPES THAT IT IS USEFUL,
# BUT I AM NOT LIABLE FOR THE BEHAVIOR OF THIS CODE IN YOUR OWN
# APPLICATION.  YOU ARE FREE TO SHARE THIS CODE SO LONG AS THE
# AUTHOR(S) AND VERSION HISTORY REMAIN INTACT.
#
#
###################################################################################################
rm(list =ls()) #Clear global environment

####################################### Homoskedastic Convergence ################################################################
#Load in workspaces saved from using multiple seeds
load("Workspace/homo1234.RData") # seed 1234
load("Workspace/homo1014.RData") # seed 1014
load("Workspace/homo111.RData")  # seed 111
load("Workspace/homo1.RData")    # seed 1
load("Workspace/homo1780.RData") # seed 1780

# mcmc is converged when the potental scale reduction factor is less than 1.1
homo = as.mcmc(homoskchain1234)
homo2 = as.mcmc(homoskchain111)
homo3 = as.mcmc(homoskchain1014)
homo4 = as.mcmc(homoskchainone)
homo5 = as.mcmc(homoskchain1780)

homolist = mcmc.list(list(homo, homo2, homo3, homo4, homo5))
gelman.diag(homolist)

# pdf(file="homo_gelman.pdf", family="Helvetica",pointsize=13)
# gelman.plot(homolist, xlab="Last iteration in chain", ylab="Shrink Factor")
# dev.off()

####################################### Heteroskedastic Convergence ################################################################
#Load in workspaces saved from using multiple seeds
load("Workspace/heter1780.RData") # seed 1780
load("Workspace/heter1234.RData") # seed 1234
load("Workspace/heter1014.RData") # seed 1014
load("Workspace/heter111.RData")  # seed 111
load("Workspace/heter1.RData")    # seed 1

# mcmc is converged when the potental scale reduction factor is less than 1.1
heter = as.mcmc(prechain1014)
heter2 = as.mcmc(prechain111)
heter3 = as.mcmc(prechain1780)
heter4 = as.mcmc(prechain1234)
heter5 = as.mcmc(prechainone)

heterlist = mcmc.list(list(heter, heter2, heter3, heter4, heter5))
gelman.diag(heterlist)

# gelman.plot(heterlist)

####################################### END ################################################################


