#Global sea-level rise code for Ruckert et al. (accepted)

README file last updated by Kelsey Ruckert, klr324-at-psu-dot-edu, Mon Nov 21 15:05:51 EST 2016

##Citation

This code is intended to accompany the results of

>Ruckert, KL, Guan, Y, Bakker AMR, Forest, FE, and Keller, K. The effects of time-varying observation errors on semi-empirical sea-level projections, (accepted at Climatic Change).

##Overview

This code requires R with the following libraries:
- adaptMCMC
- mcmc
- ncdf
- coda
- RColorBrewer
- mvtnorm
- DEoptim
- compiler

This R code is intended to help users who wish to work with the sea-level rise projections or methods shown in Ruckert et al. (accepted) in greater detail than provided in the appendix of the text. Key functionality of these scripts include:

1. Global sea-level rise projections from 1880 to 2300 with associated probabilities
2. How to fit a model to observations with AR1 residuals using Markov Chain Monte Carlo and Bootstrap
3. Produces plots from the paper

The RFILES directory contains all the scripts and data necessary to run the analysis along with a README file. _(Note that the user may have to edit the scripts according to their folder directory so that the scripts will locate the files/scripts needed to run. Additionally, the following empty folders need to be created before running the analysis: 'Workspace', 'Figures', 'SuppFigures', and 'ToyFigures'. The 'Workspace' must exist in each 'calibrate' folder. Output will be saved to these folders.)_

Instructions on how to run the scripts can be found in the README file in the RFILES directory. The main files for running the analysis with the Vermeer and Rahmtstorf (2009) model can be found in the VR_calibrate folder. The files to run with the Grinsted et al. (2010) model are in the Grinsted_calibrate folder and the main files for the Rahmstorf (2007) analysis are in the Rahmstorf_calibrate folder.

All of the **Toy** scripts are meant to show/teach how to fit a model using MCMC or the bootstrap method. Additionally, they were used to test bias in the surprise index. These files should be run line by line instead of sources. Several files require the user to designate the length of observations.

##Contact
Kelsey Ruckert: <klr324@psu.edu>  
Corresponding author: Klaus Keller at <klaus@psu.edu>
