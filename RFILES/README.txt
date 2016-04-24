  Improving the statistical method can raise the upper tail of sea-level projections Codes

  Reference:
  Ruckert, KL, Guan, Y, Forest, FE, and Keller, K. Improving the statistical
  method can raise the upper tail of sea-level projections, (submitted to ERL).

  This program is distributed in the hope that it will be useful,
  but WITH NO WARRANTY (NEITHER EXPLICITNOR IMPLICIT). We are not liable
  for the behavior of these codes in your own application. You are free
  to share this code so long as the authors(s) and version history remain
  intact.

    Kelsey Ruckert, klr324@psu.edu
    Yawen Guan, yig5031@psu.edu

 ===================================================================
| For the impatient:(NOTE:This program runs for roughly 5 hours)    |
|  1. Start R                                                       |
|  2. Type source(“Mega_Rahmstorf.R”)                               |
|  3. Type source(“PlotRuckert_etal.R”)                             |
|  4. Type source(“Toy_TestMCMC_plot.R”)                            |
 ===================================================================

Required software:
    R

Required libraries:
    mcmc
    ncdf
    coda
    mvtnorm
    DEoptim
    compiler

   Please note that the folder directory MUST be in the same format as when downloaded
   otherwise the scripts will not locate the files/scripts needed to run.

 =============================================================================
| Short description of main R scripts:
|  The main scripts were used in the analysis of the paper.
|
|  1. Mega_Rahmstorf.R: Script that runs all three model fitting methods when sourced.
|  2. bootstrap_Rahm_new.R: Bootstrap fitting of the Rahmstorf SLR model.
|  3. Rcali_homo_model_AR.R: Homoskedastic AR1 Markov Chain Monte Carlo (MCMC) fitting of  
|     the Rahmstorf SLR model.
|  4. Rcali_heter_model_AR.R: Heteroskedastic AR1 MCMC fitting of the Rahmstorf SLR model.      
|  5. converge_test.R: Tests MCMC convergence based on the Potential Scale Reduction 
|     factor.
|  6. RCP_temp_scenarios.R: Fits the SLR model with Bootstrap and MCMC methods using 
|     RCP 2.6, 4.5, 6.0, and 8.0 and plots S. Fig. 4.
|  7. Rar.R: Estimates the log likelihood of the AR1 process and simulates the lag-1 
|     autocorrelation coefficient.
|  8. Robs_likelihood_AR.R: Estimates the AR1 likelihood function for the MCMC methods.
 =============================================================================



 =============================================================================
| Short description of the toy scripts:
|  The toy scripts are meant to show/ teach how to fit a model using MCMC or the bootstrap 
|  method.
|  Additionally, they were used to test bias in the surprise index. The program requires 
|  the user to designate the length of observations and the assumptions, so the scripts 
|  can not be sourced.
|  
|  1. Toy_MCMC_Test1.R: Generates a User defined number of iid observations and fits a 
|     simple model to the observations with MCMC. 
|  2. Toy_MCMC_Test2a.R: Generates 200 homoskedastic AR observations and fits a simple 
|     two-parameter linear modelmodel to the observations with MCMC. 
|  3. Toy_MCMC_Test2b.R: Generates 200 heteroskedastic AR observations and fits a simple 
|     two-parameter linear modelmodel to the observations with MCMC.
|  4. Toy_TestMCMC_plots.R: Loads in the results from test scripts 1, 2a, and 2b to test 
|     why the surprise index fails in the low probabilities. It produces S. Fig. 6-8.
|  5. Toy_bootstrap.R: This program will run and plot the bootstrap method for a simple 
|     two-parameter linear model with AR1 residuals. 
 =============================================================================

Credits:
    SLR Model: 
    Rahmstorf S (2007) A Semi-empirical approach to projecting future sea-level rise 
        Science 315(5810) 368–370, doi:10.1126/science.1135456.


