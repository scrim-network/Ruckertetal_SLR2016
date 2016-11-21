  The effects of time-varying observation errors on semi-empirical sea-level projections

  Reference:
  Ruckert, KL, Guan, Y, Bakker, AMR, Forest, FE, and Keller, K. The effects of time-varying observation errors on semi-empirical sea-level projections, (accepted at Climatic Change).

  This program is distributed in the hope that it will be useful,
  but WITH NO WARRANTY (NEITHER EXPLICITNOR IMPLICIT). We are not liable
  for the behavior of these codes in your own application. You are free
  to share this code so long as the authors(s) and version history remain
  intact.

    Kelsey Ruckert, klr324@psu.edu
    Yawen Guan, yig5031@psu.edu
            
            Using the Vermeer & Rahmstorf 2009 model
 ===================================================================
| For the impatient:(NOTE:This program runs for roughly 5 hours)    |
|  1. cd VR_calibrate                                               |
|  2. makedir Workspace                                             |
|  3. Start R                                                       |
|  4. Type source(“cali_homo_model_AR_VR09.R”)                      |
|  5. Type source(“cali_heter_model_AR_VR09.R”)                     |
|  6. Type source(“bootstrap_VR09.R”)                               |
|  7. Type source(“PlotRuckert_etal_VR.R”)                          |
 ===================================================================
 
             Using the Grinsted et al. 2010 model
  ===================================================================
| For the impatient:(NOTE:This program runs for roughly 5 hours)    |
|  1. cd Grinsted_calibrate                                         |
|  2. makedir Workspace                                             |
|  3. Start R                                                       |
|  4. Type source(“cali_homo_model_AR_grinsted.R”)                  |
|  5. Type source(“cali_heter_model_AR_grinsted.R”)                 |
|  6. Type source(“bootstrap_grinsted.R”)                           |
|  7. Type source(“PlotRuckert_etal_Grinsted.R”)                    |
 ===================================================================

            Using the Rahmstorf 2007 model
 ===================================================================
| For the impatient:(NOTE:This program runs for roughly 5 hours)    |
|  1. cd Rahmstorf_calibrate                                        |
|  2. makedir Workspace                                             |
|  3. Start R                                                       |
|  4. Type source(“Mega_Rahmstorf.R”)                               |
|  5. Type source(“PlotRuckert_etal.R”)                             |
 ===================================================================

Required software:
    R (used R 3.2.1)

Required libraries:
    mcmc
    ncdf or ncdf4
    coda
    RColorBrewer
    mvtnorm
    DEoptim
    compiler

   Please note that scripts may need to be changed to point to the correct 
   location of the files/script in order to run.

 =============================================================================
| Short description of main R scripts:
|  The main scripts were used in the analysis of the paper.
|
|  1. bootstrap_MODELNAME.R: Bootstrap fitting of a SLR model.
|  2. cali_homo_model_AR_MODELNAME.R: Homoskedastic AR1 Markov Chain Monte Carlo (MCMC) fitting of  
|     a SLR model.
|  3. cali_heter_model_AR_MODELNAME.R: Heteroskedastic AR1 MCMC fitting of a SLR model.      
|  4. Rar.R: Estimates the log likelihood of the AR1 process and simulates the lag-1 
|     autocorrelation coefficient.
|  5. MODELNAMEobs_likelihood_AR.R: Estimates the AR1 likelihood function for the MCMC methods.
|  6. Estimate_surprise_index.R: Script/function used in part to help create the MODELNAME_surprise.csv files
|  7. RCP_temp_scenarios.R: Fits the Rahmstorf 2007 SLR model with Bootstrap and MCMC methods using 
|     RCP 2.6, 4.5, 6.0, and 8.5.
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
|     two-parameter linear model to the observations with MCMC. 
|  3. Toy_MCMC_Test2b.R: Generates 200 heteroskedastic AR observations and fits a simple 
|     two-parameter linear modelmodel to the observations with MCMC.
|  4. Toy_TestMCMC_plots.R: Loads in the results from test scripts 1, 2a, and 2b to test 
|     why the surprise index fails in the low probabilities. It produces S. Fig. 8.
|  5. Toy_bootstrap.R: This program will run and plot the bootstrap method for a simple 
|     two-parameter linear model with AR1 residuals. 
|  6. testAR1vsIID.R: This program assesses the impact of AR1 representation on model calibration.  
 =============================================================================

Credits:
    SLR Models: 
    Rahmstorf S (2007) A Semi-empirical approach to projecting future sea-level rise 
        Science 315(5810) 368–370, doi:10.1126/science.1135456.
        
    Grinsted A, Moore JC, and Jevrejeva S (2010) Reconstructing sea level from paleo and projected temperatures 200 to 2100 A.D. 
        Clim Dyn 34(4):461–472. doi:10.1007/s00382-008-0507-2
        
    Vermeer M and Rahmstorf S (2009) Global sea level linked to global temperature. Proceedings of the National 
	      Academy of Sciences 106(51):21527-21532, doi:10.1073/pnas/0907765106



