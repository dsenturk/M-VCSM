# M-VCSM
CONTENTS OF THIS FOLDER ——————————————

M_VCSM_tutorial.R :A step-by-step implementation of M-VCSM and the associated procedures described in "Multilevel Varying Coefficient Spatiotemporal Model".

M_VCSM_simulation.R : Function for simulating one data set under the simulation design described in Section 4.1.

M_VCSM_decomposition.R : Function for FPCA decomposition (Estimation Algorithm steps 1-4 in Section 2.2) of M-VCSM described in "Multilevel Varying Coefficient Spatiotemporal Model" including estimation of multilevel eigenfunctions and eigenvalues. 

M_VCSM_MCMC.R : Function for MCMC estimation (Estimation Algorithm steps 5-6) of M-VCSM model described in "Multilevel Varying Coefficient Spatiotemporal Model", including estimation of varying coefficient functions (VCFs), spatial variance parameters, region-specific variances, measurement error variance and region- and facility-specific PC scores. 

M_VCSM_inference.R : Function for obtaining inference for varying coefficient functions, time-static parameter estimates and prediction of region- and facility-specific deviations.


INTRODUCTION ——————————————

The contents of this folder allow for implementation of the M-VCSM estimation and inference described in "Multilevel Varying Coefficient Spatiotemporal Model". Users can simulate a sample data frame (M_VCSM_simulation.R) and apply the proposed estimation algorithm (M_VCSM_decomposition.R, M_VCSM_MCMC.R). Also, we include tools to perform inference on multilevel varying coefficient functions and prediction of region- and facility-specific deviations (MST_FM_inference.R), allowing users to obtain effect of region- and facility-level risk factors as well as their simultaneous confidence bands. Detailed instructions on how to perform the aforementioned procedures and visualize results are included in M_VCSM_tutorial.R.

REQUIREMENTS ——————————————

The included R programs require R 4.0.2 (R Core Team, 2020) and the packages listed in M_VCSM_tutorial.R.

INSTALLATION ——————————————

Load the R program files into the global environment and install required packages using commands in M_VCSM_tutorial.R
