# M-VCSM
CONTENTS OF THIS FOLDER ——————————————

M_VCSM_tutorial.R : A step-by-step implementation of MST-FM and the associated procedures described in "Multilevel Modeling of Spatially Nested Functional Data: Spatiotemporal Patterns of Hospitalization Rates in the U.S. Dialysis Population".

M_VCSM_simulation.R : Function for simulating one data set under the simulation design described in Web Appendix E of the supplementary materials.

M_VCSM_decomposition.R : Function for FPCA decomposition (estimation steps 1-5 in Table 1) described in "Multilevel Modeling of Spatially Nested Functional Data: Spatiotemporal Patterns of Hospitalization Rates in the U.S. Dialysis Population", including estimation of mean function, multilevel eigenfunctions and eigenvalues.

M_VCSM_MCMC.R : Function for MCMC estimation (estimation step 6 in Table 1) of MST-FM model described in "Multilevel Modeling of Spatially Nested Functional Data: Spatiotemporal Patterns of Hospitalization Rates in the U.S. Dialysis Population", including estimation of spatial variance parameters, measurement error variance and region- and facility-specific PC scores.

M_VCSM_inference.R : Function for obtaining prediction and inference for multilevel trajectories described in Section 2.2 of "Multilevel Modeling of Spatially Nested Functional Data: Spatiotemporal Patterns of Hospitalization Rates in the U.S. Dialysis Population".


INTRODUCTION ——————————————

The contents of this folder allow for implementation of the MST-FM estimation and inference described in "Multilevel Modeling of Spatially Nested Functional Data: Spatiotemporal Patterns of Hospitalization Rates in the U.S. Dialysis Population". Users can simulate a sample data frame (MST_FM_simulation.R) and apply the proposed estimation algorithm (MST_FM_decomposition.R, MST_FM_MCMC.R). Also, we include tools to perform prediction and inference on multilevel hospitalization rate trajectories (MST_FM_inference.R), allowing users to obtain region- and facility-specific predicted hospitalization rate trajectories as well as their pointwise confidence intervals. Further, for applications with small number of regions, users can obtain corrected prediction and confidence intervals for region-specific trajectories for a more precise inference (MST_FM_correctedinference.R). Detailed instructions on how to perform the aforementioned procedures, make predictions of region- and facility-level hospitalization rate trajectories and visualize results are included in MST_FM_tutorial.R.

REQUIREMENTS ——————————————

The included R programs require R 3.5.3 (R Core Team, 2019) and the packages listed in MST_FM_tutorial.R.

INSTALLATION ——————————————

Load the R program files into the global environment and install required packages using commands in MST_FM_tutorial.R
