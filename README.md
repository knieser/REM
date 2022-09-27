# Robust Expectation-Maximization (REM)

Individuals differ in many substantive ways that are not always captured through the assumed data-generating model. In an effort to address this, we present a robust estimation procedure based on the EM algorithm which we call REM (robust expectation-maximization).

The files in this repository contain the MATLAB source code used for simulation studies found in the paper, Addressing Heterogeneous Populations in Latent Variable Settings through Robust Estimation (paper link: https://doi.org/10.1037/met0000413). Further details about REM and the simulation studies can be found there. 


## MixtureModel
This folder contains code for running simulations
comparing EM and REM estimation of Gaussian mixture model parameters
under model misspecification. Simulations can be run from MixtureModel/GMM_sim_main.m 

Estimation can be run for an input  p-by-n dataset X from MixtureModel/GMM_estimates.m



## FactorAnalysis
This folder contains code for running simulations
comparing EM and REM estimation of factor structures
within heterogeneous populations. Simulations can be run from FactorAnalysis/FA_sim_main.m

Estimation can be run for an input  p-by-n dataset X from FactorAnalysis/FA_estimates.m








