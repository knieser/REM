# Robust Expectation-Maximization (REM)

Individuals differ in many substantive ways that are not always captured through the assumed data-generating model. In an effort to address this, we present a robust estimation procedure based on the EM algorithm which we call REM (robust expectation-maximization).

The files in this repository contain the MATLAB source code used for simulation studies found in the paper, Addressing Heterogeneous Populations in Latent Variable Settings through Robust Estimation. Further details about REM and the simulation studies can be found there. 


## MixtureModel
This folder contains code for running simulations
comparing EM and REM estimation of Gaussian mixture model parameters
under model misspecification.

Includes:
- GMM_sim_main
- GMM_data_sim
- GMM_estimates
- globalRobustEMAlgMixture
- RobustEMAlgMixture
- EMAlgMixture
- checkEps
- makeFigures
- PrettyFig
- printFig1
- printFig2
- printFig3


## FactorAnalysis
This folder contains code for running simulations
comparing EM and REM estimation of factor structures
within heterogeneous populations.

Includes:
- FA_sim_main
- FA_data_sim
- SigmaSim
- SigmaSimNormals
- computeCongruence
- FA_sample_data
- FA_estimates
- EMAlg
- RobustEMAlg
- checkEpsFA
- aggregate
- PrettyFig
- printFig4
- printFig5







