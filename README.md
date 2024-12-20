# Robust Expectation-Maximization (REM)

Individuals differ in many substantive ways that are not always captured through the assumed data-generating model. In an effort to address this, we present a robust estimation procedure based on the EM algorithm which we call REM (robust expectation-maximization). 

## R package

<!-- badges: start -->  
[![CRAN status](https://www.r-pkg.org/badges/version/REMLA)](https://CRAN.R-project.org/package=REMLA) 
[![](https://cranlogs.r-pkg.org/badges/grand-total/REMLA)](https://cran.r-project.org/package=REMLA)  
<!-- badges: end -->

This folder contains R package files to run exploratory and confirmatory factor analyses (EFAs and CFAs) using the REM algorithm. The following code can be used to download the latest version of the package to your RStudio from Github.

``` r
library(devtools)
devtools::install_github('knieser/REM/R_package')
library(REMLA)
```

### Vignette

[REMLA tutorial](https://knieser.github.io/REM/REMLA_tutorial.html)


## Psychological Methods paper
This folder contains the MATLAB code used for simulation studies found in the paper, [Addressing Heterogeneous Populations in Latent Variable Settings through Robust Estimation](https://doi.org/10.1037/met0000413 "https://doi.org/10.1037/met0000413"). Further details about REM and the simulation studies can be found there. 

In MATLAB, with an input p-by-n dataset X, REM estimation for 
- Gaussian mixture models can be run from `Psych_Methods_paper/MixtureModel/GMM_estimates.m` 
- Exploratory factor analysis can be run from `Psych_Methods_paper/FactorAnalysis/FA_estimates.m`

Simulations can be run from `Psych_Methods_paper/MixtureModel/GMM_sim_main.m` and `Psych_Methods_paper/MixtureModel/FA_sim_main.m` for mixture models and factor models, respectively


## Gaussian mixture modeling 
The code in `xxprelimR_code_for_mixture_models/` folder adds functions to run REM estimation with Gaussian mixture models and is currently under development.







