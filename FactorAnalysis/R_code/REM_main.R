REM_main <- function(X, k_range, delta = 0.05){
  
  source('REM_estimates.R')
  source('EMAlg.R')
  source('epsilon_search.R')
  source('checkEpsFA.R')
  source('RobustEMAlg.R')
  
  # data dimensions
  n = nrow(X)
  p = ncol(X)
  
  REM_output = vector(mode = 'list', length = length(k_range))
  
  for (k in 1:length(k_range)){
    message('..working on k = ', k_range[k])
    REM_output[[k]] = REM_estimates(X, k_range[k], delta)
    message('..done with k = ', k_range[k])
  }
  return(REM_output)
}