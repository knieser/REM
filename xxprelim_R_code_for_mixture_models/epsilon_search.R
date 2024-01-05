epsilon_search <- function(X, delta, ueps, intl_mix, intl_mu, intl_sigma){
  
  # number of groups
  k = ncol(intl_mu)
  
  # number of steps for search
  steps = 25
  
  # lower bound for epsilon
  leps = 0
  
  # initialize
  eps_range = rep(0, steps)
  chk = rep(1, steps)
  
  # binary search
  for (iter in 1:steps){
    message('.....step ', iter, ' out of ', steps)
    # take midpoint between lower and upper epsilon
    eps_range[iter] = (leps + ueps)/2
    
    # get REM estimates
    REM_output <- RobustEMAlg(X, k, eps_range[iter], intl_mix, intl_mu, intl_sigma)
    mix2 = REM_output$mix
    mu2 = REM_output$mu
    sigma2 = REM_output$sigma
    
    # check heuristic
    chk[iter] = checkEps(mix2, mu2, sigma2, eps_range[iter], 0.9)
    
    # update upper or lower epsilon
    if(chk[iter] < delta){
      leps = eps_range[iter]
    } else{
        ueps = eps_range[iter]}
  }
  
  opt_eps = (leps + ueps) / 2
  message(paste0('leps = ', leps,
                 '; ueps = ', ueps))
  return(opt_eps)
}