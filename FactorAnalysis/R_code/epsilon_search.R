epsilon_search <- function(X, delta, ueps, intl_nu, intl_lambda, intl_psi){
  
  # number of factors
  k = ncol(intl_lambda)
  
  # number of steps for search
  steps = 25
  
  # lower bound for epsilon
  leps = 0
  
  # initialize
  eps_range = rep(0, steps)
  chk = rep(1, steps)
  
  # binary search
  for (iter in 1:steps){
    message('...step ', iter)
    # take midpoint between lower and upper epsilon
    eps_range[iter] = (leps + ueps)/2
    
    # get REM estimates
    REM_output <- RobustEMAlg(X, k, eps_range[iter], intl_nu, intl_lambda, intl_psi)
    nu2 = REM_output$nu
    lambda2 = REM_output$lambda
    psi2 = REM_output$psi
    
    # check heuristic
    chk[iter] = checkEpsFA(nu2, lambda2, psi2, eps_range[iter], 0.9)
    
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