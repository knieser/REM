REM_main <- function(X, k, delta){
  
  # EM estimates
  message('...getting EM estimates')
  EM_output <- EMAlg(X, k)
  nu1 = EM_output$nu
  lambda1 = EM_output$lambda
  psi1 = EM_output$psi
  ind_lik = EM_output$ind_lik
  
  # set max epsilon
  ueps = quantile(exp(ind_lik), 0.3)
  
  # epsilon search
  message('...searching for optimal epsilon')
  opt_eps = epsilon_search(X, delta, ueps, nu1, lambda1, psi1)
  
  # REM algorithm
  message('...getting REM estimates')
  REM_output = RobustEMAlg(X, k, opt_eps, nu1, lambda1, psi1)
  
  output = list(
    EM_output = EM_output,
    REM_output = REM_output,
    epsilon = opt_eps)
  return(output)
  
}