REM_estimates <- function(X, k, delta){
  
  n = nrow(X)
  p = ncol(X)
  
  # EM estimates
  message('...getting EM estimates')
  EM_output <- EMAlg(X, k)
  nu1 = EM_output$nu
  lambda1 = EM_output$lambda
  psi1 = EM_output$psi
  ind_lik = EM_output$ind_lik
  
  # epsilon search
  message('...searching for optimal epsilon')
  ueps = quantile(exp(ind_lik), 0.3)
  opt_eps = epsilon_search(X, delta, ueps, nu1, lambda1, psi1)
  
  # REM algorithm
  message('...getting REM estimates')
  REM_output = RobustEMAlg(X, k, opt_eps, nu1, lambda1, psi1)
  
  # calculate AIC and BIC
  parms = p + k*p + p*(p+1)/2
  lik_rem = REM_output$lik_rem
  AIC_rem = -2*lik_rem + 2*parms
  BIC_rem = -2*lik_rem + log(n)*parms
  
  output = list(
    k = k,
    EM_output = EM_output,
    REM_output = REM_output,
    epsilon = opt_eps,
    AIC_rem = AIC_rem,
    BIC_rem = BIC_rem)
  return(output)
  
}