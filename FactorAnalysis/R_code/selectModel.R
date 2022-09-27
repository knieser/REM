selectModel <- function(X, k_range, delta){
  
  # data dimensions
  n = nrow(X)
  p = ncol(X)
  
  AIC = k_range
  BIC = k_range
  for (k in k_range){
    
    # REM estimates
    REM_output = REM_main(X, k, delta)
    lik_rem = REM_output$lik_rem
    
    # count parameters
    # nu (p) + lambda (k*p) + psi (p*(p+1)/2)
    parms = p + k*p + p*(p+1)/2
    
    AIC[k] = -2*lik_rem + 2*parms
    BIC[k] = -2*lik_rem + log(n)*parms
  }
  
  msg1 = paste('Number of factors based on AIC: ', AIC[which(AIC==min(AIC))])
  msg2 = paste('Number of factors based on BIC: ', BIC[which(BIC==min(BIC))])
  message(msg1)
  message(msg2)
  
  output = list(AIC = AIC, BIC = BIC)
  return(output)
}