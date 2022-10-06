REM_estimates <- function(X, k, delta, nstart){
  
  n = nrow(X)
  p = ncol(X)
  
  # initial guess of mixture proportions and sigma using kmeans
  kfit <- kmeans(X, k, nstart = 1e4)
  intl_mix <- kfit$size / n  
  intl_sigma <- array(data = NA, dim = c(p,p,k))
  for (j in 1:k){
    intl_sigma[,,j] = cov(X[kfit$cluster==j,])
  }
  
  # initial starts for mu
  intl_mu = t(X[sample(nrow(X), size = k*nstart, replace = T),])
  intl_mu = cbind(t(kfit$centers), intl_mu)
  intl_mu = array(data = intl_mu, dim = c(p, k, nstart+1))
  
  max_lik = -Inf
  max_lik_rem = -Inf
  
  # EM estimates
  message('...getting EM estimates')
  for (s in 1:nstart){
    local_EM_output <- EMAlg(X, k, intl_mix, intl_mu[,,s], intl_sigma)
    lik = sum(local_EM_output$ind_lik)
    
    if (lik > max_lik){
      message('.....new max reached')
      max_lik = lik
      EM_output <- local_EM_output
    }
  }
  
  # REM estimates
  message('...getting REM estimates')
  message('...starting REM global optimization')
  ind_lik = EM_output$ind_lik
  eps_start = quantile(exp(ind_lik), 0.05)
  for (s in 1:nstart){
    local_REM_output <- RobustEMAlg(X, k, eps_start, intl_mix, intl_mu[,,s], intl_sigma)
    lik_rem = sum(local_REM_output$ind_lik_rem)
    
    if (lik_rem > max_lik_rem){
      message('.....new max reached')
      max_lik_rem = lik_rem
      REM_start <- local_REM_output
    }
  }
  mu_start = REM_start$mu
  
  # epsilon search
  message('...searching for optimal epsilon')
  ueps = quantile(exp(ind_lik), 0.3)
  opt_eps = epsilon_search(X, delta, ueps, intl_mix, mu_start, intl_sigma)
  
  # REM estimates at optimal epsilon
  message('...getting REM estimates at optimal epsilon')
  REM_output = RobustEMAlg(X, k, opt_eps, intl_mix, mu_start, intl_sigma)
  
  # calculate AIC and BIC
  parms = k*p + k*p*(p+1)/2 + k-1
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