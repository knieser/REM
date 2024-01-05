RobustEMAlg <- function(X, k, epsilon, intl_mix, intl_mu, intl_sigma, min_var = 0.05){
  
  # data dimensions
  n = nrow(X)
  p = ncol(X)
  
  # initial guesses
  mix = intl_mix
  mu = intl_mu
  sigma = intl_sigma
  gamma = 0.5
  
  # tolerance parameters
  tol = 1e-6
  maxiter = 1e3
  
  # REM algorithm
  cond_lik = matrix(data = 0, nrow = n, ncol = k)
  obj = rep(-Inf, maxiter)
  exit = 0
  for (iter in 2:maxiter){
    for (j in 1:k){
      # check minimum variance
      diag(sigma[,,j])[diag(sigma[,,j]) < min_var] <- min_var
      
      # calculate log-likelihoods
      tryCatch({chol(sigma[,,j])},
               error = function(e){
                 message('could not apply cholesky to sigma')
                 exit = 1
               })
      if (exit == 1){break}
      
      V = chol(sigma[,,j])
      inv_V = solve(V)
      logdetV = log(det(V))
      Z = sweep(X, 2, mu[,j], "-") %*% inv_V
      cond_lik[,j] = -(1/2)*(p*log(2*pi) + 2*logdetV + apply(Z, 1, function(y) t(y) %*% y))
    }
    
    if (exit == 1){break}
    omega_num = sweep(exp(cond_lik), 2, mix, "*")
    ind_lik = log(apply(omega_num, 1, sum))

    # calculate modified log-likelihoods
    ind_lik_rem = log(gamma * exp(ind_lik) + (1-gamma) * epsilon)
    
    # calculate value for objective function
    obj[iter] = sum(ind_lik_rem)
    
    # check for convergence
    if(abs(obj[iter] - obj[iter-1]) < tol*abs(obj[iter-1])){break}
    
    # compute weights
    weights_num = log(gamma) + ind_lik
    weights_denom = ind_lik_rem
    log_weights = weights_num - weights_denom
    weights = exp(log_weights)
    
    # estimate gamma
    gamma = mean(weights)
    
    # estimate omega
    omega = omega_num / exp(ind_lik)
    
    # estimate mu, sigma, and mix
    for (j in 1:k){
      total_weight = weights * omega[,j]
      mu[,j] = (t(X) %*% total_weight) / sum(total_weight)
      X_centered = sweep(X, 2, mu[,j], "-")
      sigma[,,j] = ( t(X_centered) %*% (total_weight * X_centered) ) / sum(total_weight)
      mix[j] = sum(total_weight) / sum(weights)
    }
    
    if(iter==maxiter){message(
      paste('maxiter reached; increase maxiter', abs(obj[iter] - obj[iter-1])/abs(obj[iter-1])))}
  }
  
  output = list(mix = mix, 
                mu = mu, 
                sigma = sigma,
                posterior = omega,
                gamma = gamma,
                weights = weights,
                ind_lik = ind_lik,
                ind_lik_rem = ind_lik_rem,
                lik = sum(ind_lik),
                lik_rem = sum(ind_lik_rem))
  return(output)
}