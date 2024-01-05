EMAlg <- function(X, k, intl_mix, intl_mu, intl_sigma, min_var = 0.05){
  
  # data dimensions
  n = nrow(X)
  p = ncol(X)

  # initial guesses
  mix = intl_mix
  mu = intl_mu
  sigma = intl_sigma
  
  # tolerance parameters
  tol = 1e-6
  maxiter = 1e3

  # EM algorithm
  cond_lik = matrix(data = 0, nrow = n, ncol = k)
  obj = rep(-Inf, maxiter)

  for (iter in 2:maxiter){
    for (j in 1:k){
      # check minimum variance
      diag(sigma[,,j])[diag(sigma[,,j]) < min_var] <- min_var
      
      # calculate log-likelihoods
      V = chol(sigma[,,j])
      inv_V = solve(V)
      logdetV = log(det(V))
      Z = sweep(X, 2, mu[,j], "-") %*% inv_V
      cond_lik[,j] = -(1/2)*(p*log(2*pi) + 2*logdetV + apply(Z, 1, function(y) t(y) %*% y))
    }
    omega_num = sweep(exp(cond_lik), 2, mix, "*")
    ind_lik = log(apply(omega_num, 1, sum))
    
    # calculate value for objective function
    obj[iter] = sum(ind_lik)
    
    # check for convergence
    if(abs(obj[iter] - obj[iter-1]) < tol*abs(obj[iter-1])){break}
    
    # estimate omega
    omega = omega_num / exp(ind_lik)
    
    # estimate mu, sigma, and mix
    for (j in 1:k){
      omega_sum = sum(omega[,j])
      mu[,j] = (t(X) %*% omega[,j]) / omega_sum
      X_centered = sweep(X, 2, mu[,j], "-")
      sigma[,,j] = ( t(X_centered) %*% (omega[,j] * X_centered) ) / omega_sum
      mix[j] = omega_sum / n
    }
    
    if(iter==maxiter){message(
      paste('maxiter reached; increase maxiter', abs(obj[iter] - obj[iter-1])/abs(obj[iter-1])))}
  }
 
  output = list(mix = mix, 
                mu = mu, 
                sigma = sigma,
                posterior = omega,
                ind_lik = ind_lik,
                lik = sum(ind_lik))
  return(output)
}