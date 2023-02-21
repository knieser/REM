RobustEMAlg <- function(X, k, epsilon, intl_nu, intl_lambda, intl_psi){
  
  # data dimensions
  n = nrow(X)
  p = ncol(X)
  
  # initial guess of lambda and psi
  nu = intl_nu
  lambda = intl_lambda
  psi = intl_psi
  gamma = 0.5
  
  # tolerance parameters
  tol = 1e-6
  maxiter = 1e3
  
  # REM algorithm
  obj = rep(-Inf, maxiter)
  for (iter in 2:maxiter){
    
    # calculate log-likelihoods
    sigma = lambda %*% t(lambda) + diag(psi)
    V = chol(sigma)
    inv_V = solve(V)
    logdetV = log(det(V))
    Z = sweep(X, 2, nu, "-") %*% inv_V
    ind_lik = -(1/2)*(p*log(2*pi) + 2*logdetV + apply(Z, 1, function(y) t(y) %*% y))
    
    # calculate modified log-likelihoods
    ind_lik_rem = log(gamma*exp(ind_lik) + (1-gamma)*epsilon)
    
    # calculate value for objective function
    obj[iter] = sum(ind_lik_rem)
    
    # check for convergence
    if(abs(obj[iter] - obj[iter-1]) < tol*abs(obj[iter-1])){break}
    
    # compute weights
    weights_num = log(gamma) + ind_lik
    weights_denom = ind_lik_rem
    log_weights = weights_num - weights_denom
    weights = exp(log_weights)
    
    # set lower bound for weights
    weights[weights < 1e-30] <- 1e-30
    
    # estimate gamma
    gamma = mean(weights)

    # estimate nu
    nu = apply(weights*X, 2, sum) / (n*gamma)
    
    # estimate weighted Cxx
    X_centered = sweep(X, 2, nu, "-")
    Cxx = t(X_centered) %*% (weights * X_centered) / (n*gamma)
    
    # estimate beta
    beta = (t(lambda) %*% inv_V) %*% t(inv_V)
    
    # estimate lambda and psi
    Ezz = diag(k) - beta %*% lambda + beta %*% Cxx %*% t(beta)
    chol_Ezz = chol(Ezz)
    inv_chol_Ezz = solve(chol_Ezz)
    lambda = ((Cxx %*% t(beta)) %*% inv_chol_Ezz) %*% t(inv_chol_Ezz)
    psi = diag((diag(p) - lambda %*% beta) %*% Cxx)
    
    if(iter==maxiter){message(
      paste('maxiter reached; increase maxiter', abs(obj[iter] - obj[iter-1])/abs(obj[iter-1])))}
    
  }
  
  # store joint log-likelihood values
  lik = sum(ind_lik)
  lik_rem = sum(ind_lik_rem)
  
  # rescale estimates
  total_var = apply(lambda^2, 1, sum) + psi
  lambda = apply(lambda, 2, function(x) x / sqrt(total_var))
  psi = psi / total_var
  nu =  nu / sqrt(total_var)
  
  # rotate loadings
  lambda = varimax(lambda)
  
  output = list(nu = nu, 
                lambda = lambda$loadings, 
                psi = psi, 
                gamma = gamma,
                weights = weights,
                lik_rem = lik_rem,
                lik = lik)
  return(output)
}