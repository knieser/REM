EMAlg <- function(X, k){
  
  # data dimensions
  n = nrow(X)
  p = ncol(X)
  
  # sample mean and covariance matrix
  nu = apply(X, 2, mean)
  X_centered = sweep(X, 2, nu, "-")
  Cxx = 1/n * t(X_centered) %*% X_centered
  
  # initial guess of lambda and psi
  intl_guess <- factanal(X, k)
  lambda <- intl_guess$loadings
  psi <- intl_guess$uniquenesses
  
  # tolerance parameters
  tol = 1e-6
  maxiter = 1e3
  
  # EM algorithm
  obj = rep(-Inf, maxiter)
  for (iter in 2:maxiter){
    
    # calculate log-likelihoods
    sigma = lambda %*% t(lambda) + diag(psi)
    V = chol(sigma)
    inv_V = solve(V)
    logdetV = log(det(V))
    Z = sweep(X, 2, nu, "-") %*% inv_V
    ind_lik = -(1/2)*(p*log(2*pi) + 2*logdetV + apply(Z, 1, function(y) t(y) %*% y))
    
    # calculate value for objective function
    obj[iter] = sum(ind_lik)
    
    # check for convergence
    if(abs(obj[iter] - obj[iter-1]) < tol*abs(obj[iter-1])){break}
    
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
  
  # rescale estimates
  total_var = apply(lambda^2, 1, sum) + psi
  lambda = apply(lambda, 2, function(x) x / sqrt(total_var))
  psi = psi / total_var
  nu =  nu / sqrt(total_var)
  
  # rotate loadings
  lambda = varimax(lambda)
  
  output = list(nu = nu, lambda = lambda$loadings, psi = psi, ind_lik = ind_lik)
  return(output)
}