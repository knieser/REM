checkEps <- function(mix, mu, sigma, epsilon, gamma){
  
  # simulate data
  n = 2e4
  p = nrow(mu)
  k = ncol(mu)
  
  # random draws from uniform
  u = runif(n)
  
  # get cumulative sums to check against
  acc_mix = cumsum(mix)
  
  # calculate the randomly generated number of people in each group
  acc_m = sapply(acc_mix, function(x) sum(u < x))
  m = c(acc_m[1], diff(acc_m))
  
  # sample from MVN for each latent group
  X = matrix(data = 0, nrow = 1, ncol = p)
  for (j in 1:k){
    X0 = mvrnorm(m[j], mu = mu[,j], Sigma = sigma[,,j])
    X = rbind(X, X0)
  }
  X = X[-1,]
  
  # calculate log-likelihoods
  cond_lik = matrix(data = 0, nrow = n, ncol = k)
  for (j in 1:k){
    V = chol(sigma[,,j])
    inv_V = solve(V)
    logdetV = log(det(V))
    Z = sweep(X, 2, mu[,j], "-") %*% inv_V
    cond_lik[,j] = -(1/2)*(p*log(2*pi) + 2*logdetV + apply(Z, 1, function(y) t(y) %*% y))
  }
  omega_num = sweep(exp(cond_lik), 2, mix, "*")
  ind_lik = log(apply(omega_num, 1, sum))
  
  # calculate modified log-likelihoods
  ind_lik_rem = log(gamma * exp(ind_lik) + (1-gamma) * epsilon)
  
  # compute weights
  weights_num = log(gamma) + ind_lik
  weights_denom = ind_lik_rem
  log_weights = weights_num - weights_denom
  weights = exp(log_weights)
  
  # calculate heuristic check value
  chk = 1 - mean(weights)
  
  return(chk)
}