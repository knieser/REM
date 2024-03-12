#' Expectation Maximization Algorithm
#' @param X numeric data frame containing all the variables needed in the factor analysis
#' @param k number of factors considered in the analysis
#' @param constraints p x k matrix of zeros and ones denoting the factors (rows) and observed variables (columns)
#' @param rotation factor rotation method (either 'oblimin' or 'varimax')
#' @param ctrREM provides control parameters (default: controlREM(steps = 25, tol = 1e-6, maxiter = 1e3))
#' @returns List of outputs that contains the following components:
#'  \item{mu}{intercept}
#'  \item{lambda}{loadings}
#'  \item{psi}{variance}
#'  \iten(phi){factor covariance matrix}
#'  \item{ind_lik}{likelihood value for each individual}
#' @author Bryan Ortiz-Torres (bortiztorres@wisc.edu); Kenneth Nieser (nieser@wisc.edu)
#' @references Nieser, K. J., & Cochran, A. L. (2021). Addressing heterogeneous populations in latent variable settings through robust estimation. Psychological Methods.
#' @importFrom stats factanal quantile rnorm varimax na.omit cov2cor pnorm
#' @importFrom GPArotation oblimin
#' @noRd

EMAlg <- function(X, k, constraints, rotation, ctrREM = controlREM()){

  # data dimensions
  n = nrow(X)
  p = ncol(X)

  # sample mean and covariance matrix
  mu = apply(X, 2, mean)
  X_centered = t(apply(X, 1, function(y) y - mu))
  Cxx = 1/n * t(X_centered) %*% X_centered

  if(k == 1)
  {
    intl_guess <- factanal(X, k) #Maximum likelihood factor analysis  (r function)
    inv_rotmat <- 1
  } else {
    intl_guess <- factanal(X, k, rotation = 'oblimin') #Maximum likelihood factor analysis  (r function)
    inv_rotmat <- solve(intl_guess$rotmat)
  }

  lambda <- matrix(intl_guess$loadings, ncol = k)
  lambda[constraints==0] <- 0
  psi <- intl_guess$uniquenesses #describes variances of the error terms (unique factors $\epsilon in the paper$)
  phi <- inv_rotmat %*% t(inv_rotmat)

  # tolerance parameters
  tol = ctrREM$tol
  maxiter = ctrREM$maxiter

  # constraint patterns to be processed separately
  patterns = unique(constraints)
  pattern.id = 1:nrow(patterns)
  for (j in 1:nrow(patterns)){
    idx = which(apply(constraints, 1, function(x) return(all(x == c(patterns[j,])))))
    pattern.id[idx] = j
  }

  # EM algorithm
  obj = rep(-Inf, maxiter)
  for (iter in 2:maxiter){

    # calculate log-likelihoods
    sigma = lambda %*% phi %*% t(lambda) + diag(psi)
    V = chol(sigma)
    inv_V = solve(V)
    logdetV = log(det(V))
    Z = t(apply(X, 1, function(y) y - mu)) %*% inv_V
    ind_lik = -(1/2)*(p*log(2*pi) + 2*logdetV + apply(Z, 1, function(y) t(y) %*% y))

    # calculate value for objective function
    obj[iter] = sum(ind_lik)

    # check for convergence
    if(abs(obj[iter] - obj[iter-1]) < tol*abs(obj[iter-1])){break}

    # E step
    beta = phi %*% (t(lambda) %*% inv_V) %*% t(inv_V)
    CxxB = Cxx %*% t(beta)
    Ezz = (diag(k) - beta %*% lambda) %*% phi + beta %*% CxxB

    # M step
    # estimate the subsets of lambda and psi for each constraint pattern
    for (i in 1:nrow(patterns)){
      item_idx <- pattern.id==i
      factor_idx <- patterns[i,]!=0
      Cxx.ii = diag(Cxx)[item_idx]
      CxxB.i = matrix(CxxB[item_idx, factor_idx], nrow = length(Cxx.ii))
      Ezz.i = Ezz[factor_idx, factor_idx]
      chol_Ezz.i = chol(Ezz.i)
      inv_chol_Ezz.i = solve(chol_Ezz.i)
      lambda[item_idx, factor_idx] = (CxxB.i %*% inv_chol_Ezz.i) %*% t(inv_chol_Ezz.i)
      psi[item_idx] = Cxx.ii - diag(((CxxB.i %*% inv_chol_Ezz.i) %*% t(inv_chol_Ezz.i)) %*% t(CxxB.i))
      }
      phi = cov2cor(phi - beta %*% lambda %*% phi + beta %*% CxxB)

      if(iter==maxiter){message(
        paste('maxiter reached; consider increasing maxiter or increasing convergence tolerance'))}
  }

 # rescale estimates
  sigma = lambda %*% phi %*% t(lambda) + diag(psi)
  total_var = diag(sigma)
  lambda = apply(lambda, 2, function(x) x / sqrt(total_var))
  psi = psi / total_var
  mu =  mu / sqrt(total_var)

  # rotate loadings
  if (k > 1){
    if (rotation == "oblimin"){lambda = matrix(GPArotation::oblimin(lambda)$loadings, ncol = k)}
    if (rotation == "varimax"){lambda = matrix(varimax(lambda)$loadings, ncol = k)}
  }

  output = list(mu = mu, lambda = lambda, psi = psi, phi = phi, ind_lik = ind_lik)
  return(output)
}
