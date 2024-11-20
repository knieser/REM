#' Standard errors
#' @param data numeric data frame containing all the variables needed in the factor analysis
#' @param constraints p x k matrix of zeros and ones denoting the factors (rows) and observed variables (columns)
#' @param opt.eps optimal epsilon value from the [epsilon_search()]
#' @param EM_output list of Expectation Maximization output from [REM_estimates()]
#' @param REM_output list of Robust Expectation Maximization output from [REM_estimates()]
#' @returns List of estimates:
#'  \item{results}{data frame with list of parameter estimates, SEs, Z test statistic, p-value, and 95% confidence intervals}
#' @author Bryan Ortiz-Torres (bortiztorres@wisc.edu); Kenneth Nieser (nieser@stanford.edu)
#' @references Nieser, K. J., & Cochran, A. L. (2023). Addressing heterogeneous populations in latent variable settings through robust estimation. Psychological methods, 28(1), 39.
#' @importFrom geex m_estimate vcov
#' @noRd

calcSE <- function(data, constraints, opt.eps, EM_output, REM_output){
  # rescale data so variance = 1
  sd = apply(data, 2, sd)
  data.standard = as.data.frame(t(apply(data, 1, function(y) y/sd)))

  n = nrow(data.standard)
  p = ncol(data.standard)
  k = ncol(constraints)
  lambda.parms = sum(constraints==1)

  theta.EM = as.numeric(c(
    EM_output$mu,
    EM_output$lambda[constraints==1],
    EM_output$phi[lower.tri(EM_output$phi)],
    EM_output$psi
    ))

  theta.REM = as.numeric(c(
    REM_output$mu,
    REM_output$lambda[constraints==1],
    REM_output$phi[lower.tri(REM_output$phi)],
    REM_output$psi,
    REM_output$gamma
    ))

  theta = c(theta.EM, theta.REM)

  weights = REM_output$weights

  # estimating function
  # if weights are supplied, then it uses REM estimating equations; otherwise uses EM estimating equations
  estFUN <- function(data, constraints, REM = 0, weights = NA){
    mu = apply(data, 2, mean)
    n = nrow(data)
    p = ncol(data)
    k = ncol(constraints)
    lambda = matrix(0, nrow = p, ncol = k)
    lambda.idx = which(constraints==1)
    lambda.parms = length(lambda.idx)
    phi = matrix(0, nrow = k, ncol = k)
    d.lambda = matrix(0, nrow = n, ncol = lambda.parms)
    d.phi = matrix(0, nrow = n, ncol = k*(k-1)/2)
    d.psi = matrix(0, nrow = n, ncol = p)
    d.gamma = vector(length = n)

    if (REM == 1){
    function(theta){
      mu = theta[1:p]
      lambda.est = theta[(p+1):(p+lambda.parms)]
      lambda[lambda.idx] <- lambda.est
      phi.tri = theta[(p+lambda.parms+1):(p+lambda.parms+k*(k-1)/2)]
      diag(phi) <- 1
      phi[lower.tri(phi)] <- phi.tri
      phi[upper.tri(phi)] <- phi.tri
      psi = theta[(p+lambda.parms+k*(k-1)/2 + 1):(p+lambda.parms+k*(k-1)/2 + p)]
      gamma = theta[p+lambda.parms+k*(k-1)/2+p+1]

      # estimating equations
      Z = apply(data, 1, function(y) y - mu)
      sigma = lambda %*% phi %*% t(lambda) + diag(psi)
      inv.sigma = solve(sigma)
      d.mu = t(inv.sigma %*% Z)
      for (i in 1:n){
        D = weights[i] * inv.sigma %*% (diag(p) - (Z[,i] %*% t(Z[,i]) %*% inv.sigma))
        d.lambda[i,] = -weights[i] * (D %*% lambda %*% phi)[lambda.idx]
        d.phi[i,] = -weights[i] * (t(lambda) %*% D %*% lambda)[lower.tri(phi)]
        d.psi[i,] = -0.5 * weights[i] * diag(D)
        d.gamma[i] = (weights[i] - gamma)/(gamma)/(1-gamma)
      }

      # score function for each person and parameter
      c(d.mu, d.lambda, d.phi, d.psi, d.gamma)
    }
    } else{
      function(theta){
        mu = theta[1:p]
        lambda.est = theta[(p+1):(p+lambda.parms)]
        lambda[lambda.idx] <- lambda.est
        phi.tri = theta[(p+lambda.parms+1):(p+lambda.parms+k*(k-1)/2)]
        diag(phi) <- 1
        phi[lower.tri(phi)] <- phi.tri
        phi[upper.tri(phi)] <- phi.tri
        psi = theta[(p+lambda.parms+k*(k-1)/2 + 1):(p+lambda.parms+k*(k-1)/2 + p)]

        # estimating equations
        Z = apply(data, 1, function(y) y - mu)
        sigma = lambda %*% phi %*% t(lambda) + diag(psi)
        inv.sigma = solve(sigma)
        d.mu = t(inv.sigma %*% Z)
        for (i in 1:n){
          D = inv.sigma %*% (diag(p) - (Z[,i] %*% t(Z[,i]) %*% inv.sigma))
          d.lambda[i,] = -(D %*% lambda %*% phi)[lambda.idx]
          d.phi[i,] = -(t(lambda) %*% D %*% lambda)[lower.tri(phi)]
          d.psi[i,] = -0.5 * diag(D)
        }

        # score function for each person and parameter
        c(d.mu, d.lambda, d.phi, d.psi)
      }
    }
  }

  # use geex package to estimate standard errors
  m.results.EM <- geex::m_estimate(
    estFUN = estFUN,
    data = data.standard,
    roots = theta.EM,
    compute_roots = FALSE,
    outer_args = list(constraints = constraints))

  m.results.REM <- geex::m_estimate(
    estFUN = estFUN,
    data = data.standard,
    roots = theta.REM,
    compute_roots = FALSE,
    outer_args = list(constraints = constraints, REM = 1, weights = weights))

  se.EM = sqrt(diag(geex::vcov(m.results.EM)))
  se.REM = sqrt(diag(geex::vcov(m.results.REM)))
  se = c(se.EM, se.REM)
  name.prefix = c(paste0(rep(c('EM_', 'REM_'), each = 2*p + lambda.parms + k*(k-1)/2)), 'REM_')
  if (k > 1) {names(se) <-paste0(name.prefix, c(rep(
    c(paste0('mu', 1:p),
      paste0('lambda', which(constraints==1)),
      paste0('phi', 1:(k*(k-1)/2)),
      paste0('psi', 1:p)
      ), 2), "gamma"))
  } else{names(se) <-paste0(name.prefix, c(rep(
    c(paste0('mu', 1:p),
      paste0('lambda', which(constraints==1)),
      paste0('psi', 1:p)
    ), 2), "gamma"))}

  return(se)
}
