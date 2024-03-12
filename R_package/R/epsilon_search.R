#' Optimal Epsilon Search
#' @param X numeric data frame containing all the variables needed in the factor analysis
#' @param delta hyperparameter between 0 and 1 that captures the researcher’s tolerance of incorrectly down-weighting data from the model
#' @param constraints p x k matrix of zeros and ones denoting the factors (rows) and observed variables (columns)
#' @param rotation factor rotation method (either 'oblimin' or 'varimax')
#' @param ueps upper epsilon
#' @param intl_mu initial mu (intercept).
#' @param intl_lambda initial lambda (loadings)
#' @param intl_psi initial psi (variance)
#' @param intl_phi initial phi (factor correlation matrix).
#' @param ctrREM provides control parameters (default: controlREM(steps = 25, tol = 1e-6, maxiter = 1e3))
#' @returns optimal epsilon.
#' @author Bryan Ortiz-Torres (bortiztorres@wisc.edu); Kenneth Nieser (nieser@wisc.edu)
#' @references Nieser, K. J., & Cochran, A. L. (2021). Addressing heterogeneous populations in latent variable settings through robust estimation. Psychological Methods.
#' @noRd

epsilon_search <- function(X, delta, constraints, rotation, ueps, intl_mu, intl_lambda, intl_psi, intl_phi, ctrREM = controlREM()){

  # number of factors
  k = ncol(intl_lambda) #initial loadings

  # number of steps
  steps = ctrREM$steps

  # lower bound for epsilon
  leps = 0

  # initialize
  eps_range = rep(0, steps)
  chk = rep(1, steps) #check the heuristic

  # binary search
  for (iter in 1:steps){

    # take midpoint between lower and upper epsilon
    eps_range[iter] = (leps + ueps)/2

    # get REM estimates
    REM_output <- RobustEMAlg(X, k, eps_range[iter], constraints, rotation, intl_mu, intl_lambda, intl_psi, intl_phi, ctrREM)
    mu2 = REM_output$mu
    lambda2 = REM_output$lambda
    psi2 = REM_output$psi
    phi2 = REM_output$phi

    # check heuristic
    chk_gamma <- ctrREM$chk_gamma
    chk[iter] = checkEpsFA(mu2, lambda2, psi2, phi2, eps_range[iter], chk_gamma)

    # update upper or lower epsilon
    if(chk[iter] < delta){
      leps = eps_range[iter]
    } else{
      ueps = eps_range[iter]}
  }

  opt_eps = (leps + ueps) / 2
  return(opt_eps)
}
