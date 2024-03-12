#' Robust Estimation Maximization Estimates
#' @param X numeric data frame containing all the variables needed in the factor analysis
#' @param k number of factors considered in the model
#' @param delta hyperparameter between 0 and 1 that captures the researcher’s tolerance of incorrectly down-weighting data from the model (default = 0.05).
#' @param constraints p x k matrix of zeros and ones denoting the factors (rows) and observed variables (columns)
#' @param rotation factor rotation method (either 'oblimin' or 'varimax')
#' @param ctrREM provides control parameters (default: controlREM(steps = 25, tol = 1e-6, maxiter = 1e3))
#' @returns List of estimates from the REM algorithm:
#'  \item{k}{number of factors}
#'  \item{constraints}{matrix of zeros and ones denoting the factors (rows) and observed variables (columns)}
#'  \item{EM_output}{list of outputs. See [EMAlg()] for more details.}
#'  \item{REM_output}{list of outputs. See [REMAlg()] for more details.}
#'  \item{epsilon}{hyperparameter on the likelihood scale}
#'  \item{AIC_rem}{Akaike Information Criterion}
#'  \item{BIC_rem}{Bayesian Information Criterion}
#' @author Bryan Ortiz-Torres (bortiztorres@wisc.edu); Kenneth Nieser (nieser@wisc.edu)
#' @references Nieser, K. J., & Cochran, A. L. (2021). Addressing heterogeneous populations in latent variable settings through robust estimation. Psychological Methods.
#' @noRd

REM_estimates <- function(X, k, delta, constraints, rotation, ctrREM = controlREM()){

  n = nrow(X)
  p = ncol(X)

  # EM estimates
  message(' getting EM estimates...')
  EM_output <- EMAlg(X, k, constraints, rotation, ctrREM)
  mu1 = EM_output$mu
  lambda1 = EM_output$lambda
  psi1 = EM_output$psi
  phi1 = EM_output$phi
  ind_lik = EM_output$ind_lik
  message('done', '\n')

  # epsilon search
  message(' searching for optimal epsilon')
  max_ueps <- ctrREM$max_ueps
  ueps = quantile(exp(ind_lik), max_ueps) #the epsilon value is below the 0.30 percentile of the likelihood
  opt_eps = epsilon_search(X, delta, constraints, rotation, ueps, mu1, lambda1, psi1, phi1, ctrREM)
  message('done', '\n')

  # REM algorithm
  message(' getting REM estimates...')
  REM_output = RobustEMAlg(X, k, opt_eps, constraints, rotation, mu1, lambda1, psi1, phi1, ctrREM)
  message('done', '\n')

  # standard errors
  message(' calculating standard errors...')
  se <- calcSE(X, constraints, opt_eps, EM_output, REM_output)
  EM_output$mu.se = se[grep('^EM_mu', names(se))]
  EM_output$lambda.se = se[grep('^EM_lambda', names(se))]
  if (k > 1){EM_output$phi.se = se[grep('^EM_phi', names(se))]}
  EM_output$psi.se = se[grep('^EM_psi', names(se))]
  REM_output$mu.se = se[grep('^REM_mu', names(se))]
  REM_output$lambda.se = se[grep('^REM_lambda', names(se))]
  if (k > 1){REM_output$phi.se = se[grep('^REM_phi', names(se))]}
  REM_output$psi.se = se[grep('^REM_psi', names(se))]
  REM_output$gamma.se = se[grep('^REM_gamma', names(se))]
  message('done', '\n')

  # calculate AIC and BIC
  lambda.parms = sum(constraints==1)
  if(rotation=='varimax'){parms = p + lambda.parms + p} else{parms = p + lambda.parms + p + k*(k-1)/2}
  lik_rem = REM_output$lik_rem
  AIC_rem = -2*lik_rem + 2*parms
  BIC_rem = -2*lik_rem + log(n)*parms

  output = list(
    k = k,
    constraints = constraints,
    EM_output = EM_output,
    REM_output = REM_output,
    epsilon = opt_eps,
    AIC_rem = AIC_rem,
    BIC_rem = BIC_rem)
  return(output)
}
