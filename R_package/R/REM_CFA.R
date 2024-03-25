#' Robust Estimation Maximization Estimates for Confirmatory Factor Analysis
#' @description
#' This function uses the robust expectation maximization (REM) algorithm to estimate the parameters of a confirmatory factor analysis model as suggested by Nieser & Cochran (2021).
#' @param X data to analyze; should be a data frame or matrix
#' @param delta hyperparameter between 0 and 1 that captures the researcher’s tolerance of incorrectly down-weighting data from the model (default = 0.05).
#' @param model string variable that contains each structural equation in a new line where equalities are denoted by the symbol "~".
#' @param ctrREM control parameters (default: (steps = 25, tol = 1e-6, maxiter = 1e3, min_weights = 1e-30, max_ueps =  0.3, chk_gamma = 0.9, n = 2e4))
#' @returns REM_CFA returns an object of class "REM". The function [summary()] is used to obtain estimated parameters from the model. An object of class "REM" in Confirmatory Factor Analysis is a list of outputs with four different components: the matched call (call), estimates using traditional expectation maximization (EM_output), estimates using robust expectation maximization (REM_output), and a summary table (summary_table). The list contains the following components:
#'  \item{call}{match call}
#'  \item{model}{model frame}
#'  \item{delta}{hyperparameter between 0 and 1 that captures the researcher’s tolerance of incorrectly down-weighting data from the model}
#'  \item{k}{number of factors}
#'  \item{constraints}{p x k matrix of zeros and ones denoting the factors (rows) and observed variables (columns)}
#'  \item{epsilon}{hyperparameter on the likelihood scale}
#'  \item{AIC_rem}{Akaike Information Criterion}
#'  \item{BIC_rem}{Bayesian Information Criterion}
#'  \item{mu}{item intercepts}
#'  \item{lambda}{factor loadings}
#'  \item{psi}{unique variances of items}
#'  \item{gamma}{average weights}
#'  \item{weights}{estimated REM weights}
#'  \item{ind_lik}{likelihood value for each individual}
#'  \item{lik_rem}{joint log-likelihood evaluated at REM estimates}
#'  \item{lik}{joint log-likelihood evaluated at EM estimates}
#'  \item{summary_table}{summary of EM and REM estimates, SEs, Z statistics, p-values, and 95% confidence intervals}
#' @author Bryan Ortiz-Torres (bortiztorres@wisc.edu); Kenneth Nieser (nieser@stanford.edu)
#' @references Nieser, K. J., & Cochran, A. L. (2021). Addressing heterogeneous populations in latent variable settings through robust estimation. Psychological Methods.
#' @seealso [REM_EFA()], [summary.REMLA()]
#' @examples
#' \donttest{
#' # Creating latent model
#' library(lavaan)
#' library(GPArotation)
#' df <- HolzingerSwineford1939
#' data = df[,-c(1:6)]
#'
#' model <- "Visual  =~  x1 + x2 + x3
#'          Textual =~  x4 + x5 + x6
#'          Speed   =~  x7 + x8 + x9"
#'
#' # Modeling Confirmatory Factor Analysis
#' model_CFA = REM_CFA(X = data, delta = 0.05, model = model)
#' summary(model_CFA)
#' }
#' @importFrom stats factanal quantile rnorm varimax na.omit cov2cor pnorm
#' @importFrom GPArotation oblimin
#' @export
REM_CFA <- function(X, delta = 0.05, model = NA, ctrREM = controlREM()){

  if (any(is.na(X) == TRUE)) warning("rows with missing values were removed")
  X = na.omit(as.matrix(X))
  n = nrow(X)
  p = ncol(X)
  order <- colnames(X)

  constraints <- constraints(model, order)
  rownames(constraints) <- NULL
  colnames(constraints) <- NULL
  k = ncol(constraints)

  # error checking for constraints matrix
  try(if(nrow(constraints) != p) stop(paste0("constraints should have p = ", p, " rows")))
  try(if(ncol(constraints) != k) stop(paste0("constraints should have k = ", k, " columns")))
  try(if(any(!(constraints %in% c(0,1)))) stop(paste0("constraints should only contain 0s and 1s")))

  string <- paste('CFA with', k, 'factors', '\n')
  message(string)
  out <- REM_estimates(X, k, delta, constraints, rotation = 0, ctrREM)

  # summary table
  if (k > 1){
    theta = c(out$EM_output$mu,
              out$EM_output$lambda[constraints==1],
              out$EM_output$phi[lower.tri(out$EM_output$phi)],
              out$EM_output$psi,
              out$REM_output$mu,
              out$REM_output$lambda[constraints==1],
              out$REM_output$phi[lower.tri(out$REM_output$phi)],
              out$REM_output$psi,
              out$REM_output$gamma)
    se = c(out$EM_output$mu.se,
           out$EM_output$lambda.se,
           out$EM_output$phi.se,
           out$EM_output$psi.se,
           out$REM_output$mu.se,
           out$REM_output$lambda.se,
           out$REM_output$phi.se,
           out$REM_output$psi.se,
           out$REM_output$gamma.se)
    par = c(rep(c(paste0('mu', 1:p),
                  paste0('lambda', which(constraints==1)),
                  paste0('phi', 1:(k*(k-1)/2)),
                  paste0('psi', 1:p)), 2), "gamma")
  } else{
    theta = c(out$EM_output$mu,
              out$EM_output$lambda[constraints==1],
              out$EM_output$psi,
              out$REM_output$mu,
              out$REM_output$lambda[constraints==1],
              out$REM_output$psi,
              out$REM_output$gamma)
    se = c(out$EM_output$mu.se,
           out$EM_output$lambda.se,
           out$EM_output$psi.se,
           out$REM_output$mu.se,
           out$REM_output$lambda.se,
           out$REM_output$psi.se,
           out$REM_output$gamma.se)
    par = c(rep(c(paste0('mu', 1:p),
                  paste0('lambda', which(constraints==1)),
                  paste0('psi', 1:p)), 2), "gamma")
  }

  summary_table <- data.frame(
    method = c(rep(c('EM', 'REM'), each = 2*p + k*(k-1)/2 + sum(constraints==1)), 'REM'),
    parameter = par,
    estimates = theta,
    se = se,
    Z = theta / se,
    p.value = 2 *(1- pnorm(abs(theta/se))),
    ci.lower = theta - 1.96 * se,
    ci.upper = theta + 1.96 * se
  )

  output = list(
    call = match.call(),
    model = model,
    delta = delta,
    k = out$k,
    constraints = out$constraints,
    epsilon = out$epsilon,
    AIC_rem = out$AIC_rem,
    BIC_rem = out$BIC_rem,
    EM_output = out$EM_output,
    REM_output = out$REM_output,
    summary_table = summary_table
)

  class(output) <- "REMLA"
  return(output)
}
