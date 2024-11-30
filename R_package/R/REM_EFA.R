#' Robust Estimation Maximization for Exploratory Factor Analysis
#' @description
#' This function uses the robust expectation maximization (REM) algorithm to estimate the parameters of an exploratory factor analysis model as suggested by Nieser & Cochran (2023).
#' @param X data to analyze; should be a data frame or matrix
#' @param k_range vector of the number of factors to consider
#' @param delta hyperparameter between 0 and 1 that captures the researcher’s tolerance of incorrectly down-weighting data from the model (default = 0.05)
#' @param rotation factor rotation method (default = 'oblimin'); 'varimax' is the only other available option at this time
#' @param ctrREM control parameters (default: (steps = 25, tol = 1e-6, maxiter = 1e3, min_weights = 1e-30, max_ueps =  0.3, chk_gamma = 0.9, n = 2e4))
#' @returns REM_EFA returns an object of class "REM". The function [summary()] is used to obtain estimated parameters from the model. An object of class "REM" in Exploratory Factor Analysis is a list of outputs with four different components for each number of factor: the matched call (call), estimates using traditional expectation maximization (EM_output), estimates using robust expectation maximization (REM_output), and a summary table (summary_table). The list contains the following components:
#'  \item{call}{match call}
#'  \item{model}{model frame}
#'  \item{k}{number of factors}
#'  \item{constraints}{p x k matrix of zeros and ones denoting the factors (rows) and observed variables (columns)}
#'  \item{epsilon}{hyperparameter on the likelihood scale}
#'  \item{AIC_rem}{Akaike information criterion based on REM estimates}
#'  \item{BIC_rem}{Bayesian information criterion based on REM estimates}
#'  \item{mu}{item intercepts}
#'  \item{lambda}{factor loadings}
#'  \item{psi}{unique variances of items}
#'  \item{phi}{factor covariance matrix}
#'  \item{gamma}{average weight}
#'  \item{weights}{estimated REM weights}
#'  \item{ind_lik}{likelihood value for each individual}
#'  \item{lik_rem}{joint log-likelihood evaluated at REM estimates}
#'  \item{lik}{joint log-likelihood evaluated at EM estimates}
#'  \item{mu.se}{standard errors of items intercepts}
#'  \item{lambda.se}{standard errors of factor loadings}
#'  \item{psi.se}{standard errors of unique variances of items}
#'  \item{gamma.se}{standard error of gamma}
#'  \item{summary_table}{summary of EM and REM estimates, SEs, Z statistics, p-values, and 95% confidence intervals}
#' @returns The summary function can be used to obtain estimated parameters from the optimal model based on the BIC from the EM and REM algorithms.
#' @author Bryan Ortiz-Torres (bortiztorres@wisc.edu); Kenneth Nieser (nieser@stanford.edu)
#' @references Nieser, K. J., & Cochran, A. L. (2023). Addressing heterogeneous populations in latent variable settings through robust estimation. Psychological methods, 28(1), 39.
#' @seealso [REM_CFA()], [summary.REMLA()] for more detailed summaries, [GPArotation::oblimin()] and [varimax()] for details on the rotation
#' @examples
#' \donttest{
#' # EFA of Holzinger-Swineford dataset
#' library(lavaan)
#' df <- HolzingerSwineford1939
#' data = df[,-c(1:6)]
#'
#' model_EFA = REM_EFA(X = data, k_range = 1:3)
#' summary(model_EFA)
#' }
#' @importFrom stats factanal quantile rnorm varimax na.omit cov2cor pnorm complete.cases
#' @importFrom GPArotation oblimin
#' @export

REM_EFA <- function(X, k_range, delta = 0.05, rotation = 'oblimin', ctrREM = controlREM()){

  if (!all(sapply(X, is.numeric))) stop("The dataset should be entirely numeric.")
  if(any(k_range <= 0) || !isTRUE(all(k_range == floor(k_range)))) stop("k_range should contain positive integers only.")
  if(delta < 0 || delta > 1) stop("delta values should be between 0 and 1.")

  n.missing <- sum(!complete.cases(X))
  X = na.omit(as.matrix(X))
  if (n.missing > 0) warning(paste0(n.missing, " rows with missing values were removed."))

  n = nrow(X)
  p = ncol(X)
  if (p < 3) stop("The dataset should include at least 3 variables.")

  cl <- match.call()
  REM_output = vector(mode = 'list', length = length(k_range))

  for (k in 1:length(k_range)){
    string <- paste('EFA with', k_range[k], 'factors', '\n')
    message(string)

    constraints = matrix(1, nrow = p, ncol = k_range[k])
    out = REM_estimates(X, k_range[k], delta, constraints, rotation, ctrREM)

    # summary table
    if (k_range[k] > 1){
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
                  paste0('phi', 1:(k_range[k]*(k_range[k]-1)/2)),
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
      method = c(rep(c('EM', 'REM'), each = 2*p + k_range[k]*(k_range[k]-1)/2 + sum(constraints==1)), 'REM'),
      parameter = par,
      estimates = theta,
      se = se,
      Z = theta / se,
      p.value = 2 *(1- pnorm(abs(theta/se))),
      ci.lower = theta - 1.96 * se,
      ci.upper = theta + 1.96 * se
    )
    out$summary_table <- summary_table
    REM_output[[k]] <- out
  }
  names(REM_output) <- paste0("nf", k_range)
  REM_output$call <- cl
  REM_output$delta <- delta
  class(REM_output) <- "REMLA"
  return(REM_output)
}
