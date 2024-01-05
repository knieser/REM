#' Control parameters for REM package
#' @param steps maximum number of epsilons considered in the analysis (default = 25).
#' @param tol tolerance parameter to check for convergence (default = 1e-6).
#' @param maxiter maximum iterations allowed to calculate the estimates (default = 1e3)
#' @param min_weights Lower bound for the individual weights estimated by REM
#' @param max_ueps The largest epsilon value to check when searching for the optimal epsilon.
#' @param chk_gamma gamma value used when searching for epsilon
#' @param n The simulated n in checkEpsFA
#' @return control parameters used in the REM package (steps, tol, maxiter, min_weights,ueps,n).
#' @author Bryan Ortiz-Torres (bortiztorres@wisc.edu); Kenneth Nieser (nieser@wisc.edu)
#' @references Nieser, K. J., & Cochran, A. L. (2021). Addressing heterogeneous populations in latent variable settings through robust estimation. Psychological Methods.
#' @seealso [REM_EFA()], [REM_CFA()]
#' @export

controlREM <- function(steps = 25, tol = 1e-6, maxiter = 1e3, min_weights = 1e-30, max_ueps =  0.3, chk_gamma = 0.9, n = 2e4)
{
  output = list()

  output$steps <- steps
  output$ tol <- tol
  output$maxiter <- maxiter
  output$chk_gamma <-  chk_gamma
  output$min_weights <- min_weights
  output$max_ueps <- max_ueps
  output$n <- n

  return(output)
}
