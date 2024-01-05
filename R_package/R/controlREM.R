#' Control parameters for REM package
#' @param steps number of steps in binary search for optimal epsilon value (default = 25)
#' @param tol tolerance parameter to check for convergence of EM and REM algorithm (default = 1e-6)
#' @param maxiter maximum number iterations of EM and REM algorithm (default = 1e3)
#' @param min_weights lower bound for the individual weights estimated by REM (default = 1e-30)
#' @param max_ueps percentile of the distribution of likelihood values to use as the maximum epsilon value to consider
#' @param chk_gamma gamma value used when searching for epsilon
#' @param n sample size of simulated data used when checking heuristic criterion in the epsilon search
#' @return control parameters used in the REM package (steps, tol, maxiter, min_weights, ueps, n).
#' @author Bryan Ortiz-Torres (bortiztorres@wisc.edu); Kenneth Nieser (nieser@stanford.edu)
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
