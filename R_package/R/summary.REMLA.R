#' Summary for Robust Estimation Maximization
#' @description
#' Summary method for class "REMLA".
#'
#' @param object an object of class "REMLA", usually a result of a call to [REM_EFA].
#' @param ... further arguments passed to or from other methods.
#' @return The summary.REM function returns estimated parameters from the optimal model based on the BIC from the EM and REM algorithms.
#' @return Output include:
#'  \item{optimal}{optimal number of factors based on BIC}
#'  \item{mu}{intercept}
#'  \item{lambda}{loadings}
#'  \item{psi}{variance}
#'  \item{indk_lik}{likelihood value for each individual}
#'  \item{epsilon}{hyperparameter on the likelihood scale}
#'  \item{diff}{differences between EM and REM}
#' @author Bryan Ortiz-Torres (bortiztorres@wisc.edu); Kenneth Nieser (nieser@stanford.edu)
#' @references Nieser, K. J., & Cochran, A. L. (2023). Addressing heterogeneous populations in latent variable settings through robust estimation. Psychological methods, 28(1), 39.
#' @seealso [REM_EFA()], [REM_CFA()], [summary()].
#' @importFrom stats factanal quantile rnorm varimax
#' @export


summary.REMLA = function(object,...) {

  x = object

  if (x[["call"]][1]=="REM_CFA()") {
    idx = 1
    k = x$k
    optimal = x$k
    epsilon = x$epsilon
    gamma = x$REM_output$gamma
    EM <- cbind(x$EM_output$mu,x$EM_output$lambda,x$EM_output$psi)
    EM_phi <- x$EM_output$phi
    REM <- cbind(x$REM_output$mu,x$REM_output$lambda,x$REM_output$psi)
    REM_phi <- x$REM_output$phi
  } else {
    j = length(x) - 2
    k <- vector(length = j)
    BIC <- vector(length = j)

    for (i in 1:j) {
      k[i] = x[[i]]$k
      BIC[i] = x[[i]]$BIC_rem
    }

    idx = which.min(BIC)
    optimal = x[[idx]]$k
    epsilon = x[[idx]]$epsilon
    gamma = x[[idx]]$REM_output$gamma
    EM <- cbind(x[[idx]]$EM_output$mu,x[[idx]]$EM_output$lambda,x[[idx]]$EM_output$psi)
    EM_phi <- x[[idx]]$EM_output$phi
    REM <- cbind(x[[idx]]$REM_output$mu,x[[idx]]$REM_output$lambda,x[[idx]]$REM_output$psi)
    REM_phi <- x[[idx]]$REM_output$phi
  }

  colnames(REM) <- c("Intercept",paste("Factor",1:k[idx]),"Res Var")
  colnames(EM) <- c("Intercept",paste("Factor",1:k[idx]),"Res Var")
  colnames(EM_phi) <- paste("Factor",1:k[idx])
  colnames(REM_phi) <- paste("Factor",1:k[idx])
  rownames(EM_phi) <- paste("Factor",1:k[idx])
  rownames(REM_phi) <- paste("Factor",1:k[idx])

  if (x[["call"]][1]=="REM_CFA()") {
    cat("-----Confirmatory factor analysis using REM-----", "\n")
  } else{
      cat("-----Exploratory factor analysis using REM-----", "\n")
    }
  cat("Call:", deparse(x$call) , "\n\n")
  if (x[["call"]][1]=="REM_EFA()") {
    cat("According to BIC evaluated at REM estimates,\n the optimal number of factors is:", optimal, "\n\n")
  } else{
    cat("Model:\t", x$model, "\n\n")
  }
  cat("REM estimates w/ delta =", x$delta, "( epsilon =", epsilon, ")\n", "--------------------------", "\n")
  cat("Est. gamma =", round(gamma, 3), "\n\n")
  print(round(REM,3))
  cat("\n")
  cat("Factor correlations", "\n")
  print(round(REM_phi, 3))
  cat("\n")
  cat("EM estimates", "\n", "--------------------------", "\n")
  print(round(EM, 3))
  cat("\n")
  cat("Factor correlations", "\n")
  print(round(EM_phi,3))
  cat("\n")
  cat("Difference (EM - REM)", "\n", "--------------------------", "\n")
  print(round(EM - REM,3))
  cat("\n")
  cat("Factor correlations", "\n")
  print(round(EM_phi - REM_phi,3))

  ans <- x["call"]
  ans$optimal <-  optimal
  ans$EM <- EM
  ans$REM <- REM
  ans$epsilon <- epsilon
  ans$gamma <- gamma
  ans$diff <- as.data.frame(ans$EM-ans$REM)
  class(ans) <- "summary.REMLA"

  #return(ans)
}
