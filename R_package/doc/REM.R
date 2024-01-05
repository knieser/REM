## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width=6, fig.height=4
)

## ----setup--------------------------------------------------------------------
library(REM)

## -----------------------------------------------------------------------------
library(lavaan)
library(dplyr)

df <- HolzingerSwineford1939
head(df)

## -----------------------------------------------------------------------------
data=df[,7:15]

## -----------------------------------------------------------------------------
model_EFA = REM_EFA( X = data, k_range = 2:4, delta = 0.05)

## -----------------------------------------------------------------------------
summary(model_EFA)

## -----------------------------------------------------------------------------

hist(model_EFA[[2]]$REM_output$weights,main="REM Weight Distribution",xlab = "Weights", ylab = "Frequency")

## -----------------------------------------------------------------------------
library(dplyr)
library(ggplot2)

 x <- model_EFA
  k <- vector()
  j = length(x) - 1
  BIC <- vector(length = j)

  for (i in 1:j) {
    k[i] = x[[i]]$k
    BIC[i] = x[[i]]$BIC_rem
  }

  idx = which.min(BIC)
  optimal = x[[idx]]$k
  epsilon = x[[idx]]$epsilon
  gamma = x[[idx]]$REM_output$gamma
  EM <- cbind(x[[idx]]$EM_output$nu,x[[idx]]$EM_output$lambda,x[[idx]]$EM_output$psi)
  REM <- cbind(x[[idx]]$REM_output$nu,x[[idx]]$REM_output$lambda,x[[idx]]$REM_output$psi)

  colnames(REM) <- c("Intercept",paste("Factor",1:k[idx]),"Res Var")
  colnames(EM) <- c("Intercept",paste("Factor",1:k[idx]),"Res Var")


  ans <- x["call"]
  ans$optimal <-  optimal
  ans$EM <- EM
  ans$REM <- REM
  ans$epsilon <- epsilon
  ans$gamma <- gamma
  class(ans) <- "summary.REM"
  
diff = as.data.frame(ans$EM-ans$REM)
diff_plot <- ggplot(diff, aes(1:nrow(diff)),) +  
    geom_line(aes(y = diff$Intercept), color = "red") +
    geom_line(aes(y = diff$`Factor 1`), color = "orange") +
    geom_line(aes(y = diff$`Factor 2`), color = "green") +
    geom_line(aes(y = diff$`Factor 3`), color = "blue") +
    geom_line(aes(y = diff$`Res Var`), color = "purple") +
    geom_hline(aes(yintercept=0),colour="black", linetype="dashed") +
    xlab("Observations") +
    ylab("Differences")

diff_plot

## -----------------------------------------------------------------------------
set.seed(1)
cons = matrix(nrow = ncol(data),ncol = 3)
cons[,1] = 1
cons[,2] = as.numeric(rbinom(ncol(data),1,0.5))
cons[,3] = 1
cons

model_CFA = REM_CFA(X = data,k=3,delta = 0.05,constraints = cons)
model_CFA

