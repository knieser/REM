#' REM Startup Message
#' @importFrom utils packageVersion


REMStartupMessage <- function()
{
  msg <- c(paste0(
    " This is REM version ",
    packageVersion("REM")),"\n REM is a BETA software.","\nFor documentation, questions, and issues please see github link.",
    "\nType 'citation(\"REM\")' for citing this R package in publications.")
  return(msg)
}

.onAttach <- function(lib, pkg)
{

  msg <- REMStartupMessage()
  if(!interactive())
    msg[1] <- paste("Package 'REM' version", packageVersion("REM"))
  packageStartupMessage(msg)
  invisible()
}
