#' REM Startup Message
#' @importFrom utils packageVersion
#' @noRd


REMLAStartupMessage <- function()
{
  msg <- c(paste0(
    " This is REMLA version ",
    packageVersion("REMLA")),"\nFor documentation, questions, and issues please see github link https://github.com/knieser/REM.",
    "\nType 'citation(\"REMLA\")' for citing this R package in publications.")
  return(msg)
}

.onAttach <- function(lib, pkg)
{

  msg <- REMLAStartupMessage()
  if(!interactive())
    msg[1] <- paste("Package 'REMLA' version", packageVersion("REMLA"))
  packageStartupMessage(msg)
  invisible()
}
