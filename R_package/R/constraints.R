#' CFA Constraints Matrix
#' @description
#' Function to create constraints matrix based on user inputs for Confirmatory Factor Analysis
#' @param X numeric data frame containing all the variables needed in the factor analysis
#' @param order original order of variables in the data frame
#' @returns matrix of zeros and ones indicating the relationship that wants to be modeled
#' @author Bryan Ortiz-Torres (bortiztorres@wisc.edu); Kenneth Nieser (nieser@wisc.edu)
#' @references Nieser, K. J., & Cochran, A. L. (2021). Addressing heterogeneous populations in latent variable settings through robust estimation. Psychological Methods.
#' @noRd


constraints <- function(X,order) {


  # add condition to make matrix in the same order of variables as the original dataset. Rearrange the order at the end.
  # add stop functions to make sure that model haves the same name of variables and add all of them in the model at least once.

  lines <- strsplit(X, "\n")[[1]]

  indicators <- list()
  latent_vars <- c()

  for (line in lines) {
    if (grepl("=~", line)) {
      parts <- strsplit(line, "=~")[[1]]
      latent_var <- trimws(parts[1])
      indicators <- c(indicators, strsplit(trimws(parts[2]), "\\+")[[1]]) #change this part to allow for not spaces between plus signs
      indicators <- trimws(indicators)
      latent_vars <- c(latent_vars, latent_var)
    }
  }

  unique_indicators <- unique(indicators)

  matrix_data <- matrix(0, nrow = length(unique_indicators), ncol = length(latent_vars))
  rownames(matrix_data) <- unique_indicators
  colnames(matrix_data) <- latent_vars

  for (line in lines) {
    if (grepl("=~", line)) {
      parts <- strsplit(line, "=~")[[1]]
      latent_var <- trimws(parts[1])
      indicators <- strsplit(trimws(parts[2]), "\\+")[[1]]
      indicators <- trimws(indicators)

      for (indicator in indicators) {
        matrix_data[indicator, latent_var] <- 1
      }
    }
  }
  matrix_data2 = matrix_data[order,]
  return(matrix_data)
}
