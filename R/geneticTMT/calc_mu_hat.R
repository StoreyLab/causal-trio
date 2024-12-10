#' Population Mean Estimate
#'
#' This function computes the unbiased estimate for the population mean.
#'
#' @param Y The n-length vector of children' phenotype.
#' @param W0 A vector of indicators for assigning to group 0.
#' @param W1 A vector of indicators for assigning to group 1. 
#'
#' @return A numeric value of the population mean estimate.
#' 
#'
#' @examples 
#' trio_example <- matrix(c(1,1,2,0, 1,2,1,0, 1,1,1,0), nrow=3, byrow=TRUE) 
#' y_quant <- c(11.8,1.7,3.6,12.2)
#' list_assign <- assignment_index(trio_example)
#' mu_hat <- calc_mu_hat(y_quant, list_assign$W0, list_assign$W1)
#' 
#' @export
calc_mu_hat <- function(Y, W0, W1){
  mu_hat <- (sum(Y*W0) + sum(Y*W1)) / (sum(W0) + sum(W1))
  return(mu_hat)
}







