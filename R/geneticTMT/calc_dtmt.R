# Transmission-Mean Test Statistic
#'
#' This function computes the unbiased TMT estimand given trio genotypes and 
#' children's phenotype. The phenotype can be either quantitative or dichotomous.
#'
#' @param Y The n-length vector of children' phenotype.
#' @param is_phased Label for unphased/phased genotype. If `FALSE` (default), the
#'     child's genotype is unphased. Else, the child's genotype is phased.
#' @param trio_geno A \eqn{3 \times n}{3 x n} matrix of trio genotype if `FALSE`
#'     for is_phased. Or a \eqn{4 \times n}{4 x n} matrix of trio genotype if
#'     `TRUE` for is_phased. Must contain only values \eqn{(0, 1, 2)}{0, 1, 2}.
#'     The first row is assumed maternal, the second row is assumed paternal, and
#'     the third (or together with the fourth) row is assumed child's genotypes.
#'
#' @return A numeric value of the TMT statistic.
#' 
#' @note Undesired dimensions of the genotype matrices will result in errors.
#'
#' @examples 
#' trio_example <- matrix(c(1,1,2,0, 1,2,1,0, 1,1,1,0), nrow=3, byrow=TRUE) 
#' y_quant <- c(11.8,1.7,3.6,12.2)
#' dtmt <- calc_dtmt(y_quant, trio_example)
#' 
#' @export
calc_dtmt <- function(Y, trio_geno, is_phased=FALSE){
  # construct assignment vectors
  list_assign <- assignment_index(trio_geno, is_phased)
  W0 <- list_assign$W0
  W1 <- list_assign$W1
  # estimate mu to center trait value
  mu_hat <- calc_mu_hat(Y, W0, W1)
  # sizes of treatment and control
  N <- sum(W1) + sum(W0)
  # calculate the tmt estimand
  dtmt <- ifelse(N>0, 2/N*(sum((Y-mu_hat)*W1) - sum((Y-mu_hat)*W0)), NaN)
  return(dtmt)
}







