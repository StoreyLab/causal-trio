# Transmission-Mean Test Centered Statistic 
#'
#' This function computes the centered TMT statistic given trio genotypes and 
#' children's phenotype. The phenotype can be either quantitative or dichotomous.
#'
#' @param Y The n-length vector of children' phenotype.
#' @param is_phased Label for unphased/phased genotype. If `FALSE` (default), the
#'     child's genotype is unphased. Else, the child's genotype is phased.
#' @param mu A numeric value of the population mean. 
#' @param trio_geno A \eqn{3 \times n}{3 x n} matrix of trio genotype if `FALSE`
#'     for is_phased. Or a \eqn{4 \times n}{4 x n} matrix of trio genotype if
#'     `TRUE` for is_phased. Must contain only values \eqn{(0, 1, 2)}{0, 1, 2}.
#'     The first row is assumed maternal, the second row is assumed paternal, and
#'     the third (or together with the fourth) row is assumed child's genotypes.
#'
#' @return A numeric value of the centered TMT statistic.
#' 
#' @note Undesired dimensions of the genotype matrices will result in errors.
#'
#' @examples 
#' trio_example <- matrix(c(1,1,2,0, 1,2,1,0, 1,1,1,0), nrow=3, byrow=TRUE) 
#' list_assign <- assignment_index(trio_example, is_phased=FALSE)
#' y_quant <- c(11.8,1.7,3.6,12.2)
#' dtmt_pc <- calc_dtmt_pc(y_quant, trio_example, mu=5)
#' 
#' @export
calc_dtmt_pc <- function(Y, trio_geno, mu, is_phased=FALSE){
  # construct assignment vectors
  list_assign <- assignment_index(trio_geno, is_phased)
  W0 <- list_assign$W0
  W1 <- list_assign$W1
  # sizes of treatment and control
  N <- sum(W1) + sum(W0)
  # calculate the tmt estimand
  dtmt <- ifelse(N>0, 2/N*(sum((Y-mu)*W1) - sum((Y-mu)*W0)), NaN)
  return(dtmt)
}







