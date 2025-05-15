#' Transmission-Disequilibrium Test Statistic
#'
#' This function calculates the TDT statistic given trio genotypes and children's 
#' dichotomous phenotype.
#'
#'
#' @param Y The n-length vector of the children's dichotomous trait
#' @param is_phased Label for unphased/phased genotype. If `FALSE` (default), the
#'     child's genotype is unphased. Else, the child's genotype is phased.
#' @param trio_geno A \eqn{3 \times n}{3 x n} matrix of trio genotype if `FALSE`
#'     for is_phased. Or a \eqn{4 \times n}{4 x n} matrix of trio genotype if
#'     `TRUE` for is_phased. Must contain only values \eqn{(0, 1, 2)}{0, 1, 2}.
#'     The first row is assumed maternal, the second row is assumed paternal, and
#'     the third (or together with the fourth) row is assumed child's genotype.
#'
#' @return A numeric value of the TDT statistic.
#' 
#' @note Inputting phenotypes that are not dichotomous and undesired dimensions
#'     of the genotype matrices will result in errors.
#'
#' @examples 
#' trio_example <- matrix(c(1,1,2,0, 1,2,1,0, 1,1,1,0), nrow=3, byrow=TRUE) 
#' y_binary <- c(0,1,1,0)
#' dtdt <- calc_dtdt(y_binary, trio_example)
#' 
#' @export
calc_dtdt <- function(Y, trio_geno, is_phased=FALSE){
  if(!all(Y %in% 0:1)){stop("Y must be a binary vector")}
  list_assign <- assignment_index(trio_geno, is_phased)
  W0 <- list_assign$W0
  W1 <- list_assign$W1
  # sizes of treatment and control
  N <- sum(W1) + sum(W0)
  # assign trait value to treatment and control group
  vec_n1 <- rep(Y, list_assign$W1)
  vec_n0 <- rep(Y, list_assign$W0)
  n1 <- sum(vec_n1==1)
  n0 <- sum(vec_n0==1)
  ttdt <- ifelse(N>0, 2*(n1-n0)/N, NaN)
  return(ttdt)
}


