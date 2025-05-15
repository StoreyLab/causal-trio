#' Transmission-Disequilibrium Test
#'
#' This function calculates the \eqn{\chi^2} test statistic and p value for TDT.
#'
#'
#' @param Y The n-length vector of children' quantitative trait
#' @param is_phased Label for unphased/phased genotype. If `FALSE` (default), the
#'     child's genotype is unphased. Else, the child's genotype is phased.
#' @param trio_geno A \eqn{3 \times n}{3 x n} matrix of trio genotype if `FALSE`
#'     for is_phased. Or a \eqn{4 \times n}{4 x n} matrix of trio genotype if
#'     `TRUE` for is_phased. Must contain only values \eqn{(0, 1, 2)}{0, 1, 2}.
#'     The first row is assumed maternal, the second row is assumed paternal, and
#'     the third (or together with the fourth) row is assumed child's genotype.
#'
#' @importFrom stats pchisq
#'
#' @return A list of p value and McNemar's test statistic.
#' 
#' @note Undesired dimensions of the genotype matrices will result in errors.
#'
#' @examples 
#' trio_example <- matrix(c(1,1,2,0, 1,2,1,0, 1,1,1,0), nrow=3, byrow=TRUE) 
#' y_binary <- c(0,1,1,1)
#' tdt_list <- tdtest(y_binary, trio_example)
#' 
#' @export 
tdtest <- function(Y, trio_geno, is_phased=FALSE){
  if(!all(Y %in% 0:1)){stop("Y must be a binary vector")}
  list_assign <- assignment_index(trio_geno, is_phased)
  # assign trait value to treatment and control group
  vec_n1 <- rep(Y, list_assign$W1)
  vec_n0 <- rep(Y, list_assign$W0)
  n1 <- sum(vec_n1==1)
  n0 <- sum(vec_n0==1)
  # compute chi2
  chi2 <- (n1 - n0)^2 / (n1 + n0)
  p_value <- 1 - stats::pchisq(chi2, 1)
  tdt_list <- list(p_value, chi2)
  names(tdt_list) <- c('p','chisqr')
  return(tdt_list)
}



