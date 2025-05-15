# The Transmission-Mean Test (TMT)
#'
#' This function calculates the p-value for the TMT given trio genotypes and 
#' children's phenotype. The phenotype can be either quantitative or dichotomous.
#'
#' @param Y The n-length vector of children' quantitative trait
#' @param is_phased Label for unphased/phased genotype. If `FALSE` (default), the
#'     child's genotype is unphased. Else, the child's genotype is phased.
#' @param trio_geno A \eqn{3 \times n}{3 x n} matrix of trio genotype if `FALSE`
#'     for is_phased. Or a \eqn{4 \times n}{4 x n} matrix of trio genotype if
#'     `TRUE` for is_phased. Must contain only values \eqn{(0, 1, 2)}{0, 1, 2}.
#'     The first row is assumed maternal, the second row is assumed paternal, and
#'     the third (or together with the fourth) row is assumed child's genotypes.
#' @param test A character string of the test type. 
#' @param n_per A numeric number for the permutation times if \code{test=}'Permutation'
#'
#' @importFrom stats pnorm
#'
#' @return A numeric result for the p-value.
#' 
#' @note Undesired dimensions of the genotype matrices and any test type other 
#'     than 'TMT' or 'Permutation' for \code{test} will result in errors. 
#' 
#' @examples
#' trio_example <- matrix(c(1,1,2,0, 1,2,1,0, 1,1,1,0), nrow=3, byrow=TRUE) 
#' y_quant <- c(11.8,1.7,3.6,12.2)
#' tmt_p <- tmtest(y_quant, trio_example, test='TMT')
#' 
#' tmt_p_perm <- tmtest(y_quant, trio_example, test='Permutation', n_per=1000)
#' 
#' @export 
tmtest <- function(Y, trio_geno, is_phased=FALSE,
                    test=c('TMT','Permutation'), n_per=100){
  test <- match.arg(test)
  dtmt <- calc_dtmt(Y, trio_geno, is_phased)
  varhat <- calc_varhat_dtmt(Y, trio_geno, is_phased)
  
  if(test == 'TMT'){
    tau <- dtmt / sqrt(varhat)
    pvalue <- 2 * stats::pnorm(-abs(tau))
  } else if (test == 'Permutation'){
    if( is.na(dtmt) ){
      pvalue <- NaN
    } else {
      perm_loop <- function(i){
        set.seed(i)
        dtmt_per <- calc_dtmt(sample(Y, size=length(Y), replace=FALSE),
                              trio_geno, is_phased)
        return( abs(dtmt_per) >= abs(dtmt) )
      }
      p_perm <- lapply(c(1:n_per), FUN=perm_loop)
      p_perm <- unlist(p_perm)
      pvalue <- sum(p_perm) / length(p_perm)
    }
  }
  return(pvalue)
}






