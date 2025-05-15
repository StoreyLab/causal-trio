#' Assignment Index
#'
#' This function constructs assignment indices based on trio genotypes. 
#'
#' @param is_phased Label for unphased/phased genotype. If `FALSE` (default), the
#'     child's genotype is unphased. Else, the child's genotype is phased.
#' @param trio_geno A \eqn{3 \times n}{3 x n} matrix of trio genotype if `FALSE`
#'     for \code{is_phased}. Or a \eqn{4 \times n}{4 x n} matrix of trio genotype if
#'     `TRUE` for \code{is_phased}. Must contain only values \eqn{(0, 1, 2)}{0, 1, 2}.
#'     The first row is assumed maternal, the second row is assumed paternal, and
#'     the third (or together with the fourth) row is assumed child's genotypes.
#'
#' @return A named list of two vectors containing indicators for assigning to 
#'     treatment group 0 or treatment group 1.
#'     
#' @note Undesired dimensions of the genotype matrices will result in errors.
#' 
#' @examples
#' trio_example <- matrix(c(1,1,2,0, 1,2,1,0, 1,1,1,0), nrow=3, byrow=TRUE)
#' list_assign <- assignment_index(trio_example, is_phased=FALSE)
#' 
#' @export
assignment_index <- function(trio_geno, is_phased=FALSE){
  if(is_phased){
    # validate input
    if(nrow(trio_geno)!=4) {
      stop('undesired dimension for phased trio genotypes')
    }
    # assignment index
    z_m <- trio_geno[1,]
    z_p <- trio_geno[2,]
    a_m <- trio_geno[3,]
    a_p <- trio_geno[4,]
    w1m <- as.numeric((z_m==1) & (a_m==1))
    w0m <- as.numeric((z_m==1) & (a_m==0))
    w1p <- as.numeric((z_p==1) & (a_p==1))
    w0p <- as.numeric((z_p==1) & (a_p==0))
    w1 <- w1m + w1p
    w0 <- w0m + w0p
    # output
    list_assign <- list(w1, w0, w1m, w0m, w1p, w0p)
    names(list_assign) <- c('W1','W0','W1m','W0m','W1p','W0p')
    return(list_assign)
  } else {
    # validate input
    if(nrow(trio_geno)!=3) {
      stop('undesired dimension for unphased trio genotypes')
    }
    # assignment index
    z_m <- trio_geno[1,]
    z_p <- trio_geno[2,]
    g <- trio_geno[3,]
    # assign to control
    w0 <- as.numeric((z_m==2) & (z_p==1) & (g==1)) +
          as.numeric((z_m==1) & (z_p==2) & (g==1)) +
          as.numeric((z_m==1) & (z_p==1) & (g==1)) +
          2 * as.numeric((z_m==1) & (z_p==1) & (g==0)) +
          as.numeric((z_m==1) & (z_p==0) & (g==0)) +
          as.numeric((z_m==0) & (z_p==1) & (g==0)) 
    # assign to treatment
    w1 <- as.numeric((z_m==2) & (z_p==1) & (g==2)) +
          as.numeric((z_m==1) & (z_p==2) & (g==2)) +
          2 * as.numeric((z_m==1) & (z_p==1) & (g==2)) +
          as.numeric((z_m==1) & (z_p==1) & (g==1)) +
          as.numeric((z_m==1) & (z_p==0) & (g==1)) +
          as.numeric((z_m==0) & (z_p==1) & (g==1))
    # output
    list_assign <- list(w1, w0)
    names(list_assign) <- c('W1','W0')
    return(list_assign)
  }
}






