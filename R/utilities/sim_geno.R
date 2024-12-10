#' Simulate genetic parameters under the BN-PSD model
#'
#' This function generates genetic parameters of structured population.
#' The simulation is constructed on the BN-PSD model.
#' Given parameters p and f, simulate the individual-allele-frequency pi by BN(p,f)
#'
#' @param m The total number of genotype loci.
#' @param n The total number of individuals.
#' @param k An integer as the number of intermediate sub-populations.
#' @param Fst The scalar as the overall fixation index of the unerlying population.
#' 
#' @importFrom bnpsd admix_prop_1d_linear draw_p_anc draw_p_subpops fst_admix
#' 
#' @return A list of genetic parameters for bnpsd
#'
#' @example
#' param_list <- sim_geno_param(m=100000, n=5000, k=4, Fst=0.2)
#'
#'
#' @export
sim_geno_param <- function(m, n, k=4, Fst = 0.2){
  # simulate individual allele frequencies given a numeric number of FST
  # this FST is the overall fixation index of all subpopulations
  # draw ancestral allele frequency ############################################
  p_t <- bnpsd::draw_p_anc(m_loci=m, p_min=0.1, p_max=0.9)
  # admixed proportions and fst for each subpopulation #########################
  # sub-population index for simulating admixed proportions
  index_subpops <- 1:k
  # desired FST is Fst
  # simulate admixture proportions from linear-1D geography
  param_struct <- bnpsd::admix_prop_1d_linear(n_ind = n, k_subpops = k,
                                              bias_coeff=0.5,
                                              coanc_subpops=index_subpops,
                                              fst = Fst)
  q <- param_struct$admix_proportions
  f_Su_T <- param_struct$coanc_subpops
  FST <- bnpsd::fst_admix(q, f_Su_T)
  # draw allele frequency in the subpopulations ################################
  p_s <- bnpsd::draw_p_subpops(p_anc=p_t, inbr_subpops = f_Su_T)
  # wrap up the output
  list_out <- list(p_t, p_s, q, FST)
  names(list_out) <- c('p_t', 'p_s', 'q', 'FST')
  return(list_out)
}


#' Simulate genotype matrix given the genetic parameters under the BN-PSD model
#'
#' This function generates synthetic genotype matrix of structured population.
#' The simulation is constructed on the BN-PSD model.
#' Given parameters p and f, simulate the individual-allele-frequency pi by BN(p,f)
#' Then simulate biallelic genotype `x_ij` by Binomial(2,pi)
#'
#' @param p_s A matrix of allele frequencies among intermediate sub-populations.
#' @param q A matrix of the admixture proportions.
#' 
#' @importFrom bnpsd draw_genotypes_admix
#' 
#' @return A matrix of genotypes.
#'
#' @examples
#' param_list <- sim_geno_param(m=100000, n=5000, k=4, Fst=0.2)
#' x <- sim_geno(p_s=param_list$p_s, q=param_list$q)
#'
#'
#' @export
sim_geno <- function(p_s, q){
  x_geno <- bnpsd::draw_genotypes_admix(p_s, q)
  return(x_geno)
}



sim_allele_single <- function(geno_entry){
  if(geno_entry == 1){
    return(rbinom(n=1,size=1,prob=0.5))
  } else {
    return(geno_entry/2)
  }
}

sim_allele_row <- function(geno_row){
  allele_row <- lapply(geno_row, FUN=sim_allele_single)
  return(unlist(allele_row))
}

sim_allele <- function(geno, n_cores=18, fix_seed=FALSE, seed_start=1){
  if(fix_seed==FALSE){
    sim_allele_row <- function(i){
      allele_row <- lapply(geno[i,], FUN=sim_allele_single)
      return(unlist(allele_row))
    }
    allele <- parallel::mclapply(c(1:nrow(geno)), sim_allele_row, mc.cores=n_cores)
    allele <- matrix(unlist(allele), ncol=ncol(geno), byrow=TRUE)
    return(allele)
  } else {
    set.seed(seed_start)
    seed_list <- sample(c(1:nrow(geno)), size=nrow(geno)) + seed_start
    sim_allele_row <- function(i){
      set.seed(seed_list[i])
      allele_row <- lapply(geno[i,], FUN=sim_allele_single)
      return(unlist(allele_row))
    }
    allele <- parallel::mclapply(c(1:nrow(geno)), sim_allele_row, mc.cores=n_cores)
    allele <- matrix(unlist(allele), ncol=ncol(geno), byrow=TRUE)
    return(allele)
  }
}

calc_corr_square_snp <- function(snp_corr){
  return(snp_corr**2)
}

calc_corr_square_geno <- function(geno){
  coarse_snp <- function(snp){
    snp_tail <- snp[length(snp)]
    if(sd(snp)==0){snp_tail <- snp_tail + 0.5}
    return(snp_tail)
  }
  geno[,ncol(geno)] <- apply(geno, 1, FUN=coarse_snp)
  geno_corr <- cor(t(geno))
  geno_corr_square <- apply(geno_corr, 1:2, FUN=calc_corr_square_snp)
  return(geno_corr_square)
}


scan_ld <- function(geno, causal_loci, test_loci, thred=0.2){
  geno_t <- geno[test_loci,]
  geno_o <- geno[causal_loci,]
  coarse_snp <- function(snp){
    snp_tail <- snp[length(snp)]
    if(sd(snp)==0){snp_tail <- snp_tail + 0.5}
    return(snp_tail)
  }
  geno_t[,ncol(geno_t)] <- apply(geno_t, 1, FUN=coarse_snp)
  geno_o[,ncol(geno_o)] <- apply(geno_o, 1, FUN=coarse_snp)
  
  geno_corr <- cor(t(geno_t), t(geno_o))
  snp_ld <- function(entry){return(abs(entry)>thred)}
  geno_ld <- apply(geno_corr, 1:2, FUN=snp_ld)
  return(test_loci[which(rowSums(geno_ld)>0)])
}



