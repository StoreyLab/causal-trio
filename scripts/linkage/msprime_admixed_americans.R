#' Generate Genotypes Based on Haplotypes from msprime
#' 
#' This document load output of msprime to generate genotypes
#' 
n_fam <- 5000
m_snp <- 100000
# import parental haplotype from msprime output
hap_parent <- RcppCNPy::npyLoad('../work/linkage/admixed_americans/geno_parent.npy','integer')
# merge haplotype as genotype
id_m <- seq(1,4*n_fam-1,4) # maternal genotype indices
id_p <- id_m + 2 # paternal genotype indices
Zm <- hap_parent[,id_m] + hap_parent[,id_m+1]
Zp <- hap_parent[,id_p] + hap_parent[,id_p+1]
# check parental kinship
kin_m <- popkin::popkin(Zm)
kin_p <- popkin::popkin(Zp)
fst_m <- popkin::fst(kin_m) # around 0.17
fst_p <- popkin::fst(kin_p)