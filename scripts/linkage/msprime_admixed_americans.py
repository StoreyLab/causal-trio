import sys
import os
import msprime
import numpy as np

'''
# use this script

module load python/3.9.2
python3 msprime_admixed_americans.py 5000 100000

'''

def sim_geno_parent(n_ind, m_snp, no_chr):
    '''simulate phased parental genotype per chromosome'''
    # population structure
    T_OOA = 920
    demography = msprime.Demography()
    demography.add_population(name="AFR", description="African",    initial_size=14474)
    demography.add_population(name="EUR", description="European",   initial_size=34039, growth_rate=0.0038)
    demography.add_population(name="EAS", description="East Asian", initial_size=45852, growth_rate=0.0048)
    demography.add_population(name="AMA", description="Admixed America", initial_size=54664, growth_rate=0.05)
    demography.add_population(name="OOA", description="Bottleneck out-of-Africa", initial_size=1861)
    demography.add_population(name="AMH", description="Anatomically modern humans", initial_size=14474)
    demography.add_population(name="ANC", description="Ancestral equilibrium", initial_size=7310)
    demography.set_symmetric_migration_rate(["AFR", "EUR"], 2.5e-5)
    demography.set_symmetric_migration_rate(["AFR", "EAS"], 0.78e-5)
    demography.set_symmetric_migration_rate(["EUR", "EAS"], 3.11e-5)
    demography.add_admixture(12, derived="AMA", ancestral=["AFR", "EUR", "EAS"], proportions=[1/6, 2/6, 3/6])
    demography.add_population_split(T_OOA, derived=["EUR", "EAS"], ancestral="OOA")
    demography.add_symmetric_migration_rate_change(time=T_OOA, populations=["AFR", "OOA"], rate=15e-5)
    demography.add_population_split(2040, derived=["OOA", "AFR"], ancestral="AMH")
    demography.add_population_split(5920, derived=["AMH"], ancestral="ANC")
    demography.debug()

    sample_list = {"AMA": n_ind*2}
    ts = msprime.sim_ancestry(samples=sample_list, demography=demography, sequence_length=m_snp, random_seed=no_chr+1)
    mts = msprime.sim_mutations(ts, rate=1e-4, random_seed=no_chr+1, model='binary')
    geno = np.array([var.genotypes for var in mts.variants()])
    np.save('./geno_parent_admixed_americans_chr'+str(no_chr+1), geno)


def sim_whole_genome_parent(n_ind, m_snp):
    '''simulate phased parental genotype across all chromosomes'''
    chr_length = [24.9,24.2,19.8,19,18.2,17.1,16,14.5,13.8,13.4,13.5,13.3,11.4,10.7,10.2,9,8.3,8,5.9,6.4,4.7,5.1]
    m_list = [int(np.round(m_snp/np.sum(chr_length)*i)) for i in chr_length[0:21]]
    # reserve 90 spots for direct causal loci without causal linkage
    m_list.append(m_snp - np.sum(m_list))
    for id_chr in range(len(m_list)):
        sim_geno_parent(n_ind, m_list[id_chr], id_chr)
    # assemble simulated genotype
    geno_out = np.load('./geno_parent_admixed_americans_chr1.npy')
    for id_chr in range(21):
        temp = np.load('./geno_parent_admixed_americans_chr'+str(id_chr+2)+'.npy')
        geno_out = np.concatenate((geno_out, temp), axis=0)
    np.save('./geno_parent', geno_out)


if __name__ == '__main__':
    n_ind = int(sys.argv[1])
    m_snp = int(sys.argv[2])

    sim_whole_genome_parent(n_ind, m_snp)
