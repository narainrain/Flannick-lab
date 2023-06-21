"""
Simulation Arguments
====================================
This is an file used that contains the arguments used in the simulation pipeline
There are two functions here: arg_generate_causal_snps and arg_simulations_pipeline

The arg_generate_causal_snps is used before the arg_simulations_pipeline to generate the causal snps
that will be the input to the simulations_pipeline. The different non-zero arguments are required for
different models

Null Case (probability true betas = 0):
    p-causal-gene = 0
    p-causal-snp = 0
Infinitesimal (probability true betas = 1):
    p-causal-snp = 1
Non-Infinitesimal (0 < probability true betas < 1):
    0 < p-casual-snp < 1
Non-Infinitesimal Snp Cluster (0 < probability causal gene < 1 && 0 < probability causal snp in gene < 1):
    0 < p-causal-gene < 1
    0 < p-snps-in-causal-gene < 1
PEGS (Same as Non-Infinitesimal Snp Cluster + probability of causal snp outside causal gene):
    0 < p-causal-gene < 1
    0 < p-snps-in-causal-gene < 1
    0 < p-snps-out-causal-gene < 1
"""

# Author: Sumit Narain <narain@broadinstitute.org>
#
# License:

from optparse import OptionParser


def arg_generate_causal_snps():
    parser = OptionParser("usage: %prog [options]")
    # General options  ====================================================
    parser.add_option("", "--output-file", default="input_file")
    # the --output file is the name of the file for causal snps which will be used
    # as input for the simulations class.
    parser.add_option("", "--p-causal-gene", type="float", default=0.05)  # causal genes throughout genome
    parser.add_option("", "--p-snps-in-causal-gene", type="float", default=0.0014)  # causal snps inside the causal genes
    parser.add_option("", "--p-snps-out-causal-gene", type="float", default=0.0001)  # causal snps outside the causal genes
    parser.add_option("", "--p-causal-snps", type="float", default=1)  # causal snps throughout genome
    parser.add_option("", "--window", default=400, type="int")  # the window size for each gene
    parser.add_option("", "--gene-loc", default="./raw/gene_loc_37.tab")  # File that contains the gene location information
    parser.add_option("", "--raw-file-location", default="./raw/")  # the location of the ld files
    parser.add_option("", "--chrom-start", default="22", type="int")  # This is used in case a specific chromosome is wanted
    # the default chromosome should be 1 to cover all chromosomes from 1-22
    parser.add_option("", "--replicates", default="1", type="int")  # This is used in case a specific chromosome is wante
    (options, _) = parser.parse_args()
    return options


def arg_simulations_pipeline():
    parser = OptionParser("usage: %prog [options]")
    # General options  ====================================================
    parser.add_option("", "--bim-file", default="./raw/1000G.EUR.QC.22.bim")  # all SNPs in the region
    parser.add_option("", "--ld-file", default="./raw/22.ld")  # in row form
    parser.add_option("", "--ld-threshold", default=0, type="float")
    parser.add_option("", "--h2", default=0.1212, type="float")  # 1-h2 is added as phenotypic variance
    parser.add_option("", "--N", default=100000, type="float")  # sample size for sumstat simulations
    # N has a standard deviation of 2000
    parser.add_option("", "--total_snp_genome", default=10000000, type="int")  # total number of snps in genome
    parser.add_option("", "--beta-file", default=None)  # specify mapping from SNP to true beta
    parser.add_option("", "--beta-dist-file", default=None)  # specify a distribution of beta values to be assigned to
    # each SNP within a region. No header, columns are chrom, begin, end, p, mean, var [replicate] where if the optional
    # [1-based] replicate column is present, it will apply the region only to the specified replicate. If multiple lines
    # apply to an SNP, it will take the first one
    parser.add_option("", "--input-causal-snps", default="input_file.causal_snps") # Input file of causal snps
    parser.add_option("", "--max-component-size", type="int", default=150)  # control the maximum number of SNPs to be included in an LD-block (smaller is
    # faster)
    parser.add_option("", "--assume-r", action="store_true")
    parser.add_option("", "--num-sim", type="int", default=1)  # generate statistics for multiple independent replicates
    parser.add_option("", "--ldsc-format", action="store_true")  # output sumstats in LDSC format
    parser.add_option("", "--num-causal-snps-out",
                      default=None)  # write a file with the number of causal SNPs per iteration
    parser.add_option("", "--output-file", default="output_file") # this is used to name the output-file
    parser.add_option("", "--file-path", default=None)  # specify the path of the file for outputs
    parser.add_option("", "--p", type="float", default=0.0000007)  # Proportion of values going into the random distribution
    parser.add_option("", "--chrom", default=22, type="int")  # Tells use the chromosome number that we're examining
    (options, _) = parser.parse_args()
    return options


