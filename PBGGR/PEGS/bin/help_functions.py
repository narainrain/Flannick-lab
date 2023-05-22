import math
import sys
from optparse import OptionParser
import numpy as np
import scipy.stats
import scipy as sc
import gzip
from itertools import chain, combinations
import scipy.cluster.hierarchy as sch
import os
# import seaborn as sns
import matplotlib.pyplot as plt
from itertools import cycle
from itertools import product
from scipy.stats import multivariate_normal, invwishart
from sklearn.linear_model import LogisticRegression
import plotly.graph_objects as go
import warnings
from copy import deepcopy
import matplotlib
matplotlib.use('Qt5Agg')

#  Global variables  ---------------------------------------------------------------------------------------------------
# log levels
NONE = 0
INFO = 1
DEBUG = 2
TRACE = 3
num_lik_calc = 0
iterations = 0
# arguments
global options
# files handlers
global warnings_fh
global diagnostics_fh
# ld variables
global ld_header
global have_negative
# process gene variables
global cache
# file for ldpred results
global file_ldpred


# ----------------------------------------------------------------------------------------------------------------------


def initialize_global_var():
    global cache
    global options
    cache = {}
    global file_ldpred
    file_ldpred = open('../data/ldpred_results/ldpred.results', 'w')
    file_ldpred.write("MarkerName\tchr\tpos\tref\talt\tposterior\tp_vec\tzscore\tabs_pos\n")
    # ================================================================
    options.causal_window = options.causal_window * 1000
    a = options.heritability_causal
    M_1 = options.total_num_SNPs * options.p
    if options.true_beta:
        # options.sigma = 2.096024597659547e-6
        # options.null_sigma = 1.4318124378580449e-15
        options.sigma = (options.heritability / M_1)
        options.null_sigma = (1 - a) * (options.heritability / (options.total_num_SNPs - M_1))
    else:
        options.sigma = a * ((options.sample_size * options.heritability) / M_1)
        options.null_sigma = (1 - a) * ((options.sample_size * options.heritability) / (options.total_num_SNPs - M_1))


def arg_settings():
    parser = OptionParser("usage: %prog [options]")
    parser.add_option("", "--warnings-file", default=None)
    parser.add_option("", "--debug-level", type='int', default=1)
    parser.add_option("", "--debug-stdout", type='int', default=2)
    parser.add_option("", "--null-sigma", type=float, default=0)  # the variance of all causal snps under the null model
    parser.add_option("", "--h2-per-causal-gene", type=float, default=None)  # the variance of all causal snps
    parser.add_option("", "--sigma", type=float,
                      default=0.004039)  # the variance of causal snps; h2/M * var(pheno) 0.03146426544510455
    # use window model for causal SNPs. Argument is how far in each direction the window extends
    parser.add_option("", "--causal-window", type=int, default=10)  # in kilo bases
    parser.add_option("", "--extension", type=int, default=20000)
    # use gaussian model for causal SNPs. Arguments are window1 (distance when causality is height1), height1,
    # window2 (distance when probability of causality=height2), and height2
    parser.add_option("", "--causal-gaussian", type=float, nargs=4, action='append', default=None)
    parser.add_option("", "--causal-table", default=None)  # use table of distance to probability
    parser.add_option("", "--causal-table-dist-col", default=1)
    parser.add_option("", "--causal-table-prob-col", default=2)
    parser.add_option("", "--epsilon-ratio", type=float, default=0)  # multiply by sigma to get epsilon
    parser.add_option("", "--stab-term-frac", type=float, default=1e-4)
    parser.add_option("", "--log-diff-ignore", type=float, default=0.01)
    parser.add_option("", "--SNP-log-diff-ignore", type=float, default=0.01)
    parser.add_option("", "--max-num-causal-genes", type=int,
                      default=40)  # maximum number of causal genes in the region
    parser.add_option("", "--gene-prior", type=float, default=0.05)  # prior that each gene is causal
    parser.add_option("", "--range", default=None)  # filter analysis to a chr:start-stop range
    parser.add_option("", "--chrom", default=None)  # filter analysis to a chromosome
    parser.add_option("", "--start", type=int, default=None)  # filter analysis to a beginning position
    parser.add_option("", "--end", type=int, default=None)  # filter analysis to an ending position
    parser.add_option("", "--gene-file", default=None)  # Must have GENE, CHROM, START, STOP
    parser.add_option("", "--gene-id-col", default=None)
    parser.add_option("", "--gene-chrom-col", default=None)
    parser.add_option("", "--gene-start-col", default=None)
    parser.add_option("", "--gene-end-col", default=None)
    parser.add_option("", "--extra-snp-file", default=None)  # Must have ID, CHROM, POS
    parser.add_option("", "--extra-snp-id-col", default=2)
    parser.add_option("", "--extra-snp-chrom-col", default=1)
    parser.add_option("", "--extra-snp-pos-col", default=4)
    parser.add_option("", "--sumstats-file", default=None)  # Must have ID, BETA, and SE or PVALUE
    parser.add_option("", "--sumstats-id-col", default=None)
    parser.add_option("", "--sumstats-beta-col", default=None)
    parser.add_option("", "--sumstats-se-col", default=None)
    parser.add_option("", "--sumstats-p-col", default=None)
    parser.add_option("", "--sumstats-z-col", default=None)
    parser.add_option("", "--sumstats-n-col", default=None)
    parser.add_option("", "--sumstats-effect-allele-col", default="Allele1")
    parser.add_option("", "--sumstats-other-allele-col", default="Allele2")
    parser.add_option("", "--n-threshold", type=float, default=0)
    parser.add_option("", "--ld-file", default=None)
    parser.add_option("", "--ld-id1-col", default="SNP_A")
    parser.add_option("", "--ld-chrom1-col", default="CHR_A")
    parser.add_option("", "--ld-pos1-col", default="BP_A")
    parser.add_option("", "--ld-id2-col", default="SNP_B")
    parser.add_option("", "--ld-chrom2-col", default="CHR_B")
    parser.add_option("", "--ld-pos2-col", default="BP_B")
    parser.add_option("", "--ld-r-col", default="R2")
    parser.add_option("", "--ld-allele-file", default=None)  # Must have ID, ALLELE
    parser.add_option("", "--ld-allele-id-col", default=2)
    parser.add_option("", "--ld-allele-effect-allele-col", default=5)
    parser.add_option("", "--ld-allele-other-allele-col", default=6)
    parser.add_option("", "--min-ld-threshold", default=0, type="float")
    parser.add_option("", "--max-ld-threshold", default=1, type="float")
    parser.add_option("", "--all-positive", action="store_true")
    parser.add_option("", "--detect-flips", type='float', default=10)
    parser.add_option("", "--debug-make-positive", action="store_true")
    # control the maximum number of SNPs to be included in an LD-block (smaller is faster)
    parser.add_option("", "--max-component-size", type="int", default=100)
    parser.add_option("", "--block-size", type="int", default=1000000)
    parser.add_option("", "--out-diagnostics", default=None)  # write diagnostics
    parser.add_option("", "--additive", type="int",
                      default=0)  # model for overlapping SNPs, Additive or Max
    parser.add_option("", "--model", action="store", type="string", default="greedy")
    parser.add_option("", "--SNP-model", action="store", type="string", default="Greedy_SNP")
    parser.add_option("", "--cluster", action="store_true", default=False)
    parser.add_option("", "--dentist-file", default=None)
    parser.add_option("", "--dentist-thr", type="float", default=float('+inf'))
    parser.add_option("", "--p", type="float", default=0.001)
    parser.add_option("", "--heritability", type="float", default=0.1212)
    parser.add_option("", "--sample-size", type="float", default=84685)
    parser.add_option("", "--total_num_SNPs", type="float", default=84700000)
    parser.add_option("", "--simulate", action="store_true", default=False)  # generates simulated data
    parser.add_option("", "--create-data", action="store_true", default=False)
    parser.add_option("", "--SNPs-simulated", type="int", default=100)
    parser.add_option("", "--SNPs-in-simulated-gene", type="int", default=10)
    parser.add_option("", "--heritability-causal", type="float", default=1.0)
    parser.add_option("", "--job-id", type="int")
    parser.add_option("", "--add-cov-uncertainty", type="int", default=0)
    parser.add_option("", "--true-beta", type="int", default=0)
    parser.add_option("", "--magma-out", default="test")
    parser.add_option("", "--locus-out", default="test")
    global options
    (options, _) = parser.parse_args()
    return options


def print_model_setup(level):
    log("Model setup", level)
    log("Total num of SNPs: %s -- Heritability: %s -- Sample size: %s"
        % (options.total_num_SNPs, options.heritability, options.sample_size), level)
    log("Null sigma: %s -- Causal sigma: %s" % (options.null_sigma, options.sigma), level)
    log("Fraction of causal SNPs, i.e. p: %s" % options.p, level)
    log("Percentage of heritability explain by causal SNPs: %s" % options.heritability_causal, level)
    log("Extension: %s -- range: %s" % (options.extension, options.range), level)
    log("Causal window: %s" % options.causal_window, level)
    log("Gene prior (pi): %s" % options.gene_prior, level)
    log("Max component size: %s" % options.max_component_size, level)
    if options.SNP_model == "All_one":
        log("Using All one SNP configuration", level)
    elif options.SNP_model == "Full_SNP":
        log("Using Full SNP configuration, i.e. All possible combinations", level)
    else:
        log("Using SNP greedy method with SNP-log-diff-ignore: %s, number of log lik calculated: %s, iterations: %s"
            % (options.SNP_log_diff_ignore, num_lik_calc, iterations), level)
    if options.dentist_file is not None:
        log("SNPs filtered using DENTIST threshold: %s " % options.dentist_thr, level)
    if options.additive:
        log("Using Additive model", level)
    else:
        log("Using Maximum model", level)
    if options.model == "greedy":
        log("Using Greedy method with log-diff-ignore: %s" % options.log_diff_ignore, level)
    else:
        log("Using Full method", level)


def arg_validation():
    if options.ld_file is None:
        bail("Need --ld-file")


def bail(message):
    sys.stderr.write("%s\n" % message)
    sys.exit(1)


def open_files():
    # warning file
    global warnings_fh
    if options.warnings_file is not None:
        warnings_fh = open(options.warnings_file, 'w')
    else:
        warnings_fh = sys.stderr
    # diagnostics file
    global diagnostics_fh
    diagnostics_fh = None
    if options.out_diagnostics is not None:
        try:
            diagnostics_fh = open(options.out_diagnostics, 'w')
            diagnostics_fh.write(
                "Component\tGene_Set\tGene_Pos\tSNP\tPosition\tDistance"
                "\tZ\tScale\tLikelihood\tNull_Scale\tNull_Likelihood"
                "\tDiff_Likelihood\tFactor\tNull_Factor\n")
        except ValueError:
            bail("Failed to open diagnostics file")


def warn(message):
    if warnings_fh is not None:
        warnings_fh.write("Warning: %s\n" % message)
        warnings_fh.flush()


def log(message, level=INFO, end_char='\n'):
    log_fh = sys.stdout
    if level <= options.debug_level:
        log_fh.write("%s%s" % (message, end_char))
        log_fh.flush()


def open_gz(file):
    if file[-3:] == ".gz":
        return gzip.open(file, 'rt')
    else:
        return open(file)


def get_col(header_cols, col_name_or_index, require_match=True):
    try:
        if col_name_or_index is None:
            raise ValueError
        col = int(col_name_or_index) - 1
        if col < 0:
            bail("Column input must be 1-index not 0-index")
        return col
    except ValueError:
        matching_cols = [i for i in range(0, len(header_cols)) if header_cols[i] == col_name_or_index]
        if len(matching_cols) == 0:
            if require_match:
                bail("Could not find match for column %s in header: %s" % (col_name_or_index, "\t".join(header_cols)))
            else:
                return None
        if len(matching_cols) > 1:
            bail("Found two matches for column %s in header: %s" % (col_name_or_index, "\t".join(header_cols)))
        return matching_cols[0]


def get_value(f_cols, column):
    if column >= len(f_cols):
        raise IndexError
    else:
        return f_cols[column]


def clean_chrom(f_chrom):
    if len(f_chrom) > 3 and f_chrom[:3] == "chr":
        f_chrom = f_chrom[3:]
    return f_chrom


def set_range_restrictions():
    restrict_chrom = options.chrom
    restrict_start = options.start
    restrict_end = options.end
    if options.range:
        if restrict_chrom is None or restrict_start is None or restrict_end is None:
            restrict_range = options.range.split(':')
            if len(restrict_range) != 2:
                bail("Error: Couldn't parse range %s (must be of form chr:start-end)" % restrict_range)
            restrict_chrom = clean_chrom(restrict_range[0])
            restrict_range = restrict_range[1].split('-')
            if len(restrict_range) != 2:
                bail("Error: Couldn't parse range %s (must be of form chr:start-end)" % restrict_range)
            restrict_start = int(restrict_range[0])
            restrict_end = int(restrict_range[1])
        else:
            warn("Ignoring --range because --chrom, --start, and --end have all been specified")
    if options.cluster:
        # instance_id = int(os.environ.get('SGE_TASK_ID'))
        instance_id = options.job_id
        instance_start = ((instance_id - 1) * options.block_size)
        instance_end = (instance_id * options.block_size)
        user_range = {restrict_chrom: (instance_start, instance_end)}
        extended_range = {restrict_chrom: (instance_start - options.extension, instance_end + options.extension)}
    else:
        user_range = {restrict_chrom: (restrict_start, restrict_end)}
        extended_range = {restrict_chrom: (restrict_start - options.extension, restrict_end + options.extension)}
    return user_range, extended_range


def read_sumstats():
    snp_to_effect_allele = {}
    snp_to_other_allele = {}
    snp_to_z = {}
    snp_to_n = {}
    snp_to_beta = {}
    snp_to_se = {}
    snp_to_p = {}
    log("Reading sumstats...", DEBUG)
    # in this block we get SNPs alleles, z and n values
    with open_gz(options.sumstats_file) as sumstats_fh:
        sumstats_header = sumstats_fh.readline().strip().split()
        sumstats_id_col = get_col(sumstats_header, options.sumstats_id_col)
        sumstats_beta_col = get_col(sumstats_header, options.sumstats_beta_col, require_match=False)
        sumstats_z_col = get_col(sumstats_header, options.sumstats_z_col, require_match=False)
        if sumstats_beta_col is None and sumstats_z_col is None:
            bail("Require --sumstats-beta-col or --sumstats-z-col")
        sumstats_se_col = get_col(sumstats_header, options.sumstats_se_col, require_match=False)
        sumstats_n_col = get_col(sumstats_header, options.sumstats_n_col, require_match=False)
        sumstats_p_col = get_col(sumstats_header, options.sumstats_p_col, require_match=False)

        sumstats_effect_allele_col = get_col(sumstats_header, options.sumstats_effect_allele_col, require_match=False)
        sumstats_other_allele_col = get_col(sumstats_header, options.sumstats_other_allele_col, require_match=False)

        if (sumstats_beta_col is not None and sumstats_se_col and sumstats_p_col is None and sumstats_n_col is None) \
                or (sumstats_z_col is None and sumstats_se_col is None and sumstats_n_col is None):
            bail("Require either --sumstats-se-col or --sumstats-p-col or --sumstats-n-col")
        for line in sumstats_fh:
            cols = line.strip().split()

            z = None
            n = None
            try:
                snp = get_value(cols, sumstats_id_col)
                if sumstats_beta_col is not None:
                    beta = float(get_value(cols, sumstats_beta_col))
                    if sumstats_se_col is not None and sumstats_p_col is not None \
                            and sumstats_n_col is not None and sumstats_z_col is not None:
                        se = float(get_value(cols, sumstats_se_col))
                        p = float(get_value(cols, sumstats_p_col))
                        n = float(get_value(cols, sumstats_n_col))
                        z = float(get_value(cols, sumstats_z_col))
                    elif sumstats_p_col is not None and sumstats_n_col is not None:
                        p = float(get_value(cols, sumstats_p_col))
                        z = abs(scipy.stats.norm.ppf(p / 2))
                        se = abs(beta / z)
                        n = float(get_value(cols, sumstats_n_col))
                    else:
                        warn("Skipping line %s: not enough columns" % line.strip())
                else:  # fixme: check which data can we calculate
                    z = float(get_value(cols, sumstats_z_col))
                    if sumstats_se_col is not None:
                        se = float(get_value(cols, sumstats_se_col))
                    else:
                        n = float(get_value(cols, sumstats_n_col))
                if sumstats_effect_allele_col is not None and sumstats_other_allele_col is not None:
                    snp_to_effect_allele[snp] = get_value(cols, sumstats_effect_allele_col).upper()
                    snp_to_other_allele[snp] = get_value(cols, sumstats_other_allele_col).upper()

            except IndexError:
                warn("Skipping line %s: not enough columns" % line.strip())

            if options.debug_make_positive:
                z = abs(z)

            snp_to_beta[snp] = beta
            snp_to_se[snp] = se
            snp_to_z[snp] = z
            snp_to_n[snp] = n
            snp_to_p[snp] = p
    log("SNPs in sumstats: %s" % len(snp_to_z), DEBUG)
    return snp_to_effect_allele, snp_to_other_allele, snp_to_z, snp_to_n, snp_to_beta, snp_to_se, snp_to_p


def read_frequencies():
    snp_to_frq = {}
    log("Reading frequencies...")
    file = "../data_chrom2/euro.frq"
    with open_gz(file) as sumstats_fh:
        sumstats_header = sumstats_fh.readline().strip().split()
        sumstats_id_col = get_col(sumstats_header, "SNP")
        sumstats_maf_col = get_col(sumstats_header, "MAF", require_match=False)
        for line in sumstats_fh:
            cols = line.strip().split()
            MAF = None
            try:
                snp = get_value(cols, sumstats_id_col)
                MAF = float(get_value(cols, sumstats_maf_col))
            except IndexError:
                warn("Skipping line %s: not enough columns" % line.strip())

            snp_to_frq[snp] = MAF
    return snp_to_frq


def read_dentist(threshold):
    snp_to_p_dentist = {}
    log("Reading dentist...", DEBUG)
    with open_gz(options.dentist_file) as sumstats_fh:
        sumstats_id_col = 0
        sumstats_p_col = 2
        for line in sumstats_fh:
            cols = line.strip().split()
            try:
                snp = get_value(cols, sumstats_id_col)
                p_dentist = float(get_value(cols, sumstats_p_col))
            except IndexError:
                warn("Skipping line %s: not enough columns" % line.strip())
            if p_dentist <= threshold:
                snp_to_p_dentist[snp] = p_dentist
    return snp_to_p_dentist


#  Elimination of SNP based on the sample size
def allele_correction_1(snp_to_z, snp_to_n, max_n):
    to_delete = [snp for snp in snp_to_z if snp_to_n[snp] < options.n_threshold * max_n]
    before_snps = len(snp_to_z)
    for snp in to_delete:
        del snp_to_z[snp]
        del snp_to_n[snp]
    log("Deleted %d SNPs due to sample size filtering at %.3g of max; %d remaining" % (before_snps - len(snp_to_z),
                                                                                       options.n_threshold,
                                                                                       len(snp_to_z)))


# correction of the allele using the correlation matrix
def allele_correction_2(snp_zs, snp_ns, component_snps_obs, null_cov_matrix, status=TRACE):
    log("Allele correction 2 correlation based...", DEBUG)
    if len(snp_zs) > 1:
        for i in sorted(range(len(snp_ns)), key=lambda x: snp_ns[x]):
            log("SNP %s: z=%.3g, n=%.3g" % (component_snps_obs[i], snp_zs[i], snp_ns[i]), status)
            # calculate the conditional distribution
            sigma11 = null_cov_matrix[i, i]
            sigma12 = np.concatenate((null_cov_matrix[i, :i], null_cov_matrix[i, (i + 1):]))
            zs12 = np.concatenate((snp_zs[:i], snp_zs[i + 1:]))
            sigma22 = np.delete(np.delete(null_cov_matrix, i, axis=0), i, axis=1)
            sigma12_dot_sigma22_inv = sigma12.dot(np.linalg.inv(sigma22))
            conditional_mean = sigma12_dot_sigma22_inv.dot(zs12)
            conditional_cov = sigma11 - sigma12_dot_sigma22_inv.dot(sigma12)
            conditional_lik = scipy.stats.norm.logpdf(snp_zs[i], loc=conditional_mean, scale=np.sqrt(conditional_cov))
            conditional_flip_lik = scipy.stats.norm.logpdf(-snp_zs[i], loc=conditional_mean,
                                                           scale=np.sqrt(conditional_cov))
            if options.detect_flips is not None and conditional_flip_lik - conditional_lik > options.detect_flips:
                log("Detected flip", status)
                snp_zs[i] = -snp_zs[i]
                conditional_lik = conditional_flip_lik

            log("marginal: mean=%.3g; var=%.3g; se=%.3g; lik=%.3g" % (
                0, sigma11, np.sqrt(sigma11), scipy.stats.norm.logpdf(snp_zs[i], loc=0, scale=np.sqrt(sigma11))),
                status)
            log("conditional: mean=%.3g; var=%.3g; se=%.3g; lik=%.3g" % (
                conditional_mean, conditional_cov, np.sqrt(conditional_cov), conditional_lik), status)


def complement_allele(allele):
    if allele == 'A' or allele == 'a':
        return 'T'
    elif allele == 'C' or allele == 'c':
        return 'G'
    elif allele == 'G' or allele == 'g':
        return 'C'
    elif allele == 'T' or allele == 't':
        return 'A'
    else:
        return None


def allele_alignment(snp_to_effect_allele, snp_to_other_allele, snp_to_z, snp_to_n):
    log("Allele correction 1 using LD file for allele identification...", DEBUG)
    if len(snp_to_effect_allele) > 0 and options.ld_allele_file:
        with open_gz(options.ld_allele_file) as ld_allele_fh:
            ld_allele_header = ld_allele_fh.readline().strip().split()
            ld_allele_id_col = get_col(ld_allele_header, options.ld_allele_id_col)
            ld_allele_effect_allele_col = get_col(ld_allele_header, options.ld_allele_effect_allele_col)
            ld_allele_other_allele_col = get_col(ld_allele_header, options.ld_allele_other_allele_col)
            for line in ld_allele_fh:
                cols = line.strip().split()
                snp = get_value(cols, ld_allele_id_col)
                effect_allele = get_value(cols, ld_allele_effect_allele_col).upper()
                other_allele = get_value(cols, ld_allele_other_allele_col).upper()
                if snp in snp_to_z and snp in snp_to_effect_allele and snp in snp_to_other_allele:
                    if snp_to_effect_allele[snp] == effect_allele and snp_to_other_allele[snp] == other_allele:
                        pass
                    elif snp_to_effect_allele[snp] == other_allele and snp_to_other_allele[snp] == effect_allele:
                        snp_to_z[snp] = -snp_to_z[snp]
                    else:
                        other_allele_c = complement_allele(other_allele)
                        effect_allele_c = complement_allele(effect_allele)
                        if other_allele_c is not None and effect_allele_c is not None:
                            if snp_to_effect_allele[snp] == effect_allele_c \
                                    and snp_to_other_allele[snp] == other_allele_c:
                                pass
                            elif snp_to_effect_allele[snp] == other_allele_c \
                                    and snp_to_other_allele[snp] == effect_allele_c:
                                snp_to_z[snp] = -snp_to_z[snp]
                            else:
                                warn("Could not match alleles in LD file (%s, %s) "
                                     "to alleles in sumstats file (%s, %s) for SNP %s" % (effect_allele, other_allele,
                                                                                          snp_to_effect_allele[snp],
                                                                                          snp_to_other_allele[snp],
                                                                                          snp))
                                del snp_to_z[snp]
                                del snp_to_n[snp]
                        else:
                            warn("Could not convert alleles in LD file (%s, %s) for SNP %s"
                                 % (effect_allele, other_allele, snp))
                            del snp_to_z[snp]
                            del snp_to_n[snp]


def read_ld_1(snp_to_analyze, extended_range):
    chr_to_snp_pos = {}
    snp_to_chr_pos = {}
    chr_pos_to_snp = {}
    ld_data = []
    ld_dict = {}
    global ld_header
    global have_negative
    have_negative = False
    log("Reading LD pass 1, reading and setting SNP values...", DEBUG)
    # in this block we get the chrom and pos of SNPs
    with open_gz(options.ld_file) as ld_fh:
        ld_header = ld_fh.readline().strip().split()
        ld_id1_col = get_col(ld_header, options.ld_id1_col)
        ld_chrom1_col = get_col(ld_header, options.ld_chrom1_col)
        ld_pos1_col = get_col(ld_header, options.ld_pos1_col)
        ld_id2_col = get_col(ld_header, options.ld_id2_col)
        ld_chrom2_col = get_col(ld_header, options.ld_chrom2_col)
        ld_pos2_col = get_col(ld_header, options.ld_pos2_col)
        ld_r_col = get_col(ld_header, options.ld_r_col)
        for line in ld_fh:
            try:
                cols = line.strip().split()
                value = float(get_value(cols, ld_r_col))
                if value < 0:
                    have_negative = True
                if abs(value) < options.min_ld_threshold:
                    continue
                snp_1 = get_value(cols, ld_id1_col)
                snp_2 = get_value(cols, ld_id2_col)

                if snp_1 not in snp_to_analyze or snp_2 not in snp_to_analyze:
                    continue
                snp_1_chr = clean_chrom(get_value(cols, ld_chrom1_col))
                snp_1_pos = int(get_value(cols, ld_pos1_col))
                snp_2_chr = clean_chrom(get_value(cols, ld_chrom2_col))
                snp_2_pos = int(get_value(cols, ld_pos2_col))

                if snp_1_chr not in extended_range or snp_2_chr not in extended_range:
                    continue
                if snp_1_pos < extended_range[snp_1_chr][0] or snp_2_pos < extended_range[snp_1_chr][0]:
                    continue
                if snp_1_pos > extended_range[snp_1_chr][1] or snp_2_pos > extended_range[snp_1_chr][1]:
                    continue

                snp_to_chr_pos[snp_1] = (snp_1_chr, snp_1_pos)
                snp_to_chr_pos[snp_2] = (snp_2_chr, snp_2_pos)

            except IndexError:
                warn("Skipping line %s: not enough columns" % line.strip())
                continue

            if snp_1_chr not in chr_to_snp_pos:
                chr_to_snp_pos[snp_1_chr] = set()
            chr_to_snp_pos[snp_1_chr].add(snp_1_pos)
            chr_pos_to_snp[(snp_1_chr, snp_1_pos)] = snp_1
            if snp_2_chr not in chr_to_snp_pos:
                chr_to_snp_pos[snp_2_chr] = set()

            chr_to_snp_pos[snp_2_chr].add(snp_2_pos)
            chr_pos_to_snp[(snp_2_chr, snp_2_pos)] = snp_2
            ld_data.append((snp_1, snp_2, value))
            ld_dict[snp_1, snp_2] = value
    if not have_negative and not options.all_positive:
        bail("All r values are positive; did you pass in r2 values? If not, specify --all-positive to proceed")
    log("SNPs in ld, z-score, and dentist in the selected region: %s" % len(snp_to_chr_pos), DEBUG)
    return chr_to_snp_pos, snp_to_chr_pos, chr_pos_to_snp, ld_data, ld_dict


def read_ld_2(ld_data):
    snp_to_component = {}
    component_to_snp = {}
    component = 0
    max_component_size = 2
    log("Reading LD pass 2, components creation...", DEBUG)
    for entry in sorted(ld_data, key=lambda x: -abs(x[2])):
        snp_1 = entry[0]
        snp_2 = entry[1]
        if snp_1 not in snp_to_component and snp_2 not in snp_to_component:
            component += 1
            snp_to_component[snp_1] = component
            snp_to_component[snp_2] = component
            component_to_snp[component] = set()
            component_to_snp[component].add(snp_1)
            component_to_snp[component].add(snp_2)
        elif snp_1 in snp_to_component and snp_2 not in snp_to_component:
            if len(component_to_snp[snp_to_component[snp_1]]) < options.max_component_size:
                snp_to_component[snp_2] = snp_to_component[snp_1]
                component_to_snp[snp_to_component[snp_1]].add(snp_2)
            else:
                component += 1
                snp_to_component[snp_2] = component
                component_to_snp[component] = set()
                component_to_snp[component].add(snp_2)
        elif snp_2 in snp_to_component and snp_1 not in snp_to_component:
            if len(component_to_snp[snp_to_component[snp_2]]) < options.max_component_size:
                snp_to_component[snp_1] = snp_to_component[snp_2]
                component_to_snp[snp_to_component[snp_2]].add(snp_1)
            else:
                component += 1
                snp_to_component[snp_1] = component
                component_to_snp[component] = set()
                component_to_snp[component].add(snp_1)
        elif snp_2 in snp_to_component and snp_1 in snp_to_component \
                and snp_to_component[snp_1] != snp_to_component[snp_2]:
            if len(component_to_snp[snp_to_component[snp_1]]) + len(component_to_snp[snp_to_component[snp_2]]) <= \
                    options.max_component_size:
                component_1 = snp_to_component[snp_1]
                component_2 = snp_to_component[snp_2]
                for snp in component_to_snp[component_2]:
                    snp_to_component[snp] = component_1
                component_to_snp[component_1] = component_to_snp[component_1].union(component_to_snp[component_2])
                component_to_snp.pop(component_2)

        if len(component_to_snp[snp_to_component[snp_1]]) > max_component_size:
            max_component_size = len(component_to_snp[snp_to_component[snp_1]])
    return snp_to_component, component_to_snp


def read_extra_snp(chr_to_snp_pos, snp_to_chr_pos, chr_pos_to_snp, snp_to_component, component_to_snp, component,
                   extended_range):
    if options.extra_snp_file:
        log("Reading extra SNPs...")
        with open_gz(options.extra_snp_file) as extra_snp_fh:
            extra_snp_header = extra_snp_fh.readline().strip().split()
            extra_snp_id_col = get_col(ld_header, options.extra_snp_id_col)
            extra_snp_chrom_col = get_col(ld_header, options.extra_snp_chrom_col)
            extra_snp_pos_col = get_col(ld_header, options.extra_snp_pos_col)

            for line in extra_snp_fh:
                cols = line.strip().split()
                try:
                    snp = get_value(cols, extra_snp_id_col)
                    chr = clean_chrom(get_value(cols, extra_snp_chrom_col))
                    pos = int(get_value(cols, extra_snp_pos_col))
                except IndexError:
                    warn("Skipping line %s: not enough columns" % line.strip())
                    continue

                if chr not in extended_range:
                    continue
                if pos < extended_range[chr][0]:
                    continue
                if pos > extended_range[chr][1]:
                    continue

                if snp not in snp_to_chr_pos:
                    snp_to_chr_pos[snp] = (chr, pos)
                    if chr not in chr_to_snp_pos:
                        chr_to_snp_pos[chr] = set()
                    chr_to_snp_pos[chr].add(pos)
                    if (chr, pos) not in chr_pos_to_snp:
                        chr_pos_to_snp[(chr, pos)] = snp

                    if snp not in snp_to_component:
                        component += 1
                        snp_to_component[snp] = component
                        component_to_snp[component] = set()
                        component_to_snp[component].add(snp)

    for chrom in chr_to_snp_pos:
        chr_to_snp_pos[chrom] = sorted(list(chr_to_snp_pos[chrom]))


def set_correlation(snp_to_index, ld_dict):
    cor = np.identity(len(snp_to_index), dtype=np.float64)
    index_list = list(snp_to_index)
    for i in range(len(index_list)):
        for j in range(i):
            snp_1 = index_list[i]
            snp_2 = index_list[j]
            if (snp_1, snp_2) in ld_dict:
                cor[snp_to_index[snp_1], snp_to_index[snp_2]] = ld_dict[snp_1, snp_2]
                cor[snp_to_index[snp_2], snp_to_index[snp_1]] = ld_dict[snp_1, snp_2]
            elif (snp_2, snp_1) in ld_dict:
                cor[snp_to_index[snp_1], snp_to_index[snp_2]] = ld_dict[snp_2, snp_1]
                cor[snp_to_index[snp_2], snp_to_index[snp_1]] = ld_dict[snp_2, snp_1]
    return cor


def set_corr_matrix(snp_to_index, cor):
    global ld_header
    log("assigning correlation values to the correlation matrix of each component..")
    with open_gz(options.ld_file) as ld_fh:
        ld_header = ld_fh.readline().strip().split()
        ld_id1_col = get_col(ld_header, options.ld_id1_col)
        ld_id2_col = get_col(ld_header, options.ld_id2_col)
        ld_r_col = get_col(ld_header, options.ld_r_col)
        for line in ld_fh:
            try:
                cols = line.strip().split()
                value = float(get_value(cols, ld_r_col))
                if abs(value) < options.min_ld_threshold:
                    continue

                snp_1 = get_value(cols, ld_id1_col)
                snp_2 = get_value(cols, ld_id2_col)

            except IndexError:
                warn("Skipping line %s: not enough columns" % line.strip())

            if snp_1 in snp_to_index and snp_2 in snp_to_index:
                if options.debug_make_positive:
                    value = abs(value)
                cor[snp_to_index[snp_1], snp_to_index[snp_2]] = value
                cor[snp_to_index[snp_2], snp_to_index[snp_1]] = value


def gene_simulation(my_range):
    size = options.SNPs_in_simulated_gene
    chrom = None
    range_start = None
    range_end = None
    for key in my_range:
        chrom = key
        range_start = my_range[key][0]
        range_end = my_range[key][1]
    gene = 'Test_gene'
    mid = range_start + ((range_end - range_start) / 2)
    start = mid - 1
    end = mid + size
    gene_to_log_bf = {gene: 0}  # list of genes inside the range to which bf will be calculated
    gene_extend_to_chrom_start_stop = {gene: (chrom, start, end)}
    return gene_to_log_bf, gene_extend_to_chrom_start_stop


def read_gene_file(chrom_data_ranges):
    gene_to_log_bf = {}  # list of genes inside the range to which bf will be calculated
    gene_extend_to_chrom_start_stop = {}
    log("Reading gene file...", DEBUG)
    with open_gz(options.gene_file) as gene_fh:
        gene_header = gene_fh.readline().strip().split()
        gene_id_col = get_col(gene_header, options.gene_id_col)
        gene_chrom_col = get_col(gene_header, options.gene_chrom_col)
        gene_start_col = get_col(gene_header, options.gene_start_col)
        gene_end_col = get_col(gene_header, options.gene_end_col)
        for line in gene_fh:
            cols = line.strip().split()
            try:
                gene = get_value(cols, gene_id_col)
                chrom = clean_chrom(get_value(cols, gene_chrom_col))
                start = int(get_value(cols, gene_start_col))
                end = int(get_value(cols, gene_end_col))
            except IndexError:
                warn("Skipping line %s: not enough columns" % line.strip())

            # ignore genes that can't contribute any causal variants
            # FIXME: is also important to review if the gene is positive or negative
            if chrom in chrom_data_ranges:
                if (chrom_data_ranges[chrom][0] - options.extension) < start and \
                        (chrom_data_ranges[chrom][1] + options.extension) >= end:
                    gene_extend_to_chrom_start_stop[gene] = (chrom, start, end)
                if chrom_data_ranges[chrom][0] < start and chrom_data_ranges[chrom][1] >= end:
                    gene_to_log_bf[gene] = 0
    return gene_to_log_bf, gene_extend_to_chrom_start_stop


# Calculates the SNP link to a gene using table, gaussian or window
def get_snp_prob(f_gene, f_snp, snp_to_chr_pos, gene_to_chrom_start_stop):
    f_snp_chrom, f_snp_pos = snp_to_chr_pos[f_snp]
    f_causal_chrom, f_causal_start, f_causal_end = gene_to_chrom_start_stop[f_gene]
    # The other snp to gene functions were eliminated
    f_causal_start -= options.causal_window
    f_causal_end += options.causal_window
    if f_snp_chrom == f_causal_chrom and f_causal_end > f_snp_pos > f_causal_start:
        f_snp_prob = options.sigma
    else:
        f_snp_prob = options.null_sigma
    return f_snp_prob


def get_snp_dist(f_gene, f_snp, snp_to_chr_pos, gene_to_chrom_start_stop):
    f_snp_chrom, f_snp_pos = snp_to_chr_pos[f_snp]
    f_gene_chrom, f_gene_start, f_gene_end = gene_to_chrom_start_stop[f_gene]
    if f_snp_chrom == f_gene_chrom:
        if f_gene_start <= f_snp_pos <= f_gene_end:
            return 0
        else:
            if f_snp_pos < f_gene_start:
                return abs(f_gene_start - f_snp_pos)
            else:
                return abs(f_snp_pos - f_gene_start)
    else:
        return float('inf')


def close_files():
    global warnings_fh
    if warnings_fh is not None:
        warnings_fh.close()

    global diagnostics_fh
    if diagnostics_fh is not None:
        diagnostics_fh.close()
    global file_ldpred
    if file_ldpred is not None:
        file_ldpred.close()


def select_highest_hit(sigmas, snp_zs, has_z):
    new_sigmas = np.zeros(len(sigmas)) + 0.001
    has_z_pos = np.zeros(len(snp_zs), dtype=int)
    j = 0
    for i in range(len(sigmas)):
        if has_z[i]:
            has_z_pos[j] = i
            j += 1
    max_snp_z_score = -np.inf
    max_snp_z_score_pos = -1
    for i in range(len(has_z_pos)):
        if sigmas[has_z_pos[i]] == 1:
            if abs(snp_zs[i]) > max_snp_z_score:
                max_snp_z_score = snp_zs[i]
                max_snp_z_score_pos = has_z_pos[i]
    if max_snp_z_score_pos != -1:
        new_sigmas[max_snp_z_score_pos] = 1
    return new_sigmas


def process_gene_set_true_beta(gene_set, component_snps_obs, gene_to_snp_probs, cor_matrix_obs, snp_zs):
    mean = np.zeros(len(snp_zs))
    if len(gene_set) == 0:
        sigmas = np.ones(len(snp_zs)) * options.null_sigma
        sigma_var_matrix = np.diag(sigmas)
        f_loglik = scipy.stats.multivariate_normal.logpdf(snp_zs, mean=mean, cov=sigma_var_matrix)
    else:
        sigmas = np.ones(len(snp_zs)) * options.null_sigma
        for i_gene in gene_set:
            if i_gene not in gene_to_snp_probs:
                continue
            gene_sigma = gene_to_snp_probs[i_gene]
            if options.additive:
                sigmas = np.sum([sigmas, gene_sigma], axis=0)
            else:
                sigmas = np.maximum(sigmas, gene_sigma)  # Element-wise maximum of array elements
        sigma_var_matrix = np.diag(sigmas)
        f_loglik = scipy.stats.multivariate_normal.logpdf(snp_zs, mean=mean, cov=sigma_var_matrix)
    hits = np.asarray([1 if x == options.sigma else 0 for x in sigmas])
    return gene_set, f_loglik, hits.astype(int)


def process_gene_set(gene_set, component_snps_obs, gene_to_snp_probs, cor_matrix_obs, snp_zs):
    if options.true_beta:
        return process_gene_set_true_beta(gene_set, component_snps_obs, gene_to_snp_probs, cor_matrix_obs, snp_zs)
    else:
        return process_gene_set_z(gene_set, component_snps_obs, gene_to_snp_probs, cor_matrix_obs, snp_zs)


# Calculates the likelihood
def process_gene_set_z(gene_set, component_snps_obs, gene_to_snp_probs, cor_matrix_obs, snp_zs):
    frac = options.p
    mean = np.zeros(len(snp_zs))
    if len(gene_set) == 0:
        sigmas = np.ones(len(snp_zs)) * options.null_sigma
        if options.add_cov_uncertainty:
            f_loglik = ln_prob_z(snp_zs, sigmas, cor_matrix_obs)
        else:
            sigma_var_matrix = np.diag(sigmas)
            cov_matrix = cor_matrix_obs + np.dot(cor_matrix_obs, np.dot(sigma_var_matrix, cor_matrix_obs))
            f_loglik = scipy.stats.multivariate_normal.logpdf(snp_zs, mean=mean, cov=cov_matrix)
        hits = np.asarray([1 if x == options.sigma else 0 for x in sigmas])
    else:
        sigmas = np.zeros(len(component_snps_obs))
        for i_gene in gene_set:
            if i_gene not in gene_to_snp_probs:
                continue
            gene_sigma = gene_to_snp_probs[i_gene]
            if options.additive:
                sigmas = np.sum([sigmas, gene_sigma], axis=0)
            else:
                sigmas = np.maximum(sigmas, gene_sigma)  # Element-wise maximum of array elements
        hits = np.asarray([1 if x == options.sigma else 0 for x in sigmas])
        hit_num = sum(hits)
        aux = 0
        idx = np.zeros(hit_num)
        for j in range(len(component_snps_obs)):
            if hits[j] == 1:
                idx[aux] = j
                aux += 1
        idx = idx.astype(int)
        if options.SNP_model == "All_one":
            if options.add_cov_uncertainty:
                f_loglik = ln_prob_z(snp_zs, sigmas, cor_matrix_obs)
            else:
                sigma_var_matrix = np.diag(sigmas)
                cov_matrix = cor_matrix_obs + np.dot(cor_matrix_obs, np.dot(sigma_var_matrix, cor_matrix_obs))
                f_loglik = scipy.stats.multivariate_normal.logpdf(snp_zs, mean=mean, cov=cov_matrix)
        elif options.SNP_model == "Full_SNP":
            skip = True
            # calculate the likelihood of each SNP configuration
            ln_lik_com = np.zeros((2 ** hit_num))
            j = 0
            for combi in product(range(2), repeat=hit_num):
                if skip:
                    skip = False
                    continue
                c = np.zeros(len(component_snps_obs))
                initial_sigma_vec = np.ones(len(component_snps_obs)) * options.null_sigma
                for i in range(hit_num):
                    c[idx[i]] = combi[i]
                    initial_sigma_vec[idx[i]] = (combi[i] * options.sigma) + ((1 - combi[i]) * options.null_sigma)
                sigma_var_matrix = np.diag(initial_sigma_vec)
                causal_snps = sum(c)
                cov_matrix = cor_matrix_obs + np.dot(cor_matrix_obs, np.dot(sigma_var_matrix, cor_matrix_obs))
                cov_matrix = matrix_stabilization(cov_matrix)
                factor = (causal_snps * math.log(frac)) + ((hit_num - causal_snps) * math.log(1 - frac))
                if options.add_cov_uncertainty:
                    f_loglik = ln_prob_z(snp_zs, np.multiply(sigmas, c), cor_matrix_obs)
                else:
                    f_loglik = scipy.stats.multivariate_normal.logpdf(snp_zs, mean=mean, cov=cov_matrix)
                ln_lik_com[j] = factor + f_loglik
                j += 1
            max_ln = max(ln_lik_com)
            f_loglik = max_ln + np.log(sum([np.exp(x - max_ln) for x in ln_lik_com]))
        else:
            f_loglik = greedy_SNP_configuration(hits, snp_zs, cor_matrix_obs, sigmas, frac)
    return gene_set, f_loglik, hits.astype(int)


def greedy_SNP_configuration(hits, snp_zs, cor_matrix_obs, sigmas, frac):
    global num_lik_calc
    global iterations
    mean = np.zeros(len(snp_zs))
    ln_lik_vec = []
    keep_conf = []
    my_idx = set()
    num_hits = sum(hits)
    for i in range(len(hits)):
        if hits[i] == 1:
            my_idx.add(i)
            # keep_conf.append(set(tuple([i])))
    # Calculation of the all zero configuration
    sigma_var_matrix = np.diag(np.ones(len(hits)) * options.null_sigma)
    cov_matrix = cor_matrix_obs + np.dot(cor_matrix_obs, np.dot(sigma_var_matrix, cor_matrix_obs))
    factor = num_hits * math.log(1 - frac)
    if options.add_cov_uncertainty:
        ln_lik = factor + ln_prob_z(snp_zs, np.ones(len(hits)) * options.null_sigma, cor_matrix_obs)
    else:
        ln_lik = factor + scipy.stats.multivariate_normal.logpdf(snp_zs, mean=mean, cov=cov_matrix)
    ln_lik_vec.append(ln_lik)
    # ------------------------------------------
    all_conf = [set()]
    aux = [set()]
    for i in range(num_hits):  # greedy approach
        my_max = -np.inf
        for conf in aux:
            initial_causal_vec = np.ones(len(hits)) * options.null_sigma
            for element in conf:
                initial_causal_vec[element] = options.sigma
            for val in (my_idx - conf):
                new_conf = conf.union([val])
                iterations += 1
                if new_conf not in all_conf:
                    causal_vec = initial_causal_vec.copy()
                    causal_vec[val] = options.sigma
                    sigma_var_matrix = np.diag(causal_vec)
                    cov_matrix = cor_matrix_obs + np.dot(cor_matrix_obs, np.dot(sigma_var_matrix, cor_matrix_obs))
                    causal_snps = len(new_conf)
                    factor = (causal_snps * math.log(frac)) + ((num_hits - causal_snps) * math.log(1 - frac))
                    if options.add_cov_uncertainty:
                        ln_lik = factor + ln_prob_z(snp_zs, causal_vec, cor_matrix_obs)
                    else:
                        ln_lik = factor + scipy.stats.multivariate_normal.logpdf(snp_zs, mean=mean, cov=cov_matrix)
                    ln_lik_vec.append(ln_lik)
                    num_lik_calc += 1
                    if ln_lik >= my_max - options.SNP_log_diff_ignore:
                        keep_conf.append(new_conf)
                        all_conf.append(new_conf)
                        if ln_lik > my_max:
                            my_max = ln_lik
        aux = list(keep_conf)
        keep_conf = []
    max_ln = max(ln_lik_vec)
    f_loglik = max_ln + np.log(sum([np.exp(x - max_ln) for x in ln_lik_vec]))
    return f_loglik


def calculate_likelihood_greedy(genes_in_component, component_snps_obs,
                                gene_to_snp_probs, cor_matrix_obs, snp_zs, gene_to_log_bf):
    log("Genes in the component: %s" % genes_in_component, DEBUG)
    gene_to_alt = {'sigmas': {}, 'lik': {}, 'tot': {}, 'gene_sets': {}}
    gene_to_null = {'sigmas': {}, 'lik': {}, 'tot': {}, 'gene_sets': {}}

    max_overall_lik = -np.inf

    keep_gene_sets = set()
    gene_set = tuple("")
    for num_causal_genes in range(min(options.max_num_causal_genes, len(genes_in_component)) + 1):
        gene_set_logliks = []
        # calculate likelihood of the selected gene sets
        if num_causal_genes == 0:
            results = process_gene_set(gene_set, component_snps_obs, gene_to_snp_probs, cor_matrix_obs, snp_zs)
            gene_set_logliks.append(results)
        else:
            prev_gene_sets = keep_gene_sets
            keep_gene_sets = set()
            processed_gene_sets = set()
            for new_gene in genes_in_component:
                for prev_gene_set in prev_gene_sets:
                    if new_gene in prev_gene_set:
                        continue
                    gene_set = tuple(sorted((new_gene,) + prev_gene_set))
                    if gene_set in processed_gene_sets:
                        continue
                    processed_gene_sets.add(gene_set)
                    results = process_gene_set(gene_set, component_snps_obs, gene_to_snp_probs, cor_matrix_obs,
                                               snp_zs)
                    gene_set_logliks.append(results)

        # selecting the highest log lik
        for gene_set_loglik in sorted(gene_set_logliks, key=lambda x: x[1], reverse=True):
            gene_set = gene_set_loglik[0]
            loglik = gene_set_loglik[1]
            sigmas = gene_set_loglik[2]
            prob_set = (options.gene_prior ** len(gene_set)) * (1 - options.gene_prior) ** (
                    len(genes_in_component) - len(gene_set))
            adjusted_loglik = loglik + np.log(prob_set)
            # if gene_set != tuple():
            #     log("%s" % sum(sigmas.astype(int)), DEBUG, end_char='\t')
            log("Using %s, Causal SNPs %s, adjusted_loglik: %s, loglik: %s, prob_set: %s" % (
                str(gene_set), sum(sigmas.astype(int)), adjusted_loglik, loglik, prob_set), DEBUG)
            # keep for next time only if (a) adjusted likelihood within log_magnitude of best we've seen so far
            # and (b) if this gene set has better likelihood than what we saw last time
            # (did we even marginally benefit by adding this gene set)

            if adjusted_loglik >= max_overall_lik - options.log_diff_ignore:
                keep_gene_sets.add(gene_set)  # gene set that stay
                if adjusted_loglik > max_overall_lik:
                    max_overall_lik = adjusted_loglik  # setting the new max

        # Creation of the set for null and alternative hypothesis
        for gene_set_loglik in sorted(gene_set_logliks, key=lambda x: x[1], reverse=True):
            gene_set = gene_set_loglik[0]
            loglik = gene_set_loglik[1]
            sigmas = gene_set_loglik[2]
            prob_set = (options.gene_prior ** len(gene_set)) * (1 - options.gene_prior) ** (
                    len(genes_in_component) - len(gene_set))
            adjusted_loglik = loglik + np.log(prob_set)

            for gene in genes_in_component:
                if gene not in gene_to_log_bf:
                    continue
                if gene in gene_set:
                    add_gene(gene, gene_to_alt, adjusted_loglik, prob_set, gene_set, sigmas)

                else:
                    add_gene(gene, gene_to_null, adjusted_loglik, prob_set, gene_set, sigmas)

    return gene_to_alt, gene_to_null


def powerset(iterable):
    "powerset([1,2,3]) --> () (1,) (2,) (3,) (1,2) (1,3) (2,3) (1,2,3)"
    s = list(iterable)
    return list(chain.from_iterable(combinations(s, r) for r in range(len(s) + 1)))


def add_gene(gene, hypothesis, adjusted_loglik, prob_set, gene_set, sigmas):
    if gene not in hypothesis['lik']:
        hypothesis['lik'][gene] = []
        hypothesis['tot'][gene] = []
        hypothesis['gene_sets'][gene] = []
        hypothesis['sigmas'][gene] = []
    hypothesis['lik'][gene].append(adjusted_loglik)
    hypothesis['tot'][gene].append(prob_set)
    hypothesis['gene_sets'][gene].append(gene_set)
    hypothesis['sigmas'][gene].append(sigmas)


def calculate_null_prior(gene, genes_in_analysis):
    power_set = powerset(genes_in_analysis)
    prior = 0
    for gene_set in power_set:
        if gene not in gene_set:
            prior += (options.gene_prior ** len(gene_set)) * (1 - options.gene_prior) ** (
                    len(genes_in_analysis) - len(gene_set))
    return prior


def calculate_likelihood_full(genes_in_component, component_snps_obs,
                              gene_to_snp_probs, cor_matrix_obs, snp_zs, gene_to_log_bf):
    gene_to_alt = {'sigmas': {}, 'lik': {}, 'tot': {}, 'gene_sets': {}}
    gene_to_null = {'sigmas': {}, 'lik': {}, 'tot': {}, 'gene_sets': {}}
    power_set = powerset(genes_in_component)
    for gene_set in power_set:
        results = process_gene_set(gene_set, component_snps_obs, gene_to_snp_probs, cor_matrix_obs, snp_zs)
        loglik = results[1]
        sigmas = results[2]
        prob_set = (options.gene_prior ** len(gene_set)) * (1 - options.gene_prior) ** (
                len(genes_in_component) - len(gene_set))
        adjusted_loglik = loglik + np.log(prob_set)
        log("Process %s, Causal SNPs %s, adjusted_loglik: %s, loglik: %s, prob_set: %s" % (
            str(gene_set), sum(sigmas.astype(int)), adjusted_loglik, loglik, prob_set), DEBUG)
        for gene in genes_in_component:
            if gene not in gene_to_log_bf:
                continue
            if gene in gene_set:
                add_gene(gene, gene_to_alt, adjusted_loglik, prob_set, gene_set, sigmas)
            else:
                add_gene(gene, gene_to_null, adjusted_loglik, prob_set, gene_set, sigmas)

    return gene_to_alt, gene_to_null


def calculate_likelihood_single_gene(genes_in_component, component_to_snp, component, gene_to_snp_probs, gene_to_sigma,
                                     snp_ns, null_cov_matrix, cor_matrix, has_z, snp_zs, gene_to_log_bf):
    gene_to_alt = {'sigmas': {}, 'lik': {}, 'tot': {}, 'gene_sets': {}}
    gene_to_null = {'sigmas': {}, 'lik': {}, 'tot': {}, 'gene_sets': {}}

    empty_set = tuple()
    results = process_gene_set(empty_set, component_to_snp, component, gene_to_snp_probs, gene_to_sigma,
                               snp_ns, null_cov_matrix, cor_matrix, has_z, snp_zs)

    loglik = results[1]
    sigmas_empty = results[2]
    prob_set_empty = (options.gene_prior ** len(empty_set)) * (1 - options.gene_prior) ** (
            len(genes_in_component) - len(empty_set))
    adjusted_loglik_empty = loglik + np.log(prob_set_empty)

    for gene in genes_in_component:
        results = process_gene_set((gene,), component_to_snp, component, gene_to_snp_probs, gene_to_sigma,
                                   snp_ns, null_cov_matrix, cor_matrix, has_z, snp_zs)

        loglik = results[1]
        sigmas = results[2]
        prob_set = (options.gene_prior ** len(gene)) * (1 - options.gene_prior) ** (
                len(genes_in_component) - len(gene))
        adjusted_loglik = loglik + np.log(prob_set)

        add_gene(gene, gene_to_null, adjusted_loglik_empty, prob_set_empty, empty_set, sigmas_empty)
        add_gene(gene, gene_to_alt, adjusted_loglik, prob_set, gene, sigmas)

    return gene_to_alt, gene_to_null


def matrix_stabilization(matrix):
    # stabilization achieved by multiply the diagonal by a small factor to avoid singular matrix
    matrix_size = matrix.shape[0]  # this is a square matrix
    null_cov_diag = np.diag(matrix)
    stab_matrix = np.diag(np.ones(matrix_size)) * (np.mean(np.diag(matrix)) * options.stab_term_frac)
    null_cov_matrix = matrix + stab_matrix

    # eigen value decomposition to assure the resulting correlation matrix is positive semi-definite
    try:
        l, Q = np.linalg.eigh(null_cov_matrix)
    except np.linalg.linalg.LinAlgError:
        log("")
        raise

    l[l < 0] = 0
    null_cov_matrix = np.dot(Q, np.dot(np.diag(l), Q.transpose())) + stab_matrix

    # transformation to restore the values changed in the stabilization step
    null_cov_norm_matrix = np.diag(np.sqrt(null_cov_diag / np.diag(null_cov_matrix)))
    null_cov_matrix = np.dot(null_cov_norm_matrix, np.dot(null_cov_matrix, null_cov_norm_matrix))
    return null_cov_matrix


'''
def plot_data(chr_pos_to_snp, snp_to_chr_pos, snp_to_p, gene_to_chrom_start_stop, snp_to_analyze):
    snps = [snp for snp in snp_to_chr_pos if snp in snp_to_analyze]  # set with SNPs with position and z-score
    pos = [snp_to_chr_pos[snp] for snp in snps]
    sorted_keys = sorted(pos)
    x = [i[1] for i in sorted_keys]
    y = [snp_to_p[chr_pos_to_snp[j]] for j in sorted_keys]
    y = -np.log10(y)
    plt.scatter(x, np.abs(y), color='#88c999')
    x = [gene_to_chrom_start_stop[gene][1] for gene in gene_to_chrom_start_stop]
    y = [x for x in range(len(x))]
    lab = [gene for gene in gene_to_chrom_start_stop]
    for xs, ys in zip(x, y):
        plt.scatter(xs, (ys - len(lab)), color='hotpink')
        plt.annotate(lab[ys], (xs, (ys - len(lab))))
    plt.legend()
    plt.grid(False)
    f = plt.figure()
    f.set_figwidth(40)
    f.set_figheight(10)
    plt.show()
    log("plot")
'''


def hit_count(gene, gene_to_alt_num_snps, gene_to_alt, gene_to_null_num_snps, gene_to_null):
    if gene in gene_to_alt_num_snps:
        for sigma in gene_to_alt['sigmas'][gene]:  # for each gene set in the alternative hypothesis
            gene_to_alt_num_snps[gene] += sum(sigma)
    else:
        gene_to_alt_num_snps[gene] = 0
        for sigma in gene_to_alt['sigmas'][gene]:  # for each gene set in the alternative hypothesis
            gene_to_alt_num_snps[gene] += sum(sigma)

    if gene in gene_to_null_num_snps:
        for sigma in gene_to_null['sigmas'][gene]:  # for each gene set in the alternative hypothesis
            gene_to_null_num_snps[gene] += sum(sigma)
    else:
        gene_to_null_num_snps[gene] = 0
        for sigma in gene_to_null['sigmas'][gene]:  # for each gene set in the alternative hypothesis
            gene_to_null_num_snps[gene] += sum(sigma)


def select_highest_component(component_to_cor):
    max_value = 0
    name = -1
    for component in component_to_cor:
        if component_to_cor[component].shape[0] > max_value:
            max_component = component_to_cor[component]
            max_value = component_to_cor[component].shape[0]
            name = component
    component_to_return = {name: max_component}
    return component_to_return


def select_smallest_component(component_to_cor):
    min_value = np.inf
    name = -1
    for component in component_to_cor:
        if component_to_cor[component].shape[0] < min_value:
            min_component = component_to_cor[component]
            min_value = component_to_cor[component].shape[0]
            name = component
    component_to_return = {name: min_component}
    return component_to_return


def save_data_for_FINEMAP_ld(cor_matrix_obs):
    np.savetxt("chrom2.ld", cor_matrix_obs, fmt='%.2f')


def save_data_for_COJO_ma(snp_to_effect_allele, snp_to_n, snp_to_p, component_snps_obs,
                          snp_to_other_allele, snp_to_beta, snp_to_se, snp_to_frq):
    file = open('LRPPRC.ma', 'w')
    file.write("SNP\tA1\tA2\tfreq\tbeta\tse\tp\tN")
    # making sure the order in the cor matrix is the same as in the z file

    for snp in component_snps_obs:  # add only SNPs with LD and z-score
        if snp in snp_to_frq:
            allele1 = snp_to_effect_allele[snp]
            allele2 = snp_to_other_allele[snp]
            maf = str(snp_to_frq[snp])
            beta = str(snp_to_beta[snp])
            se = str(snp_to_se[snp])
            p = str(scipy.stats.norm.sf(abs(float(beta) / float(se))) * 2)
            n = str(snp_to_n[snp])
            line = snp + "\t" + allele1 + "\t" + allele2 + "\t" + maf + "\t" + beta + "\t" + se + "\t" + p + "\t" + n
            file.write("\n" + line)
    log("File created Cojo format")
    file.close()


def read_bim_file():
    log("Reading bim file...")
    file = "../data/g1000_eur.bim"
    snp_to_chr_pos = {}
    with open_gz(file) as sumstats_fh:
        for line in sumstats_fh:
            cols = line.strip().split()
            try:
                snp = get_value(cols, 1)
                chr = int(get_value(cols, 0))
                pos = int(get_value(cols, 3))
            except IndexError:
                warn("Skipping line %s: not enough columns" % line.strip())
            snp_to_chr_pos[snp] = (chr, pos)
    return snp_to_chr_pos


def save_data_for_magma(snp_to_effect_allele, snp_to_n, snp_to_p,
                        snp_to_other_allele, snp_to_beta, snp_to_se, snp_to_analyze):
    name = options.magma_out  # + "_{}.magma".format(options.range.split(':')[0])
    file = open(name, 'w')
    if options.chrom == "1":
        file.write("SNP\tA1\tA2\tb\tse\tP\tN")
    for snp in snp_to_analyze:  # the snp needs to be in the 3 dataset, gwas, frequency and position
        allele1 = str(snp_to_other_allele[snp])
        allele2 = str(snp_to_effect_allele[snp])
        beta = str(snp_to_beta[snp])
        se = str(snp_to_se[snp])
        p = str(snp_to_p[snp])
        n = str(snp_to_n[snp])
        line = snp + "\t" + allele1 + "\t" + allele2 + "\t" + beta + "\t" + se + "\t" + p + "\t" + n
        file.write("\n" + line)
    file.close()


def save_data_for_locus_gwas(snp_to_effect_allele, snp_to_n, snp_to_p,
                             snp_to_other_allele, snp_to_beta, snp_to_se, snp_to_chr_pos):
    name = options.locus_out  # + "_{}.locus".format(options.range.split(':')[0])
    file = open(name, 'w')
    if options.chrom == "1":
        file.write("SNP\tChr\tbp\tA1\tA2\tfreq\tb\tse\tp\tN")
    # making sure the order in the cor matrix is the same as in the z file
    ordered_snp = sorted(snp_to_chr_pos.items(), key=lambda x: (x[1], x[0]), reverse=False)
    for snp, value in ordered_snp:  # the snp needs to be in the 3 dataset, gwas, frequency and position
        if snp in snp_to_chr_pos:
            chrom = str(snp_to_chr_pos[snp][0])
            pos = str(snp_to_chr_pos[snp][1])
            allele1 = str(snp_to_other_allele[snp])
            allele2 = str(snp_to_effect_allele[snp])
            maf = "0.1"
            beta = str(snp_to_beta[snp])
            se = str(snp_to_se[snp])
            p = str(snp_to_p[snp])
            n = str(snp_to_n[snp])
            line = snp + "\t" + chrom + "\t" + pos + "\t" + allele1 + "\t" + allele2 + "\t" + maf + "\t" + beta + "\t" + se + "\t" + p + "\t" + n
            file.write("\n" + line)

    file.close()


def save_data_for_FINEMAP_z(component, component_to_snp, snp_to_chr_pos, snp_to_effect_allele,
                            snp_to_other_allele, snp_to_beta, snp_to_se, snp_to_index):
    file = open('chrom2.z', 'w')
    file.write("rsid chromosome position allele1 allele2 maf beta se")
    # making sure the order in the cor matrix is the same as in the z file
    index_to_snp = {}
    for snp in component_to_snp[component]:
        index_to_snp[snp_to_index[snp]] = snp
    ordered_snp = dict(sorted(index_to_snp.items()))
    for snp in ordered_snp.values():  # not all snp in component have z-score, allele, etc.
        if snp in snp_to_effect_allele:
            chrom = str(snp_to_chr_pos[snp][0])
            pos = str(snp_to_chr_pos[snp][1])
            allele1 = snp_to_effect_allele[snp]
            allele2 = snp_to_other_allele[snp]
            maf = "0.03"  # FIXME: find a dataset with maf
            beta = str(snp_to_beta[snp])
            se = str(snp_to_se[snp])
            line = snp + " " + chrom + " " + pos + " " + allele1 + " " + allele2 + " " + maf + " " + beta + " " + se
            file.write("\n" + line)

    file.close()


def filter_snps(snp_to_chr_pos, snp_to_p_dentist):
    snp_to_filter = {}
    for snp in snp_to_chr_pos:
        if snp in snp_to_p_dentist:
            snp_to_filter[snp] = 1
    return snp_to_filter


def print_results(gene_to_log_bf, gene_to_null_num_snps, gene_to_alt_num_snps):
    genes_in_analysis = sorted(gene_to_log_bf, key=lambda k: gene_to_log_bf[k], reverse=True)
    for gene in genes_in_analysis:
        pi = 1 - options.gene_prior  # calculate_null_prior(gene, genes_in_analysis)  # Prior probability of the null
        # log("Prior probability that H_1 is true given that one of H_0 or H_1 is true: %s" % (1-pi))
        po = pi / (1 - pi)  # posterior odds
        bf = np.exp(gene_to_log_bf[gene])
        p_null = (bf * po) / ((bf * po) + 1)
        # log("%s \tPosterior probability of H_1: %s \tCausal SNPs alt: %d \tCausal SNPs null: "
        #      "%d \tLogBF: %.10g\tBF: %.10g" % (gene, 1 - p_null, gene_to_alt_num_snps[gene],
        #                                        gene_to_null_num_snps[gene], gene_to_log_bf[gene], bf))
        # log("%s \t%s \t%s" % (gene, 1 - p_null, gene_to_log_bf[gene]))
        log("%s \t%s \t%s \t%s" % (gene, 1 - p_null, gene_to_log_bf[gene], gene_to_alt_num_snps[gene]))
        # log("%s \t%s \t%s" % (gene, 1, gene_to_log_bf[gene]))


def print_bf_calculation(gene, gene_to_alt, gene_to_null, null, alt):
    log("Gene: %s" % gene, DEBUG)

    sorted_alt_indices = sorted(range(len(gene_to_alt['lik'][gene])), key=lambda i: gene_to_alt['lik'][gene][i],
                                reverse=True)
    sorted_null_indices = sorted(range(len(gene_to_null['lik'][gene])), key=lambda i: gene_to_null['lik'][gene][i],
                                 reverse=True)

    sorted_alt_lik = sorted(gene_to_alt['lik'][gene], reverse=True)
    sorted_null_lik = sorted(gene_to_null['lik'][gene], reverse=True)

    log("Alt adjusted_loglik: %s, %s" % (len(sorted_alt_lik), sorted_alt_lik[:10]), DEBUG)
    log("Alt Gene Sets: %s, %s" % (
        len(gene_to_alt['gene_sets'][gene]), [gene_to_alt['gene_sets'][gene][i] for i in sorted_alt_indices]), DEBUG)
    log("Alt Prob: %s" % [gene_to_alt['tot'][gene][i] for i in sorted_alt_indices[:10]], DEBUG)
    log("Null adjusted_loglik: %s, %s" % (len(sorted_null_lik), sorted_null_lik[:10]), DEBUG)
    log("Null Gene Sets: %s, %s" % (
        len(gene_to_null['gene_sets'][gene]), [gene_to_null['gene_sets'][gene][i] for i in sorted_null_indices]), DEBUG)
    log("Alt Prob: %s" % [gene_to_null['tot'][gene][i] for i in sorted_null_indices[:10]], DEBUG)

    alt_sigmas = [sum(gene_to_alt['sigmas'][gene][i]) for i in
                  sorted(range(len(gene_to_alt['sigmas'][gene])), reverse=True)]
    null_sigmas = [sum(gene_to_null['sigmas'][gene][i]) for i in
                   sorted(range(len(gene_to_null['sigmas'][gene])), reverse=True)]

    log("Total causal SNPs in Null Gene Sets: %s, Causal SNPs per set %s" % (sum(null_sigmas), null_sigmas), DEBUG)
    log("Total causal SNPs in Alt  Gene Sets: %s, Causal SNPs per set %s" % (sum(alt_sigmas), alt_sigmas), DEBUG)
    log("Null - Alt: %s - %s" % (null, alt), DEBUG)
    log("logBF Increment: %s" % (null - alt), DEBUG)
    # log("%s" % (null - alt), INFO)


def bf_calculation(gene_to_log_bf, gene_to_null_h, gene_to_alt_h, gene_to_alt_num_snps, gene_to_null_num_snps,
                   gene_to_alt, gene_to_null):
    for gene in gene_to_alt['lik']:
        max_log_alt = max(gene_to_alt['lik'][gene])
        max_log_null = max(gene_to_null['lik'][gene])
        log_sum_alt = max_log_alt + np.log(sum([np.exp(x - max_log_alt) for x in gene_to_alt['lik'][gene]]))
        log_sum_null = max_log_null + np.log(sum([np.exp(x - max_log_null) for x in gene_to_null['lik'][gene]]))

        alt = log_sum_alt - np.log(sum(gene_to_alt['tot'][gene]))
        null = log_sum_null - np.log(sum(gene_to_null['tot'][gene]))

        gene_to_log_bf[gene] += null - alt
        gene_to_null_h[gene] += null
        gene_to_alt_h[gene] += alt

        # counting the number of hits inside genes in the gene set for each hypothesis
        hit_count(gene, gene_to_alt_num_snps, gene_to_alt, gene_to_null_num_snps, gene_to_null)

        print_bf_calculation(gene, gene_to_alt, gene_to_null, null, alt)


def print_collected_statistics(snp_to_index, component_snps_obs, snp_to_chr_pos, snp_zs, snp_ns):
    log("SNPs inside the region with LD and z-score in the component: %s" % len(component_snps_obs), DEBUG)
    # log("%s" % len(component_snps_obs), INFO, end_char='\t')
    log("SNPs: %s" % sorted(component_snps_obs, key=lambda x: snp_to_chr_pos[x][1]), DEBUG)
    log("Positions: %s" % [snp_to_chr_pos[x][1] for x in
                           sorted(component_snps_obs, key=lambda x: snp_to_chr_pos[x][1])], TRACE)
    log("SNP zs: %s" % [snp_zs[i] for i in
                        sorted(range(len(component_snps_obs)), key=lambda i: snp_to_chr_pos[component_snps_obs[i]][1])],
        TRACE)
    log("SNP ns: %s" % [snp_ns[i] for i in
                        sorted(range(len(component_snps_obs)), key=lambda i: snp_to_chr_pos[component_snps_obs[i]][1])],
        TRACE)


def calculate_likelihood(genes_in_component, component_snps_obs, gene_to_snp_probs, cor_matrix_obs, snp_zs,
                         gene_to_log_bf):
    if options.model == "greedy":
        return calculate_likelihood_greedy(genes_in_component, component_snps_obs, gene_to_snp_probs, cor_matrix_obs,
                                           snp_zs, gene_to_log_bf)
    else:
        return calculate_likelihood_full(genes_in_component, component_snps_obs, gene_to_snp_probs, cor_matrix_obs,
                                         snp_zs, gene_to_log_bf)


def distance(snp, genes_in_component, snp_to_chr_pos, gene_to_chrom_start_stop):
    gene_to_distance = {}
    for gene in genes_in_component:
        gene_to_distance[gene] = abs(float(snp_to_chr_pos[snp][1]) - float(gene_to_chrom_start_stop[gene][1]))
    return list(sorted(gene_to_distance.items(), key=lambda item: item[1]))


def gene_to_snp_70(genes_in_component, snp_to_index, snp_to_chr_pos, gene_to_chrom_start_stop, has_z):
    snp_to_gene = {}
    for snp in snp_to_index:
        ordered_genes = distance(snp, genes_in_component, snp_to_chr_pos, gene_to_chrom_start_stop)
        not_assigned = True
        for gene in ordered_genes:
            if abs(np.random.normal(0, 1, 1)) <= 1.036433:
                snp_to_gene[snp] = gene
                not_assigned = False
                break
        if not_assigned:
            snp_to_gene[snp] = ordered_genes[-1]
    gene_to_snp_probs = {}
    for gene in genes_in_component:
        gene_to_snp_probs[gene] = np.zeros(len(snp_to_index))
        for snp in snp_to_index:
            if snp_to_gene[snp][0] == gene:
                gene_to_snp_probs[gene][snp_to_index[snp]] = options.sigma
        gene_to_snp_probs[gene] = gene_to_snp_probs[gene][has_z]
    return gene_to_snp_probs


def gene_to_snp(genes_in_component, snp_to_index, snp_to_chr_pos, gene_to_chrom_start_stop, has_z):
    gene_to_snp_probs = {}
    for gene in genes_in_component:
        log("SNPs in gene %s: " % gene, TRACE)
        gene_to_snp_probs[gene] = np.zeros(len(snp_to_index))
        for snp in snp_to_index:
            snp_prob = get_snp_prob(gene, snp, snp_to_chr_pos, gene_to_chrom_start_stop)
            gene_to_snp_probs[gene][snp_to_index[snp]] = snp_prob
            if snp_prob == options.sigma:
                log("%s, " % snp, TRACE)
        gene_to_snp_probs[gene] = gene_to_snp_probs[gene][has_z]
    return gene_to_snp_probs


def gene_to_snp_distance_function(genes_in_component, snp_to_index, snp_to_chr_pos, gene_to_chrom_start_stop, has_z):
    gene_to_snp_distance = {}
    for gene in genes_in_component:
        log("SNPs in gene %s: " % gene, TRACE)
        gene_to_snp_distance[gene] = np.zeros(len(snp_to_index))
        for snp in snp_to_index:
            snp_dist = get_snp_dist(gene, snp, snp_to_chr_pos, gene_to_chrom_start_stop)
            gene_to_snp_distance[gene][snp_to_index[snp]] = snp_dist
        gene_to_snp_distance[gene] = gene_to_snp_distance[gene][has_z]
    return gene_to_snp_distance


def add_genes_to_component(genes, snp_to_index, snp_to_chr_pos):
    genes_in_component = set()
    for gene in genes:
        causal_chrom, causal_start, causal_end = genes[gene]
        causal_start -= options.causal_window
        causal_end += options.causal_window
        for snp in snp_to_index:
            snp_chrom, snp_pos = snp_to_chr_pos[snp]
            if snp_chrom == causal_chrom and causal_start < snp_pos < causal_end:
                genes_in_component.add(gene)
                break
    return genes_in_component


def create_index(snp_to_analyze):
    snp_to_index = {}
    index = 0
    for snp in sorted(snp_to_analyze):  # snp_to_chr_pos contains all SNPs inside the region with LD
        snp_to_index[snp] = index
        index += 1
    return snp_to_index


def extract_snps(snp_to_index, snp_to_z, snp_to_n, snp_to_beta, snp_to_se):
    snp_betas = np.zeros(len(snp_to_index))
    snp_ses = np.zeros(len(snp_to_index))
    snp_zs = np.zeros(len(snp_to_index))
    snp_ns = np.zeros(len(snp_to_index))
    component_snps_obs = np.array([None] * len(snp_to_index))
    has_z = np.array([False] * len(snp_to_index))
    for snp in snp_to_index:  # SNPs with LD
        if snp in snp_to_z:  # SNPs with z-score
            snp_betas[snp_to_index[snp]] = snp_to_beta[snp]
            snp_ses[snp_to_index[snp]] = snp_to_se[snp]
            snp_zs[snp_to_index[snp]] = snp_to_z[snp]
            snp_ns[snp_to_index[snp]] = snp_to_n[snp]
            component_snps_obs[snp_to_index[snp]] = snp
            has_z[snp_to_index[snp]] = True
    # subset to observed
    snp_betas = snp_betas[has_z]
    snp_ses = snp_ses[has_z]
    snp_zs = snp_zs[has_z]
    snp_ns = snp_ns[has_z]  # will contain all SNPs with LD and z-score inside the given region
    component_snps_obs = component_snps_obs[has_z]  # this will remove all the None rows
    # plt.hist(snp_zs)  # plot of the distribution of the analyzed SNPs
    # plt.show()
    return snp_zs, snp_ns, component_snps_obs, has_z, snp_betas, snp_ses


def create_index_to_snp(snp_to_analyze):
    index_to_snp = {}
    index = 0
    for snp in snp_to_analyze:
        index_to_snp[index] = snp
        index += 1
    return index_to_snp


def matrix_correlation(ld_dict, all_snp):
    index_to_snp = create_index_to_snp(all_snp)
    correlation = np.identity(len(all_snp))
    for i in range(len(all_snp)):
        for j in range(i):
            snp1 = index_to_snp[i]
            snp2 = index_to_snp[j]
            if (snp1, snp2) in ld_dict:
                correlation[i, j] = ld_dict[snp1, snp2]
                correlation[j, i] = ld_dict[snp1, snp2]
            elif (snp2, snp1) in ld_dict:
                correlation[i, j] = ld_dict[snp2, snp1]
                correlation[j, i] = ld_dict[snp2, snp1]
    print("correlation")
    return correlation, index_to_snp


def cluster_corr(corr_array, index_to_snp, inplace=False):
    pairwise_distances = sch.distance.pdist(corr_array)
    linkage = sch.linkage(pairwise_distances, method='complete')
    cluster_distance_threshold = pairwise_distances.max() / 2
    idx_to_cluster_array = sch.fcluster(linkage, cluster_distance_threshold,
                                        criterion='distance')
    idx = np.argsort(idx_to_cluster_array)
    cluster_to_snps = {}
    for i in range(len(idx_to_cluster_array)):
        cluster = idx_to_cluster_array[i]
        if cluster not in cluster_to_snps:
            cluster_to_snps[cluster] = set()
        cluster_to_snps[cluster].add(index_to_snp[i])
    if not inplace:
        corr_array = corr_array.copy()

    return corr_array[idx, :][:, idx], cluster_to_snps


'''
def components_from_cluster(ld_dict, snp_to_analyze):
    correlation, index_to_snp = matrix_correlation(ld_dict, snp_to_analyze)
    plt.figure(1)
    sns.heatmap(correlation)
    plt.show()
    # reordered_corr, cluster_to_snps = cluster_corr(correlation, index_to_snp)
    # plt.figure(2)
    # sns.heatmap(reordered_corr)
    # plt.show()
    # return cluster_to_snps
'''


def add_dentist(snp_to_z):
    if options.dentist_file is not None:
        snp_to_p_dentist = read_dentist(options.dentist_thr)  # SNP with values below the parameter are eliminated
        snps = filter_snps(snp_to_z, snp_to_p_dentist)
    else:
        snps = snp_to_z
    log("SNPs with z-score and dentist %s" % len(snps), DEBUG)
    return snps


def create_data():
    user_range, extended_range = set_range_restrictions()
    # Reading and setting data -----------------------------------------------------------------------------------------
    gene_to_log_bf = {}
    # Setting values from the GWAS study all SNPs
    snp_to_effect_allele, snp_to_other_allele, snp_to_z, snp_to_n, snp_to_beta, snp_to_se, snp_to_p = read_sumstats()
    # Allele correction using the alleles defined in the LD file
    allele_alignment(snp_to_effect_allele, snp_to_other_allele, snp_to_z, snp_to_n)

    # read Dentist file   ----------------------------------------------------------------------------------------------
    snp_to_analyze = add_dentist(snp_to_z)
    # ------------------------------------------------------------------------------------------------------------------

    # Sample size filter  ----------------------------------------------------------------------------------------------
    if len(snp_to_n) != 0:
        max_n = max(snp_to_n.values())
        to_delete = [snp for snp in snp_to_analyze if snp_to_n[snp] < options.n_threshold * max_n]
        before_snps = len(snp_to_analyze)
        for snp in to_delete:
            del snp_to_analyze[snp]
        log("Deleted %d SNPs due to sample size filtering at %.3g of max; %d remaining" % (
            before_snps - len(snp_to_analyze), options.n_threshold, len(snp_to_analyze)), DEBUG)
    # ------------------------------------------------------------------------------------------------------------------

    # fill all the dictionaries with the SNPs inside the extended region from the LD file
    chr_to_snp_pos, snp_to_chr_pos, chr_pos_to_snp, ld_data, ld_dict = read_ld_1(snp_to_analyze, extended_range)

    # create data for Locus Zoom   -------------------------------------------------------------------------------------
    # snp_to_chr_pos = read_bim_file()
    save_data_for_magma(snp_to_effect_allele, snp_to_n, snp_to_p,
                        snp_to_other_allele, snp_to_beta, snp_to_se, snp_to_chr_pos)
    save_data_for_locus_gwas(snp_to_effect_allele, snp_to_n, snp_to_p,
                             snp_to_other_allele, snp_to_beta, snp_to_se, snp_to_chr_pos)
    # ------------------------------------------------------------------------------------------------------------------

    # create data for external analysis  -------------------------------------------------------------------------------
    # snp_to_frq = read_frequencies()
    # save_data_for_COJO_ma(snp_to_effect_allele, snp_to_n, snp_to_p, snp_to_beta,
    #                       snp_to_other_allele, snp_to_beta, snp_to_se, snp_to_frq)
    # ------------------------------------------------------------------------------------------------------------------

    return gene_to_log_bf


def get_data():
    user_range, extended_range = set_range_restrictions()
    # Reading and setting data -----------------------------------------------------------------------------------------
    # Get chromosomes and position of the genes to be analyzed and defines genes able to get BF score
    gene_to_log_bf, gene_extend_to_chrom_start_stop = read_gene_file(user_range)
    # Setting values from the GWAS study all SNPs
    snp_to_effect_allele, snp_to_other_allele, snp_to_z, snp_to_n, snp_to_beta, snp_to_se, snp_to_p = read_sumstats()
    # Allele correction using the alleles defined in different file
    allele_alignment(snp_to_effect_allele, snp_to_other_allele, snp_to_z, snp_to_n)

    # read Dentist file   ----------------------------------------------------------------------------------------------
    snp_to_analyze = add_dentist(snp_to_z)
    # ------------------------------------------------------------------------------------------------------------------

    # Sample size filter  ----------------------------------------------------------------------------------------------
    if len(snp_to_n) != 0:
        max_n = max(snp_to_n.values())
        to_delete = [snp for snp in snp_to_analyze if snp_to_n[snp] < options.n_threshold * max_n]
        before_snps = len(snp_to_analyze)
        for snp in to_delete:
            del snp_to_analyze[snp]
        log("Deleted %d SNPs due to sample size filtering at %.3g of max; %d remaining" % (
            before_snps - len(snp_to_analyze), options.n_threshold, len(snp_to_analyze)), DEBUG)
    # ------------------------------------------------------------------------------------------------------------------

    # fill all the dictionaries with the SNPs inside the extended region from the LD file
    chr_to_snp_pos, snp_to_chr_pos, chr_pos_to_snp, ld_data, ld_dict = read_ld_1(snp_to_analyze, extended_range)

    # Using LD information to create components  -----------------------------------------------------------------------
    # component_to_snp = components_from_cluster(ld_dict, snp_to_chr_pos)
    snp_to_component, component_to_snp = read_ld_2(ld_data)
    del ld_data
    return gene_to_log_bf, component_to_snp, gene_extend_to_chrom_start_stop, snp_to_chr_pos, \
           snp_to_z, snp_to_n, ld_dict, snp_to_se, snp_to_beta


def simulate_GWAS():
    snp_to_z = {}
    snp_to_n = {}
    snp_to_beta = {}
    snp_to_se = {}
    snp_to_p = {}
    size = int(options.SNPs_simulated)
    cor_matrix_obs = np.diag(np.ones(size))
    mean = np.zeros(size)
    snp_zs = np.random.multivariate_normal(mean, cor_matrix_obs, 1).reshape([size])
    for i in range(size):
        rsid = "rs" + str(i)
        snp_to_z[rsid] = snp_zs[i]
        snp_to_n[rsid] = 10
        snp_to_beta[rsid] = .012
        snp_to_se[rsid] = 0.12
        snp_to_p[rsid] = 0.12
    return snp_to_z, snp_to_n, snp_to_beta, snp_to_se, snp_to_p


def simulate_ld(snp_to_z, extended_range):
    ld_data = []
    ld_dict = {}
    snp_to_chr_pos = {}
    chrom = None
    range_start = None
    range_end = None
    for key in extended_range:
        chrom = key
        range_start = extended_range[key][0]
        range_end = extended_range[key][1]
    mid = range_start + ((range_end - range_start) / 2)
    for snp in snp_to_z:
        snp_to_chr_pos[snp] = (chrom, mid)
        mid += 1
        ld_dict[snp, snp] = 1
        ld_data.append((snp, snp, 1))
    return snp_to_chr_pos, ld_data, ld_dict


def simulate_components(snp_to_z):
    component_size = options.max_component_size
    snp_to_component = {}
    component_to_snp = {}
    component = 1
    aux = 0
    for snp in snp_to_z:
        if aux == 0:
            component_to_snp[component] = set()
        aux += 1
        snp_to_component[snp] = component
        component_to_snp[component].add(snp)
        if aux == component_size:
            component += 1
            aux = 0
    return snp_to_component, component_to_snp


def simulate_data():
    user_range, extended_range = set_range_restrictions()
    # Get chromosomes and position of the genes to be analyzed and defines genes able to get BF score
    gene_to_log_bf, gene_extend_to_chrom_start_stop = gene_simulation(user_range)
    # Setting values from the GWAS study all SNPs
    snp_to_z, snp_to_n, snp_to_beta, snp_to_se, snp_to_p = simulate_GWAS()
    # fill all the dictionaries with the SNPs inside the extended region from the LD file
    snp_to_chr_pos, ld_data, ld_dict = simulate_ld(snp_to_z, extended_range)
    # create the components
    snp_to_component, component_to_snp = simulate_components(snp_to_z)
    del ld_data
    return gene_to_log_bf, component_to_snp, gene_extend_to_chrom_start_stop, snp_to_chr_pos, \
           snp_to_z, snp_to_n, ld_dict


def ln_prob_z(z_scores, causal_vector, corr_matrix):  # sample size of the reference panel
    num_samples_mc = 2000  # number of samples used for the integration
    size = corr_matrix.shape[0]
    # ==========================
    sigma_var_matrix = np.diag(causal_vector)
    ln_probs = []
    for i in range(num_samples_mc):
        try:
            val = invwishart.rvs(df=size, scale=corr_matrix, size=1)
            cov_matrix = val + np.dot(val, np.dot(sigma_var_matrix, val))
            cov_matrix = matrix_stabilization(cov_matrix)
            normal = multivariate_normal.logpdf(z_scores, np.zeros(size), cov_matrix)
            ln_probs.append(normal)
        except np.linalg.LinAlgError as e:
            log("Singular matrix", TRACE)
    max_ln_prob = max(ln_probs)
    sum_ln_prob = max_ln_prob + np.log(sum([np.exp(x - max_ln_prob) for x in ln_probs]))
    return sum_ln_prob + np.log(1 / num_samples_mc)


def ldpred_gene_window_prior(beta, n, correlation, p, h_norm, M, snp_to_index, snp_zs, snp_to_chr_pos, gene_to_snp_probs,
                       gene_to_snp_distance, gene_extend_to_chrom_start_stop, burn_in=50, gibbs_iter=100):
    """
    calculates gene probabilities using gibbs sampling, All the parameters are block wise

    :param beta normalized marginal effect sizes
    :param n sample size
    :param correlation LD matrix
    :param p fraction of causal SNPs
    :param h_norm heritability
    :param M total number of SNPs
    :param snp_to_index rsid of the SNPs in the component
    :param snp_zs used to compare posterior beta
    :param snp_to_chr_pos used to plot by genomic position
    :param gene_to_snp_probs used to plot in genomic order
    :param gene_to_snp_distance the distance of the snp to the gene
    :param gene_extend_to_chrom_start_stop position of the genes
    :param burn_in number of iterations dump from the mean
    :param gibbs_iter number of samples used for gibbs calculations
    """
    # ======== initialization ===================================
    exit_status = "success"
    total_iter = burn_in + gibbs_iter
    epsilon = 0.0001
    pi = options.gene_prior
    snps_in_component = beta.shape[0]
    num_genes = len(gene_to_snp_probs)
    norm_beta = beta  # / se)  # * math.sqrt(n)  used when the marginal betas aren't normalized
    shrinkage = 1 / (1 + ((M * p) / (n * h_norm)))
    p_gene = p
    epsilon_gene = epsilon
    mean = 50000
    s = 100
    # priors definition =================================
    # prior for window
    prior_mu = 50000
    prior_std = 25000
    window_prior_strength = 1000
    # variables for sigma update
    prior_sigma = h_norm / (M * p)
    prior_alpha = 10  # confidence on the prior sigma
    prior_beta = prior_sigma * prior_alpha
    posterior_sigma = prior_sigma
    # variables for p update
    prior_strength = 10
    prior_p_alpha = prior_strength
    prior_p_beta = (1 / p) * prior_strength
    # ======== variables for SNP to gene ================
    I = deepcopy(gene_to_snp_probs)
    # initialize I with gene to SNP probabilities
    for gene in I:
        for j in range(snps_in_component):
            if I[gene][j] == options.sigma:
                I[gene][j] = 1
            else:
                I[gene][j] = 0
    window_mu = np.zeros(num_genes)
    window_s = np.zeros(num_genes)
    window_mu_sum = np.zeros(num_genes)
    window_s_sum = np.zeros(num_genes)
    # ======== variables for plotting ===================
    window_mu_vec_print = np.zeros((total_iter, num_genes))
    window_s_vec_print = np.zeros((total_iter, num_genes))
    print_beta = np.zeros(total_iter)
    max_beta_idx = np.argmax(abs(beta))
    pi_vec_print = np.zeros((total_iter, num_genes))
    pi_avg_print = np.zeros((gibbs_iter, num_genes))
    p_vec_print = np.zeros(total_iter)
    p_in_vec_print = np.zeros(total_iter)
    p_out_vec_print = np.zeros(total_iter)
    heritability_vec_print = np.zeros(total_iter)
    posterior_sigma_vec_print = np.zeros(total_iter)
    ldpred_heritability_vec = np.zeros(total_iter)
    pi_pos_vec = np.zeros(total_iter)
    shrinkage_vec_print = np.zeros(total_iter)
    # ======== infinitesimal model =============================
    factor = np.diag(np.ones(snps_in_component) * (M / (n * h_norm)))
    inf_beta = np.matmul(np.linalg.inv(factor + correlation), norm_beta)
    # ======== gibbs sampler ===================================
    # variables to update inside the gibbs loop
    omega_vec = np.zeros(snps_in_component)
    omega_cap_vec = np.zeros(snps_in_component)
    joint_beta = np.zeros(snps_in_component)
    p_vec = np.zeros(snps_in_component)
    sampled_beta = np.copy(norm_beta)
    c_vec = np.zeros(snps_in_component)
    pi_vec = np.zeros(num_genes)
    pi_vec_sum = np.zeros(num_genes)
    # G is initialized so all genes are causal
    sampled_G = {}
    for gene in gene_to_snp_probs:
        sampled_G[gene] = 1  # np.random.binomial(1, 0.05, 1)[0]  # is better to initialize all genes to 1
    # ======== aux variables ===========================
    gene_list = list(gene_to_snp_probs)
    plt.ion()
    fig, ax = plt.subplots()
    x, y = [], []
    plot_scatter = ax.scatter(x, y)
    plot_line, = ax.plot(x, y)
    plt.xlim(0, 10)
    plt.ylim(0, 1)
    plt.draw()
    plt.show()
    # ======= gibbs iterations =========================
    for k in range(burn_in + gibbs_iter):
        log("iter: {}".format(k))
        # ==UPDATE== SNP causal status, true betas
        for j in range(snps_in_component):
            # calculation of beta tilde j (joint beta) ========================
            corr_no_j = np.delete(correlation[j], j, 0)
            sampled_beta_no_j = np.delete(sampled_beta, j, 0)
            joint_beta[j] = norm_beta[j] - np.dot(sampled_beta_no_j.T, corr_no_j)
            # Compute pj and sample cj  ==================================
            if S2CG(j, I, sampled_G):
                # p_vec[j] = 1 / (1 + (((1 - p) / p) * math.sqrt(1 + ((n * h_norm) / (M * p))) * math.exp(
                #     ((-1 / 2) * n * joint_beta[j] ** 2) /
                #     (1 + ((M * p) / (n * h_norm))))))
                p_vec[j] = 1 / (1 + (((1 - p) / p) * math.sqrt(1 + (n * posterior_sigma))
                                     * math.exp(
                            (-1 / 2) * ((n * joint_beta[j] ** 2) / (1 + (1 / (n * posterior_sigma)))))))
            else:
                p_vec[j] = epsilon
            c_vec[j] = np.random.binomial(1, p_vec[j], 1)[0]
            # Sample bj  ======================================
            mean = joint_beta[j] * shrinkage
            if c_vec[j]:
                variance = (1 / n) * shrinkage
                sampled_beta[j] = np.random.normal(mean, math.sqrt(variance), 1)
            else:
                sampled_beta[j] = 0.0
            omega_vec[j] = p_vec[j] * mean
        # ==UPDATE== gene causal status
        for l in range(num_genes):
            gene = gene_list[l]
            iter_sampled_g = sampled_G.copy()
            iter_sampled_g[gene] = 1  # setting the current gene to one
            multiplication_1 = 1
            for j in range(snps_in_component):
                if S2CG(j, I, iter_sampled_g):
                    if c_vec[j] == 1:
                        multiplication_1 *= p_gene
                    else:
                        multiplication_1 *= (1 - p_gene)
                else:
                    if c_vec[j] == 1:
                        multiplication_1 *= epsilon_gene
                    else:
                        multiplication_1 *= (1 - epsilon_gene)
            iter_sampled_g[gene] = 0  # setting the current gene to zero
            multiplication_0 = 1
            for j in range(snps_in_component):
                if S2CG(j, I, iter_sampled_g):
                    if c_vec[j] == 1:
                        multiplication_0 *= p_gene
                    else:
                        multiplication_0 *= (1 - p_gene)
                else:
                    if c_vec[j] == 1:
                        multiplication_0 *= epsilon_gene
                    else:
                        multiplication_0 *= (1 - epsilon_gene)
            pi_denominator = ((multiplication_1 * pi) + (multiplication_0 * (1 - pi)))
            if pi_denominator != 0:
                pi_vec[l] = (multiplication_1 * pi) / pi_denominator
            else:
                pi_vec[l] = 0
            sampled_G[gene] = np.random.binomial(1, pi_vec[l], 1)[0]
        # ==UPDATE== p, heritability, p_in, p_out
        update = True
        if update:
            # update gene probability
            L_c = np.count_nonzero(sampled_G.values())
            pi = np.random.beta(1 + L_c, 19 + len(gene_list) - L_c)  # 1 and 19 set a week prior of 0.05
            # update of fraction of causal SNPs in the region
            M_c = np.count_nonzero(sampled_beta)
            p = np.random.beta(prior_p_alpha + M_c, prior_p_beta + snps_in_component - M_c)
            # update the number of SNPs in the region
            M = snps_in_component
            # calculation of fraction of causal SNPs p_in inside of causal genes
            M_c_in = 0
            M_c_out = 0
            snps_inside_causal_genes = 0
            for j in range(snps_in_component):
                if S2CG(j, I, sampled_G):  # inside the genes
                    snps_inside_causal_genes += 1
                    if sampled_beta[j] != 0:
                        M_c_in += 1
                else:  # outside genes
                    if sampled_beta[j] != 0:
                        M_c_out += 1
            snps_outside_causal_genes = snps_in_component - snps_inside_causal_genes
            p_in = np.random.beta(1 + M_c_in, 1 + snps_inside_causal_genes - M_c_in)
            p_out = np.random.beta(1 + M_c_out, 1 + snps_outside_causal_genes - M_c_out)
            # print("{},{}".format(p_in, p_out))
            p_gene = p_in
            epsilon_gene = p_out
            # epsilon_gene = 0.0001  # p_out
            # print(epsilon_gene)
            # update of heritability
            ldpred_h_norm = np.matmul(np.matmul(sampled_beta.T, correlation), sampled_beta)
            pos_alpha = (prior_alpha + (M / 2))
            pos_beta = (prior_beta + (np.matmul(sampled_beta.T, sampled_beta) / 2))
            posterior_sigma = sc.stats.invgamma.rvs(a=pos_alpha, scale=pos_beta, size=1)[0]
            h_norm = np.matmul(sampled_beta.T, sampled_beta)  # posterior_sigma * M  # heritability is no longer used in gibbs it's used to compare
            # =======================================================================
        # ==UPDATE== window and S2G link
        window_update = True
        if window_update:
            # updating window mean and s
            regression_X = np.empty(0)
            regression_y = np.empty(0)
            for gene in gene_list:
                # gene_X = np.concatenate((gene_to_snp_distance[gene], np.negative(gene_to_snp_distance[gene]))).reshape(-1, 1)
                # gene_y = np.concatenate((I[gene], np.ones(snps_in_component)))
                # gene_X = gene_to_snp_distance[gene].reshape(-1, 1)
                # gene_y = I[gene]
                gene_X = np.concatenate((gene_to_snp_distance[gene], np.zeros(snps_in_component)))
                gene_y = np.concatenate((I[gene], np.ones(snps_in_component)))
                regression_X = np.concatenate((regression_X, gene_X))
                regression_y = np.concatenate((regression_y, gene_y))
            regression_X = regression_X.reshape(-1, 1)
            with warnings.catch_warnings():
                warnings.filterwarnings('error')
                try:
                    clf = LogisticRegression(random_state=0).fit(regression_X, regression_y)
                    mean = - clf.intercept_ / clf.coef_
                    s = 1 / clf.coef_
                    if abs(s) < 1300:
                        s = -1300
                        log("s reached the minimum possible s ={}, s was set to -1500".format(s), DEBUG)
                    window_mu = float(mean) * np.ones(num_genes)  # put the same window in all genes
                    window_s = float(s) * np.ones(num_genes)
                except Warning:
                    log("Window not updated for gene {}, failed to converge, using previous value".format(gene))
            if True:  # used to check the evolution of the window for a gene
                num_points = 1000
                x_axis = np.linspace(min(regression_X), max(regression_X), num_points)
                print(s)
                logit = 1 / (1 + np.exp(-(x_axis - mean) / s)).reshape(num_points, )
                x = regression_X.reshape(regression_X.shape[0], ).tolist()
                y = regression_y.tolist()
                plt.xlim(min(regression_X), max(regression_X))
                plot_scatter.set_offsets(np.c_[x, y])
                plot_line.set_data(x_axis, logit)
                fig.canvas.draw_idle()
                plt.pause(0.001)
                # plt.waitforbuttonpress()

            # update I
            for gene in I:
                l = gene_list.index(gene)
                for j in range(snps_in_component):
                    dis_j = gene_to_snp_distance[gene][j]
                    gene_mu = window_mu[l]
                    gene_s = window_s[l]
                    exponent = np.float64(-(dis_j - gene_mu)/gene_s)
                    # exponent = np.float64(-(dis_j - window_mu[l]) / window_s[l])
                    if abs(exponent) > 100:
                        prior_p_I_l_j_1 = 0
                    else:
                        prior_p_I_l_j_1 = 1 / (1 + np.exp(exponent))

                    I_aux = I.copy()
                    I_aux[gene][j] = 1
                    c_likelihood_1 = 1
                    if S2CG(j, I_aux, sampled_G):
                        if c_vec[j] == 1:
                            c_likelihood_1 *= p_gene
                        else:
                            c_likelihood_1 *= (1 - p_gene)
                    else:
                        if c_vec[j] == 1:
                            c_likelihood_1 *= epsilon_gene
                        else:
                            c_likelihood_1 *= (1 - epsilon_gene)
                    I_aux[gene][j] = 0
                    c_likelihood_0 = 1
                    if S2CG(j, I_aux, sampled_G):
                        if c_vec[j] == 1:
                            c_likelihood_0 *= p_gene
                        else:
                            c_likelihood_0 *= (1 - p_gene)
                    else:
                        if c_vec[j] == 1:
                            c_likelihood_0 *= epsilon_gene
                        else:
                            c_likelihood_0 *= (1 - epsilon_gene)
                    link_numerator = c_likelihood_1 * prior_p_I_l_j_1
                    link_denominator = (c_likelihood_0 * (1 - prior_p_I_l_j_1)) + link_numerator
                    if link_denominator == 0 or math.isnan(link_denominator):
                        posterior_p_I_l_j_1 = 0
                    else:
                        # print("num: {}, denom: {}".format(link_numerator, link_denominator))
                        posterior_p_I_l_j_1 = link_numerator / link_denominator
                    # sampling new I_j_l
                    I[gene][j] = np.random.binomial(1, posterior_p_I_l_j_1, 1)[0]
        # ======= Sums after burn-in for average calculations =============
        if k >= burn_in:  # used to calculate the results by eliminated the first iterations
            omega_cap_vec += omega_vec
            pi_vec_sum += pi_vec
            window_mu_sum += window_mu
            window_s_sum += window_s
            # setting the averages  =================
            pi_avg_print[k - burn_in] = pi_vec_sum / (k - burn_in + 1)
            # shrinkage = shrinkage_sum / (k - burn_in + 1)
        # ======= Update of the traces ===================================
        pi_vec_print[k] = pi_vec
        print_beta[k] = omega_vec[max_beta_idx]
        p_vec_print[k] = p
        heritability_vec_print[k] = h_norm
        posterior_sigma_vec_print[k] = posterior_sigma
        ldpred_heritability_vec[k] = ldpred_h_norm
        p_in_vec_print[k] = p_gene
        p_out_vec_print[k] = epsilon_gene
        pi_pos_vec[k] = pi
        shrinkage_vec_print[k] = (n * posterior_sigma) / ((n * posterior_sigma) + 1)
        window_mu_vec_print[k] = window_mu
        window_s_vec_print[k] = window_s
        # === Collapse detector ================
        limit = 1e-07  # 1e-07
        if posterior_sigma < limit:
            log("chain collapsed", INFO)
            exit_status = "collapsed"
            break
    # === end of gibbs =================================
    # === calculation of posterior betas and gene prob  =================
    posterior_beta = omega_cap_vec / gibbs_iter
    posterior_gene_prob = pi_vec_sum / gibbs_iter
    # ==== plotting results ================================
    if options.debug_level >= 2 and exit_status == "success":
        plt.figure(1)
        plt.subplot(2, 2, 1)
        plt.plot(print_beta, linewidth=1)
        plt.title("=TRACE= convergence of strongest SNP", fontsize=8)
        plt.subplot(2, 2, 2)
        plt.scatter(inf_beta, posterior_beta, s=2)
        plt.title("infinitesimal vs Gibbs results", fontsize=8)
        # sorting snps to plot by genomic position
        sort_snps = dict(sorted(snp_to_chr_pos.items(), key=lambda item: item[1]))
        to_plot = np.zeros(snps_in_component)
        j = 0
        for snp in sort_snps:
            if snp in snp_to_index:
                to_plot[j] = posterior_beta[snp_to_index[snp]]
                j += 1
        plt.subplot(2, 2, 3)
        plt.scatter(np.linspace(0, 1, snps_in_component), to_plot, s=2)
        plt.title("=AVG minus burn-in= gibbs results", fontsize=8)
        # plot of gene probabilities
        ax = plt.subplot(2, 2, 4)
        plt.title("=TRACE= gene posteriors", fontsize=8)
        gene_file = open('../data/ldpred_results/gene.results', 'w')
        NUM_COLORS = num_genes
        cm = plt.get_cmap('gist_rainbow')
        line_type = ["-", "--", "-.", ":"]
        linecycler = cycle(line_type)
        ax.set_prop_cycle('color', [cm(1. * i / NUM_COLORS) for i in range(NUM_COLORS)])
        for l in range(num_genes):
            ax.plot(pi_vec_print[:, l], next(linecycler), label=gene_list[l])
            gene_file.write("{}\t{}\t{}\t{}\n".format(gene_list[l], posterior_gene_prob[l], window_mu_sum[l] / gibbs_iter, window_s_sum[l] / gibbs_iter))
        # plt.legend(loc='upper center', bbox_to_anchor=(0.5, 1.1), ncol=3, fancybox=True, shadow=True)
        gene_file.close()
        plt.tight_layout(pad=0.5)
        # plotting other probabilities
        plt.figure(2)
        plt.subplot(3, 1, 1)
        plt.plot(p_vec_print)
        plt.title("=TRACE= Prior: {:.4f}. fraction of causal SNPs p in the region, mean {:.4f}".format(options.p, p_vec_print.mean()), fontsize=8)
        plt.subplot(3, 1, 2)
        plt.plot(p_in_vec_print)
        plt.title("=TRACE= fraction of causal SNPs p_in given the gene is causal, mean {:.4f}".format(p_in_vec_print.mean()),
                  fontsize=8)
        plt.subplot(3, 1, 3)
        plt.plot(p_out_vec_print)
        plt.title(
            "=TRACE= fraction of causal SNPs p_out given the gene is not causal, mean {:.4f}".format(p_out_vec_print.mean()),
            fontsize=8)
        plt.tight_layout(pad=0.5)
        plt.figure(3)
        plt.subplot(3, 1, 1)
        plt.plot(heritability_vec_print)
        plt.title("=TRACE= My heritability, mean {:.8f}".format(heritability_vec_print.mean()), fontsize=8)
        plt.subplot(3, 1, 2)
        plt.plot(ldpred_heritability_vec)
        plt.title("=TRACE= LDpred heritability, mean {:.6f}".format(ldpred_heritability_vec.mean()), fontsize=8)
        plt.subplot(3, 1, 3)
        plt.plot(pi_pos_vec)
        plt.title("=TRACE= Posterior gene probability, mean {:.6f}".format(pi_pos_vec.mean()), fontsize=8)
        plt.tight_layout(pad=0.5)
        plt.figure(4)
        plt.subplot(3, 1, 1)
        plt.plot(shrinkage_vec_print)
        plt.title("=TRACE= shrinkage , final: {:.4f}, used: {:.4f}".format(shrinkage_vec_print[total_iter-1], shrinkage), fontsize=8)
        plt.subplot(3, 1, 2)
        plt.plot(posterior_sigma_vec_print)
        plt.title("=TRACE= Prior sigma: {:.8f}. Posterior sigma mean {:.8f}".format(prior_sigma, posterior_sigma_vec_print.mean()), fontsize=8)
        plt.tight_layout(pad=0.5)
        plt.figure(5)
        plt.plot(pi_avg_print)
        plt.title("=AVG= Gene posterior prob")
        plot_genes_pos(gene_extend_to_chrom_start_stop, window_mu_sum / gibbs_iter, gene_list)
        plot_window_trace = False
        if plot_window_trace:
            figure = 5
            for gene in gene_list:
                l = gene_list.index(gene)
                plt.figure(figure)
                figure += 1
                ax1 = plt.subplot(2, 1, 1)
                plt.title("window mu, gene {}, mean {}".format(gene, np.mean(window_mu_vec_print[:, l])), fontsize=8)
                ax2 = plt.subplot(2, 1, 2)
                plt.title("window s, gene {}, mean {}".format(gene, np.mean(window_s_vec_print[:, l])), fontsize=8)
                line_type = next(linecycler)
                ax1.plot(window_mu_vec_print[:, l], line_type, label=gene_list[l])
                ax2.plot(window_s_vec_print[:, l], line_type, label=gene_list[l])
                plt.legend(loc='lower right', bbox_to_anchor=(0.5, 1.1), ncol=3, fancybox=True, shadow=True)
                plt.tight_layout(pad=0.5)
        plt.show()
    # saving data ======================================================
    for snp in snp_to_index:
        idx = snp_to_index[snp]
        chr = snp_to_chr_pos[snp][0]
        pos = snp_to_chr_pos[snp][1]
        line = "{}\t{}\t{}\tA\tG\t{}\t{}\t{}\t{}\n".format(snp, chr, pos, posterior_beta[idx], p_vec[idx],
                                                           snp_zs[idx], abs(posterior_beta[idx]))
        file_ldpred.write(line)
    print(posterior_sigma_vec_print)
    plt.waitforbuttonpress()
    return posterior_beta, exit_status


def plot_genes_pos(gene_pos, window, gene_list):
    open_data = []
    close_data = []
    high_data = []
    low_data = []
    gene_names = []
    sort_genes = dict(sorted(gene_pos.items(), key=lambda item: item[1]))
    for gene in sort_genes:
        if gene in gene_list:
            l = gene_list.index(gene)
            gene_names.append(gene)
            f_gene_chrom, f_gene_start, f_gene_end = gene_pos[gene]
            open_data.append(f_gene_start)
            close_data.append(f_gene_end)
            high_data.append(f_gene_end + window[l])
            low_data.append(f_gene_start - window[l])
    fig = go.Figure(data=[go.Candlestick(x=gene_names,
                                         open=open_data, high=high_data,
                                         low=low_data, close=close_data)])

    fig.show()


def S2CG(j, I, G):
    """
    returns True if SNP j is in a causal gene, False otherwise
    """
    for gene in I:
        if G[gene] == 1:  # if is a causal gene
            if I[gene][j] == 1:  # and snpj is link to gene
                # for the moment we only care if the snp is in one of the genes in the set we don't account
                # for a SNP overlapping more than one gene
                return True
    return False


def m_snps_sim(correlation):
    M = correlation.shape[0]
    n = 100000
    maf_v = np.random.uniform(0, 1, M)
    samples = list()
    for j in range(M):
        aux = np.random.binomial(2, maf_v[j], n)
        mean = np.mean(aux)
        std = np.std(aux)
        # samples.append(aux)
        samples.append((aux - mean) / std)
    genotype = np.asarray(samples)

    L = np.linalg.cholesky(correlation)
    corr_x = np.matmul(L, genotype)
    for j in range(M):
        corr_x[j] = (corr_x[j] - corr_x[j].mean()) / corr_x[j].std()

    # assigment of heritability and sigma_e assuming var[y] = 1
    h_2 = 0.01212
    sigma_e = 1 - h_2
    # true_beta = np.random.normal(0, np.sqrt(h_2/M), M)
    true_beta = np.ones(M) * 0.01
    true_beta[231] = 0.04  # manually adding a causal SNP
    Xb = np.matmul(corr_x.T, true_beta)
    noise = np.random.normal(0, np.sqrt(sigma_e), n)
    y = Xb + noise
    y = (y - y.mean())

    h_2_calc = Xb.var() / (Xb.var() + noise.var())
    heritability = sum(true_beta ** 2)

    noise_var_calc = noise.var() / (Xb.var() + noise.var())

    marginal_beta = np.zeros(len(maf_v))
    marginal_var = np.zeros(len(maf_v))
    for j in range(M):
        marginal_beta[j] = np.matmul(corr_x[j].reshape((1, n)), y.reshape((n, 1))) / n
        marginal_var[j] = sum((y - (corr_x[j] * marginal_beta[j])) ** 2) / (n - 2)
    marginal_se = np.sqrt(marginal_var / n)
    marginal_std = np.sqrt(marginal_var)

    marginal_z = marginal_beta / marginal_se

    plt.scatter(true_beta, marginal_beta)
    plt.title("true betas vs marginal betas")
    plt.show()

    marginal_p = np.zeros(M)
    for j in range(len(maf_v)):
        marginal_p[j] = sc.stats.norm.sf(abs(marginal_z[j])) * 2

    rs_ids = ['rs{}'.format(val) for val in range(M)]
    file_name = '../data/simulations/ldpred_sim_M_{}'.format(M)
    file = open(file_name, 'w')
    file.write("MarkerName\tbeta\tse\tz\tp\tWeight\n")
    for j in range(len(maf_v)):
        file.write(
            "{}\t{}\t{}\t{}\t{}\t{}\n".format(rs_ids[j], marginal_beta[j], marginal_se[j], marginal_z[j], marginal_p[j],
                                              n))

    print("h_2_calc: {} h sum of square beta: {}".format(h_2_calc, heritability))
    file.close()
    return rs_ids, marginal_beta, marginal_se, marginal_z, marginal_p, n, M, true_beta, h_2_calc, correlation
