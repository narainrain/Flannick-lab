from optparse import OptionParser
import sys
import gzip
import matplotlib.pyplot as plt
import numpy as np
from math import *
from scipy.optimize import minimize
import pickle as pl
import math
import random
import scipy as sc


# constant variables
NONE = 0
INFO = 1
DEBUG = 2
TRACE = 3


def arg_settings():
    parser = OptionParser("usage: %prog [options]")
    # General options  ====================================================
    parser.add_option("", "--betas-out-file", default=None)
    parser.add_option("", "--genes-out-file", default=None)
    parser.add_option("", "--data-out-file", default=None)
    parser.add_option("", "--tmp-folder", default=None)
    parser.add_option("", "--out-folder", default=None)
    parser.add_option("", "--debug-level", type='int', default=0)
    # Gibbs options  ====================================================
    parser.add_option("", "--num-chains", type=int, default=20)
    parser.add_option("", "--burn-in", type=int, default=50)
    parser.add_option("", "--max-iter", type=int, default=1000)
    parser.add_option("", "--convergence-thr", type=float, default=1.04)
    # Region definition  ====================================================
    parser.add_option("", "--chrom", default=None)  # filter analysis to a chromosome
    parser.add_option("", "--start", type=int, default=None)  # filter analysis to a beginning position
    parser.add_option("", "--end", type=int, default=None)  # filter analysis to an ending position
    parser.add_option("", "--range", default=None)  # filter analysis to a chr:start-stop range
    parser.add_option("", "--cluster", action="store_true", default=False)
    parser.add_option("", "--job-id", type="int")
    parser.add_option("", "--block-size", type="int", default=1000000)
    parser.add_option("", "--extension", type=int, default=500000)
    # Gene file Must have GENE, CHROM, START, STOP ===========================
    parser.add_option("", "--gene-file", default=None)
    parser.add_option("", "--gene-id-col", default=None)
    parser.add_option("", "--gene-chrom-col", default=None)
    parser.add_option("", "--gene-start-col", default=None)
    parser.add_option("", "--gene-end-col", default=None)
    # Sumstats file Must have ID, BETA, and SE or PVALUE ======================
    parser.add_option("", "--sumstats-file", default=None)
    parser.add_option("", "--tabix", default=None)
    parser.add_option("", "--sumstats-id-col", default=None)
    parser.add_option("", "--sumstats-beta-col", default=None)
    parser.add_option("", "--sumstats-se-col", default=None)
    parser.add_option("", "--sumstats-p-col", default=None)
    parser.add_option("", "--sumstats-z-col", default=None)
    parser.add_option("", "--sumstats-n-col", default=None)
    parser.add_option("", "--sumstats-effect-allele-col", default="Allele1")
    parser.add_option("", "--sumstats-other-allele-col", default="Allele2")
    # ld file ==================================================================
    parser.add_option("", "--ld-file", default=None)
    parser.add_option("", "--ld-folder", default=None)
    parser.add_option("", "--ld-id1-col", default="SNP_A")
    parser.add_option("", "--ld-chrom1-col", default="CHR_A")
    parser.add_option("", "--ld-pos1-col", default="BP_A")
    parser.add_option("", "--ld-id2-col", default="SNP_B")
    parser.add_option("", "--ld-chrom2-col", default="CHR_B")
    parser.add_option("", "--ld-pos2-col", default="BP_B")
    parser.add_option("", "--ld-r-col", default="R2")
    parser.add_option("", "--all-positive", action="store_true")
    parser.add_option("", "--max-component-size", type="int", default=1000)
    parser.add_option("", "--stab-term-frac", type=float, default=1e-4)  # used to correct ld matrix
    parser.add_option("", "--min-ld-threshold", default=0, type="float")
    parser.add_option("", "--max-ld-threshold", default=1, type="float")
    # Extra file with trusted alleles to fix sumstats alleles ======================
    parser.add_option("", "--ld-allele-file", default=None)  # Must have ID, ALLELE
    parser.add_option("", "--ld-allele-id-col", default=2)
    parser.add_option("", "--ld-allele-effect-allele-col", default=5)
    parser.add_option("", "--ld-allele-other-allele-col", default=6)
    # ld allele correction =========================================================
    parser.add_option("", "--detect-flips", type='float', default=10)
    # Dentist filtering ============================================================
    parser.add_option("", "--dentist-file", default=None)
    parser.add_option("", "--dentist-folder", default=None)
    parser.add_option("", "--dentist-thr", type="float", default=float('+inf'))
    # Sample size filtering  =======================================================
    parser.add_option("", "--n-threshold", type=float, default=0)
    # Model options ================================================================
    parser.add_option("", "--p", type="float", default=0.001)
    parser.add_option("", "--heritability", type="float", default=0.1212)
    parser.add_option("", "--sample-size", type="float", default=84685)
    parser.add_option("", "--total_num_SNPs", type="float", default=84700000)  # number of SNPs in 1000G
    parser.add_option("", "--gene-prior", type=float, default=0.05)
    parser.add_option("", "--causal-window", type=int, default=10)  # in kilo bases
    # options for window-model: one_window, multiple_windows, bayesian_one_window, bayesian_multiple_windows
    parser.add_option("", "--window-model", default="one_window")
    (options, _) = parser.parse_args()
    options.causal_window = options.causal_window * 1000  # changed to kilo bases
    return options


def bail(message):
    sys.stderr.write("%s\n" % message)
    sys.exit(1)


def clean_chrom(f_chrom):
    if len(f_chrom) > 3 and f_chrom[:3] == "chr":
        f_chrom = f_chrom[3:]
    return f_chrom


def open_gz(file):
    if file[-3:] == ".gz":
        return gzip.open(file, 'rt')
    else:
        return open(file)


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
                bail("Could not find match for column %s in header: %s" % (
                    col_name_or_index, "\t".join(header_cols)))
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


def S2CG(j, I, G):
    """
    returns True if SNP j is in a causal gene, False otherwise
    """
    for l in range(I.shape[0]):
        if G[l] == 1:  # if is a causal gene
            if I[l][j] == 1:  # and snpj is link to gene
                # for the moment we only care if the snp is in one of the genes in the set we don't account
                # for a SNP overlapping more than one gene
                return True
    return False


def create_component_index(snps_in_component):
    snp_to_index = {}
    index = 0
    for snp in sorted(snps_in_component):  # snp_to_chr_pos contains all SNPs inside the region with LD
        snp_to_index[snp] = index
        index += 1
    return snp_to_index


# @jit(nopython=True) Optimized for: MULTI-CHAIN, NUMBA
def sigmoid_f(x_value):
    """
    sigmoid function
    """
    return 1 / (1 + np.exp(-x_value))


def derivative(w, data, points, prior_beta0, prior_beta1, strength_beta0, strength_beta1):
    """
    Derivative of the negative posterior probability of the logistic parameters
    """
    m0 = np.asarray([prior_beta0, prior_beta1])
    s0_inv = np.asarray([[strength_beta0, 0], [0, strength_beta1]])
    aux = np.matmul(s0_inv, w - m0)
    sigmoid_res = sigmoid_f(w[0] + w[1] * data)
    total_sum = 0
    for j in range(len(points)):
        total_sum += (sigmoid_res[j] - points[j]) * data[j]
    return total_sum - aux


def negative_window_pos_apx(w, data, points, prior_beta0, prior_beta1, strength_beta0, strength_beta1):
    """
    Definition of the approximate negative posterior probability of the logistic parameters
    """
    m0 = np.asarray([prior_beta0, prior_beta1])
    s0_inv = np.asarray([[strength_beta0, 0], [0, strength_beta1]])
    aux = -(1 / 2) * np.matmul(np.matmul((w - m0).T, s0_inv), w - m0)
    total_sum = 0
    sigmoid_res = sigmoid_f(w[0] + w[1] * data)
    for n in range(len(points)):
        log_0 = log(sigmoid_res[n])
        log_value = 1 - sigmoid_res[n]
        if log_value == 0:
            log_value = 3.05902227e-07
        log_1 = log(log_value)
        total_sum += (points[n] * log_0) + ((1 - points[n]) * log_1)
    return -(aux + total_sum)


def bayesian_logistic_regression(data, points, prior_center, prior_decay, strength_beta0, strength_beta1,
                                 regression_method):
    """
    Returns the posterior center and decay of a Bayesian logistic regression
    """
    p_beta0 = - prior_center / prior_decay
    p_beta1 = 1 / prior_decay
    w0 = np.asarray([p_beta0, p_beta1])
    if regression_method == "BFGS":
        res = minimize(negative_window_pos_apx, w0, method='BFGS', jac=derivative,
                       args=(data, points, p_beta0, p_beta1, strength_beta0, strength_beta1),
                       options={'disp': True})
    else:
        res = minimize(negative_window_pos_apx, w0, method='nelder-mead',
                       args=(data, points, p_beta0, p_beta1, strength_beta0, strength_beta1),
                       options={'disp': True})
    pos_center = - res.x[0] / res.x[1]
    pos_decay = 1 / res.x[1]
    return pos_center, pos_decay


def plot_data(file=None, data=None):
    if file is not None:
        with open(file, 'rb') as f:
            data = pl.load(f)
    stop = data.gibbs_iter + data.burn_in
    # =======================================
    n = data.options.sample_size
    sigma_2_e = 1 - data.options.heritability
    norm_beta = (data.marginal_betas / data.marginal_se) / np.sqrt(n)
    # creation of z -scores ====================
    marginal_z_scores = data.marginal_betas / data.marginal_se
    norm_z_scores = norm_beta * np.sqrt(n)
    gibbs_z_scores = data.posterior_beta * np.sqrt(n)
    # creation of p-values ====================
    marginal_p_values = np.zeros(data.num_snps)
    norm_p_values = np.zeros(data.num_snps)
    gibbs_p_values = np.zeros(data.num_snps)
    for j in range(data.num_snps):
        marginal_p_values[j] = sc.stats.norm.sf(abs(marginal_z_scores[j])) * 2
        norm_p_values[j] = sc.stats.norm.sf(abs(norm_z_scores[j])) * 2
        gibbs_p_values[j] = sc.stats.norm.sf(abs(gibbs_z_scores[j])) * 2
    # ==================================================================================================================
    plt.figure(1)
    plt.subplot(2, 2, 1)
    plt.plot(data.print_beta[0:stop], linewidth=1)
    plt.title("=TRACE= convergence of strongest SNP chain 0", fontsize=8)
    plt.subplot(2, 2, 2)
    plt.scatter(data.inf_beta, data.posterior_beta, s=2)
    plt.axhline(y=0, color="grey", linestyle="--", linewidth=1)
    plt.axvline(color="grey", linestyle="--", linewidth=1)
    plt.title("infinitesimal vs Gibbs results all chains", fontsize=8)
    # sorting snps to plot by genomic position
    sort_snps = dict(sorted(data.snp_to_chr_pos.items(), key=lambda item: item[1]))
    to_plot = np.zeros(data.num_snps)
    snp_loc_to_plot = np.zeros(data.num_snps)
    marginal_betas_to_plot = np.zeros(data.num_snps)
    norm_betas_to_plot = np.zeros(data.num_snps)
    marginal_p_values_to_plot = np.zeros(data.num_snps)
    norm_p_values_to_plot = np.zeros(data.num_snps)
    gibbs_p_values_to_plot = np.zeros(data.num_snps)
    j = 0
    for snp in sort_snps:
        if snp in data.snps_in_region_index:
            idx = data.snps_in_region_index[snp]
            snp_loc_to_plot[j] = data.snp_to_chr_pos[snp][1]
            to_plot[j] = data.posterior_beta[idx]
            marginal_betas_to_plot[j] = data.marginal_betas[idx]
            norm_betas_to_plot[j] = norm_beta[idx]
            marginal_p_values_to_plot[j] = marginal_p_values[idx]
            norm_p_values_to_plot[j] = norm_p_values[idx]
            gibbs_p_values_to_plot[j] = gibbs_p_values[idx]
            j += 1
    plt.subplot(2, 2, 3)
    plt.scatter(snp_loc_to_plot, to_plot, s=2)
    plt.title("=AVG minus burn-in= gibbs results all chains", fontsize=8)
    # plot of gene probabilities
    plt.subplot(2, 2, 4)
    plt.plot(data.max_beta_pos_prob[0:stop], linewidth=1)
    plt.title("=TRACE= strongest SNP posterior prob chain 0", fontsize=8)
    plt.tight_layout(pad=0.5)
    # plotting other probabilities
    plt.figure(2)
    plt.subplot(3, 1, 1)
    plt.plot(data.p_vec_print.T[0:stop])
    plt.title("=TRACE= Prior: {:.2e}. "
              "fraction of causal SNPs p in the region, "
              "mean all chains {:.2e}".format(data.options.p, data.p_vec_print.T[0:stop].mean()),
              fontsize=8)
    plt.subplot(3, 1, 2)
    plt.plot(data.p_in_vec_print.T[0:stop])
    plt.title("=TRACE= Prior: {:.2e}. "
              "fraction of causal SNPs p_in given the gene is causal, "
              "mean {:.2e}".format(data.options.p, data.p_in_vec_print.T[0:stop].mean()),
              fontsize=8)
    plt.subplot(3, 1, 3)
    plt.plot(data.p_out_vec_print.T[0:stop])
    plt.title("=TRACE= Prior: {:.2e}. "
              "fraction of causal SNPs p_out given the gene is not causal, "
              "mean {:.2e}".format(data.epsilon, data.p_out_vec_print.T[0:stop].mean()),
              fontsize=8)
    plt.tight_layout(pad=0.5)
    plt.figure(3)
    plt.subplot(3, 1, 1)
    plt.plot(data.heritability_vec_print.T[0:stop])
    plt.title("=TRACE= My heritability, mean {:.2e}".format(data.heritability_vec_print.T[0:stop].mean()), fontsize=8)
    plt.subplot(3, 1, 2)
    plt.plot(data.ldpred_heritability_vec.T[0:stop])
    plt.title("=TRACE= LDpred heritability, mean {:.2e}".format(data.ldpred_heritability_vec.T[0:stop].mean()), fontsize=8)
    plt.subplot(3, 1, 3)
    plt.plot(data.pi_pos_vec.T[0:stop])
    plt.title("=TRACE= Prior: {:.2e}, Posterior gene probability, mean {:.2e}".format(data.options.gene_prior,
                                                                                  data.pi_pos_vec.T[0:stop].mean()), fontsize=8)
    plt.tight_layout(pad=0.5)
    plt.figure(4)
    plt.subplot(4, 1, 1)
    plt.plot(data.shrinkage_vec_print.T[0:stop])
    prior_shrinkage = (1 / (1 + ((data.options.total_num_SNPs * data.options.p) /
                                 (data.options.sample_size * data.options.heritability))))
    inf_shrinkage = (1 / (1 + (data.options.total_num_SNPs /
                               (data.options.sample_size * data.options.heritability))))
    plt.title("=TRACE= Shrinkage Inf: {:.2e}, Prior: {:.2e}, mean: {:.2e}".format(inf_shrinkage,
                                                                                  prior_shrinkage,
                                                                                  data.shrinkage_vec_print.T[0:stop].mean()),
              fontsize=8)
    plt.subplot(4, 1, 2)
    plt.plot(data.posterior_sigma_vec_print.T[0:stop])
    prior_sigma = data.options.heritability / (data.options.total_num_SNPs * data.options.p)
    inf_sigma = data.options.heritability / data.options.total_num_SNPs
    plt.title("=TRACE= Sigma inf: {:.2e}, Prior sigma: {:.2e}. "
              "Posterior sigma mean {:.2e}".format(inf_sigma, prior_sigma, data.posterior_sigma_vec_print.T[0:stop].mean()),
              fontsize=8)
    plt.subplot(4, 1, 3)
    plt.plot(data.window_mu_vec_print.T[0:stop])
    plt.title("=TRACE= Prior window: Random. "
              "Posterior window mean {:.2e}".format(data.window_mu_vec_print.T[0:stop].mean()),
              fontsize=8)
    plt.subplot(4, 1, 4)
    plt.plot(data.window_s_vec_print.T[0:stop])
    plt.title("=TRACE= Prior decay: {:.2e}. "
              "Posterior decay mean {:.2e}".format(data.prior_window_decay, data.window_s_vec_print.T[0:stop].mean()),
              fontsize=8)
    plt.tight_layout(pad=0.5)
    plt.figure(5)
    plt.plot(data.pi_avg_print[0:data.gibbs_iter])
    plt.title("=AVG= Gene posterior prob all chains")
    # self.plot_genes_pos(window_mu_sum / gibbs_iter)
    plt.figure(6)
    plt.subplot(4, 1, 1)
    gene = list(data.genes_extension_index.keys())[0]
    l = data.genes_extension_index[gene]
    plt.plot(data.pos_prob_pi_trace[:, l, :].T[0:stop])
    plt.axhline(y=0.05, color="grey", linestyle="--", linewidth=1)
    plt.title("Trace of gene {} on index 0 in all chains".format(gene))
    plt.subplot(4, 1, 2)
    plt.plot(data.within.T[0:stop])
    plt.title("Within chain variance All genes")
    plt.subplot(4, 1, 3)
    plt.plot(data.between.T[0:stop])
    plt.title("Between chain variance All genes")
    plt.subplot(4, 1, 4)
    plt.plot(data.ratio_var.T[1:stop])
    plt.title("Ratio variance (PSRF) All genes")
    plt.tight_layout(pad=0.5)
    plt.figure(7)
    plt.subplot(3, 2, 1)
    plt.scatter(snp_loc_to_plot, marginal_betas_to_plot, s=2)
    plt.title("Marginal betas", fontsize=8)
    plt.subplot(3, 2, 2)
    plt.scatter(snp_loc_to_plot, -np.log10(marginal_p_values_to_plot), s=2)
    plt.title("Marginal p-values", fontsize=8)
    plt.subplot(3, 2, 3)
    plt.scatter(snp_loc_to_plot, norm_betas_to_plot, s=2)
    plt.title("Norm betas", fontsize=8)
    plt.subplot(3, 2, 4)
    plt.scatter(snp_loc_to_plot, -np.log10(norm_p_values_to_plot), s=2)
    plt.title("norm p-values", fontsize=8)
    plt.subplot(3, 2, 5)
    plt.scatter(snp_loc_to_plot, to_plot, s=2)
    plt.title("Gibbs betas", fontsize=8)
    plt.subplot(3, 2, 6)
    plt.scatter(snp_loc_to_plot, -np.log10(gibbs_p_values_to_plot), s=2)
    plt.title("Gibbs p-values", fontsize=8)
    # plt.figure(7)
    # for l in range(data.num_genes):
    #     fig = plt.subplot(data.num_genes, 1, l+1)
    #     fig.plot(data.trace_sampled_G[l, 0:stop], linewidth=3)
    #     plt.grid(True)
    #     fig.xaxis.set_ticklabels([])
    #     fig.yaxis.set_ticklabels([])
    #     fig.xaxis.set_ticks_position('none')
    plt.show()


# @jit(nopython=True, parallel=True)  # Optimized for: ONE-CHAIN, NUMBA
def logistic_regression(data, points, max_iter=100, threshold=1e-4):
    points = points.reshape((data.shape[0],))
    w = np.asarray([0, 0], dtype=np.float64)  # initial value for the parameters
    x = np.concatenate((np.ones((data.shape[0], 1)), data), axis=1)
    y = sigmoid_f(x.dot(w))  # sigmoid_f(np.matmul(x, w))
    log_ll_old = np.sum(points * np.log(y) + ((1 - points) * np.log(1 - y)))
    for i in range(max_iter):
        R = y * (1 - y)
        # check if matrix is singular
        m = (x.T * R).dot(x)
        if np.linalg.det(m) == 0:
            m = matrix_stabilization(m)
        if np.linalg.det(m) == 0:
            break
        w = w - np.linalg.inv(m).dot(x.T).dot(y - points)
        y = sigmoid_f(x.dot(w))
        log_ll_new = np.sum(points * np.log(y) + ((1 - points) * np.log(1 - y)))
        tot = log_ll_new - log_ll_old
        log_ll_old = log_ll_new
        if tot < threshold:
            break
    center = -w[0] / w[1]
    decay = 1 / w[1]
    return center, decay


def update_window(param):
    link, distance, chain = param
    max_iter = 100
    threshold = 1e-4
    num_genes = distance.shape[0]
    num_snps = distance.shape[1]
    regression_x = np.empty(0)
    regression_y = np.empty(0)
    for l in range(num_genes):
        gene_X = np.concatenate((distance[l], np.zeros(num_snps)))
        gene_y = np.concatenate((link[chain, :, l], np.ones(num_snps)))
        regression_x = np.concatenate((regression_x, gene_X))
        regression_y = np.concatenate((regression_y, gene_y))
    points = regression_y.reshape((regression_y.shape[0],))
    w = np.asarray([0, 0])  # initial value for the parameters
    regression_x = regression_x.reshape(-1, 1)
    x = np.concatenate((np.ones((regression_x.shape[0], 1)), regression_x), axis=1)
    y = 1 / (1 + np.exp(-x.dot(w)))
    log_ll_old = np.sum(points * np.log(y) + ((1 - points) * np.log(1 - y)))
    for i in range(max_iter):
        R = y * (1 - y)
        w = w - np.linalg.inv((x.T * R).dot(x)).dot(x.T).dot(y - points)
        y = 1 / (1 + np.exp(-x.dot(w)))
        log_ll_new = np.sum(points * np.log(y) + ((1 - points) * np.log(1 - y)))
        tot = log_ll_new - log_ll_old
        log_ll_old = log_ll_new
        if tot < threshold:
            break
    center = -w[0] / w[1]
    decay = 1 / w[1]
    return [center, decay]


# THIS FUNCTION IS NOT WORKING CORRECTLY
# @jit(nopython=True, parallel=True)  # Optimized for: MULTI-CHAIN, NUMBA
def logistic_regression_mc(x, points, max_iter=100, threshold=1e-4):
    num_chains = x.shape[0]
    converge = np.zeros(num_chains, dtype=np.bool_)
    w = np.zeros((num_chains, 2))  # initial value for the parameters
    y = sigmoid_f(((x.T * w[:, 1]) + w[:, 0]).T)
    log_ll_old = np.sum(points * np.log(y) + ((1 - points) * np.log(1 - y)), axis=1)
    for i in range(max_iter):
        R = y * (1 - y)
        a = R.sum(axis=1)
        b = (x * R).sum(axis=1)
        d = ((x ** 2) * R).sum(axis=1)
        zed = (y - points)
        z = zed.sum(axis=1)
        zx = (zed * x).sum(axis=1)
        A = 1 / (a * d - (b ** 2))
        w[:, 0] = w[:, 0] - (A * ((d * z) - (b * zx)))
        w[:, 1] = w[:, 1] - (A * ((a * zx) - (b * z)))
        y = sigmoid_f(((x.T * w[:, 1]) + w[:, 0]).T)
        log_ll_new = np.sum(points * np.log(y) + ((1 - points) * np.log(1 - y)), axis=1)
        tot = log_ll_new - log_ll_old
        # print("old: ", log_ll_old)
        # print("new: ", log_ll_new)
        # print("tot: ", tot)
        log_ll_old = log_ll_new
        converge = converge + (tot < threshold)
        if converge.all():
            break
    center = -w[:, 0] / w[:, 1]
    decay = 1 / w[:, 1]
    return center, decay


# @jit(nopython=True, parallel=True)  # Optimized for: ONE-CHAIN, NUMBA
def matrix_stabilization(matrix):
    # stabilization achieved by multiply the diagonal by a small factor to avoid singular matrix
    matrix_size = matrix.shape[0]  # this is a square matrix
    null_cov_diag = np.diag(matrix)
    stab_term_frac = 1e-4
    stab_matrix = np.diag(np.ones(matrix_size)) * (np.mean(np.diag(matrix)) * stab_term_frac)
    null_cov_matrix = matrix + stab_matrix
    # eigen value decomposition to assure the resulting correlation matrix is positive semi-definite
    try:
        l, Q = np.linalg.eigh(null_cov_matrix)
    except np.linalg.linalg.LinAlgError:
        return null_cov_matrix
    l[l < 0] = 0
    null_cov_matrix = np.dot(Q, np.dot(np.diag(l), Q.transpose())) + stab_matrix
    # transformation to restore the values changed in the stabilization step
    null_cov_norm_matrix = np.diag(np.sqrt(null_cov_diag / np.diag(null_cov_matrix)))
    null_cov_matrix = np.dot(null_cov_norm_matrix, np.dot(null_cov_matrix, null_cov_norm_matrix))
    return null_cov_matrix


# @jit(nopython=True, parallel=True)  # Optimized for: ONLY NUMBA
def update_link_func_numba(num_chains, num_genes, num_snps, I_in, sampled_G_in, sampled_c, p_in, p_out,
                           snp_to_gene_dist, window_mu, window_s):
    I = I_in.copy().astype(np.float64)
    sampled_G = sampled_G_in.copy().astype(np.float64)
    random_gene_list = np.asarray(list(range(num_genes)))
    np.random.shuffle(random_gene_list)
    random_gene_list = list(random_gene_list)
    # random.shuffle(random_gene_list)
    for l in random_gene_list:
        I[:, :, l] = 1
        indicator_all = np.zeros((num_chains, num_snps), dtype=np.float64)
        for chain in range(num_chains):
            indicator_all[chain] = I[chain].dot(sampled_G[chain]).reshape((-1,))
        indicator_all = indicator_all.astype(np.bool_)
        indicator_all = indicator_all.astype(np.float64)
        c_likelihood_1 = link_likelihood_func(indicator_all, sampled_c, p_in, p_out)
        I[:, :, l] = 0
        for chain in range(num_chains):
            indicator_all[chain] = I[chain].dot(sampled_G[chain]).reshape((-1,))
        indicator_all = indicator_all.astype(np.bool_)
        indicator_all = indicator_all.astype(np.float64)
        c_likelihood_0 = link_likelihood_func(indicator_all, sampled_c, p_in, p_out)
        dist_to_l = np.repeat(snp_to_gene_dist[l, :].reshape(-1, 1), num_chains).reshape((num_snps, num_chains))
        exponent = (-(dist_to_l - window_mu) / window_s).T
        prior_p_I_l_j_1 = (1 / (1 + np.exp(exponent)))  # * 0.9875  # reduce the strength of prior
        link_numerator = np.multiply(c_likelihood_1, prior_p_I_l_j_1)
        link_denominator = np.multiply(c_likelihood_0, (1 - prior_p_I_l_j_1)) + link_numerator
        posterior_p_I_l_j_1 = np.divide(link_numerator, link_denominator)
        # I[:, :, l] = np.random.binomial(1, posterior_p_I_l_j_1)
        for chain in range(num_chains):
            for j in range(num_snps):
                I[chain, j, l] = np.random.binomial(1, posterior_p_I_l_j_1[chain, j])
    return I.copy().astype(np.bool_)


# @jit(nopython=True, parallel=True)  # Optimized for: MULTI-CHAIN, NUMBA
def link_likelihood_func(indicator_all, sampled_c, p_in, p_out):
    aux1 = np.multiply((1 - sampled_c).T, (1 - p_in)).T
    aux2 = np.multiply(sampled_c.T, p_in).T
    aux3 = np.multiply((1 - sampled_c).T, (1 - p_out)).T
    aux4 = np.multiply(sampled_c.T, p_out).T
    aux5 = np.multiply(indicator_all, aux2 + aux1)
    aux6 = np.multiply((1 - indicator_all), aux4 + aux3)
    c_likelihood = aux5 + aux6
    return c_likelihood


# Optimized for: ONLY MULTI-CHAIN
def update_link_func(num_chains, num_genes, num_snps, I, sampled_G, sampled_c, p_in, p_out,
                     snp_to_gene_dist, window_mu, window_s):
    random_gene_list = list(range(num_genes))
    random.shuffle(random_gene_list)
    for l in random_gene_list:
        I[:, :, l] = 1
        indicator_all = np.matmul(I, sampled_G).reshape((num_chains, num_snps))
        c_likelihood_1 = link_likelihood_func(indicator_all, sampled_c, p_in, p_out)
        I[:, :, l] = 0
        indicator_all = np.matmul(I, sampled_G).reshape((num_chains, num_snps))
        c_likelihood_0 = link_likelihood_func(indicator_all, sampled_c, p_in, p_out)
        dist_to_l = np.repeat(snp_to_gene_dist[l, :].reshape(-1, 1), num_chains).reshape((num_snps, num_chains))
        exponent = (-(dist_to_l - window_mu) / window_s).T
        prior_p_I_l_j_1 = (1 / (1 + np.exp(exponent)))  # * 0.9985  # reduce the strength of prior
        link_numerator = np.multiply(c_likelihood_1, prior_p_I_l_j_1)
        link_denominator = np.multiply(c_likelihood_0, (1 - prior_p_I_l_j_1)) + link_numerator
        posterior_p_I_l_j_1 = np.divide(link_numerator, link_denominator)
        posterior_p_I_l_j_1[np.isnan(posterior_p_I_l_j_1)] = 0  # change all NANs to 0
        I[:, :, l] = np.random.binomial(1, posterior_p_I_l_j_1)
    return I


# @jit(nopython=True, parallel=True)  # Optimized for: ONE-CHAIN, NUMBA
def allele_correction_2(snp_zs, snp_betas, snp_ns, null_cov_matrix, detect_flips):
    num_snps = len(snp_zs)
    if num_snps > 1:
        sorted_idx_by_n = np.argsort(snp_ns)
        for i in sorted_idx_by_n:
            # calculate the conditional distribution
            sigma11 = null_cov_matrix[i, i]
            sigma12 = np.concatenate((null_cov_matrix[i, :i], null_cov_matrix[i, (i + 1):]))
            zs12 = np.concatenate((snp_zs[:i], snp_zs[i + 1:]))
            # sigma22 = np.delete(np.delete(null_cov_matrix, i, axis=0), i, axis=1)
            condition = np.ones((num_snps, num_snps), dtype=np.float64)
            condition[:, i] = np.zeros(num_snps, dtype=np.float64)
            condition[i, :] = np.zeros(num_snps, dtype=np.float64)
            sigma22 = np.extract(condition, null_cov_matrix).reshape((num_snps - 1, num_snps - 1))
            sigma12_dot_sigma22_inv = sigma12.dot(np.linalg.inv(sigma22))
            conditional_mean = sigma12_dot_sigma22_inv.dot(zs12)
            conditional_cov = sigma11 - sigma12_dot_sigma22_inv.dot(sigma12)
            conditional_lik = numpy_logpdf(snp_zs[i], loc=conditional_mean, scale=np.sqrt(conditional_cov))
            conditional_flip_lik = numpy_logpdf(-snp_zs[i], loc=conditional_mean,
                                                scale=np.sqrt(conditional_cov))
            if detect_flips is not None and conditional_flip_lik - conditional_lik > detect_flips:
                snp_zs[i] = -snp_zs[i]
                snp_betas[i] = -snp_betas[i]
                # conditional_lik = conditional_flip_lik  # this update is not used
    return snp_zs, snp_betas


# @jit(nopython=True)  # Optimized for: ONLY NUMBA
def numpy_logpdf(x, loc, scale):
    var = scale ** 2
    denom = (2 * math.pi * var) ** .5
    num = np.exp(-(x - loc) ** 2 / (2 * var))
    return np.log(num / denom)


# NOT Optimized
def update_gene_with_plot(num_chains, num_snps, genes_extension, genes_extension_index, I, sampled_G, sampled_c,
                          sampled_beta, p_in, p_out, pi_prob, pos_prob_pi, snp_to_chr_pos, snps_in_region_index):
    r = list(genes_extension)
    np.random.shuffle(r)
    for gene in r:
        l = genes_extension_index[gene]
        # ==================================================================================================
        casual_genes_names = []
        plot_gene_update = False
        # causal genes
        if plot_gene_update:
            for gene_loop in genes_extension:
                l_loop = genes_extension_index[gene_loop]
                if sampled_G[0, l_loop]:
                    casual_genes_names.append(gene_loop)
        # ==================================================================================================
        sampled_G[:, l] = 1  # setting the current gene to one for all chains
        indicator_all_1 = np.matmul(I, sampled_G).reshape((num_chains, num_snps))
        snp_prob_all_1 = np.multiply(indicator_all_1,
                                     np.multiply(sampled_c, p_in[:, np.newaxis]) + np.multiply(
                                         (1 - sampled_c), (1 - p_in)[:, np.newaxis])) + \
                         np.multiply((1 - indicator_all_1),
                                     np.multiply(sampled_c, p_out[:, np.newaxis]) + np.multiply(
                                         (1 - sampled_c), (1 - p_out)[:, np.newaxis]))
        # multiplication_1 = snp_prob_all_1.prod(axis=1)
        x_1 = np.log(snp_prob_all_1).sum(axis=1) + np.log(pi_prob)
        sampled_G[:, l] = 0  # setting the current gene to zero for all chains
        indicator_all_0 = np.matmul(I, sampled_G).reshape((num_chains, num_snps))
        # ==================================================================================================
        if plot_gene_update:
            # sorting snps to plot by genomic position
            plt.figure(12)
            plt.clf()
            sort_snps = dict(sorted(snp_to_chr_pos.items(), key=lambda item: item[1]))
            to_plot = []
            snp_loc_to_plot = []
            to_plot_gene = []
            snp_loc_to_plot_gene = []
            to_plot_taken = []
            snp_loc_to_plot_taken = []
            selector = np.bitwise_and(indicator_all_1[0], indicator_all_0[0])
            for snp in sort_snps:
                if snp in snps_in_region_index:
                    idx = snps_in_region_index[snp]
                    if selector[idx]:
                        snp_loc_to_plot_taken.append(snp_to_chr_pos[snp][1])
                        to_plot_taken.append(sampled_beta[0, idx])
                    elif indicator_all_1[0, idx]:
                        snp_loc_to_plot_gene.append(snp_to_chr_pos[snp][1])
                        to_plot_gene.append(sampled_beta[0, idx])
                    else:
                        snp_loc_to_plot.append(snp_to_chr_pos[snp][1])
                        to_plot.append(sampled_beta[0, idx])
            # taken SNPs
            plt.scatter(snp_loc_to_plot_taken, to_plot_taken, s=2, color='r', marker='v')
            # outside gene SNPs
            if len(snp_loc_to_plot) != 0:
                plt.scatter(snp_loc_to_plot, to_plot, s=2, color='b', marker='^')
            # inside gene SNPs
            if len(snp_loc_to_plot_gene) != 0:
                plt.scatter(snp_loc_to_plot_gene, to_plot_gene, s=2, color='g', marker='^')
            plt.title("Available SNPs for gene: {}. Causal genes {}".format(gene, casual_genes_names),
                      fontsize=15)
            plt.show()
        # ==================================================================================================
        snp_prob_all_0 = np.multiply(indicator_all_0,
                                     np.multiply(sampled_c, p_in[:, np.newaxis]) + np.multiply(
                                         (1 - sampled_c), (1 - p_in)[:, np.newaxis])) + \
                         np.multiply((1 - indicator_all_0),
                                     np.multiply(sampled_c, p_out[:, np.newaxis]) + np.multiply(
                                         (1 - sampled_c), (1 - p_out)[:, np.newaxis]))
        # multiplication_0 = snp_prob_all_0.prod(axis=1)
        x_0 = np.log(snp_prob_all_0).sum(axis=1) + np.log(1 - pi_prob)
        c = np.max(np.column_stack((x_0, x_1)), axis=1)
        LSE = c + np.log(np.exp(x_1 - c) + np.exp(x_0 - c))
        LSE_prob = np.exp(x_1 - LSE)
        pos_prob_pi[:, l] = LSE_prob
        sampled_G[:, l] = np.random.binomial(1, LSE_prob, num_chains).reshape((num_chains, 1))
    return sampled_G


# Optimized for: ONLY MULTI-CHAIN
def update_gene(num_chains, num_snps, genes_extension, genes_extension_index, I, sampled_G, sampled_c,
                p_in, p_out, pi_prob, pos_prob_pi):
    r = list(genes_extension)
    np.random.shuffle(r)
    for gene in r:
        l = genes_extension_index[gene]
        # ==================================================================================================
        sampled_G[:, l] = 1  # setting the current gene to one for all chains
        indicator_all = np.matmul(I, sampled_G).reshape((num_chains, num_snps))
        snp_prob_all_1 = link_likelihood_func(indicator_all, sampled_c, p_in, p_out)
        x_1 = np.log(snp_prob_all_1).sum(axis=1) + np.log(pi_prob)
        # ==================================================================================================
        sampled_G[:, l] = 0  # setting the current gene to zero for all chains
        indicator_all = np.matmul(I, sampled_G).reshape((num_chains, num_snps))
        snp_prob_all_0 = link_likelihood_func(indicator_all, sampled_c, p_in, p_out)
        x_0 = np.log(snp_prob_all_0).sum(axis=1) + np.log(1 - pi_prob)
        # ==================================================================================================
        c = np.max(np.column_stack((x_0, x_1)), axis=1)
        LSE = c + np.log(np.exp(x_1 - c) + np.exp(x_0 - c))
        LSE_prob = np.exp(x_1 - LSE)
        pos_prob_pi[:, l] = LSE_prob
        sampled_G[:, l] = np.random.binomial(1, LSE_prob, num_chains).reshape((num_chains, 1))
    return sampled_G


# @jit(nopython=True, parallel=True)  # Optimized for: ONLY NUMBA
def update_gene_numba(num_chains, num_genes, num_snps, I_in, sampled_G_in, sampled_c, p_in, p_out,
                      pi_prob, pos_prob_pi):
    I = I_in.copy().astype(np.float64)
    sampled_G = sampled_G_in.copy().astype(np.float64)
    random_gene_list = np.asarray(list(range(num_genes)))
    np.random.shuffle(random_gene_list)
    random_gene_list = list(random_gene_list)
    for l in random_gene_list:
        # ==================================================================================================
        sampled_G[:, l] = 1  # setting the current gene to one for all chains
        # indicator_all_1 = np.matmul(I, sampled_G).reshape((num_chains, num_snps))
        indicator_all = np.zeros((num_chains, num_snps), dtype=np.float64)
        for chain in range(num_chains):
            indicator_all[chain] = I[chain].dot(sampled_G[chain]).reshape((-1,))
        indicator_all = indicator_all.astype(np.bool_)
        indicator_all = indicator_all.astype(np.float64)
        snp_prob_all_1 = link_likelihood_func(indicator_all, sampled_c, p_in, p_out)
        # snp_prob_all_1 = np.multiply(indicator_all_1,
        #                              np.multiply(sampled_c, p_in[:, np.newaxis]) + np.multiply(
        #                                  (1 - sampled_c), (1 - p_in)[:, np.newaxis])) + \
        #                  np.multiply((1 - indicator_all_1),
        #                              np.multiply(sampled_c, p_out[:, np.newaxis]) + np.multiply(
        #                                  (1 - sampled_c), (1 - p_out)[:, np.newaxis]))
        # multiplication_1 = snp_prob_all_1.prod(axis=1)
        x_1 = np.log(snp_prob_all_1).sum(axis=1) + np.log(pi_prob)
        sampled_G[:, l] = 0  # setting the current gene to zero for all chains
        for chain in range(num_chains):
            indicator_all[chain] = I[chain].dot(sampled_G[chain]).reshape((-1,))
        indicator_all = indicator_all.astype(np.bool_)
        indicator_all = indicator_all.astype(np.float64)
        # indicator_all_0 = np.matmul(I, sampled_G).reshape((num_chains, num_snps))
        # ==================================================================================================
        snp_prob_all_0 = link_likelihood_func(indicator_all, sampled_c, p_in, p_out)
        # snp_prob_all_0 = np.multiply(indicator_all_0,
        #                              np.multiply(sampled_c, p_in[:, np.newaxis]) + np.multiply(
        #                                  (1 - sampled_c), (1 - p_in)[:, np.newaxis])) + \
        #                  np.multiply((1 - indicator_all_0),
        #                              np.multiply(sampled_c, p_out[:, np.newaxis]) + np.multiply(
        #                                  (1 - sampled_c), (1 - p_out)[:, np.newaxis]))
        # multiplication_0 = snp_prob_all_0.prod(axis=1)
        x_0 = np.log(snp_prob_all_0).sum(axis=1) + np.log(1 - pi_prob)
        # c1 = np.max(np.column_stack((x_0, x_1)), axis=1)
        c = np.zeros(num_chains)
        vectors = np.column_stack((x_0, x_1))
        for chain in range(num_chains):
            c[chain] = np.max(vectors[chain])
        LSE = c + np.log(np.exp(x_1 - c) + np.exp(x_0 - c))
        LSE_prob = np.exp(x_1 - LSE)
        pos_prob_pi[:, l] = LSE_prob
        for chain in range(num_chains):
            sampled_G[chain, l] = np.random.binomial(1, LSE_prob[chain])
        # sampled_G[:, l] = np.random.binomial(1, LSE_prob, num_chains).reshape((num_chains, 1))
    return sampled_G.copy().astype(np.bool_)


def plot_window_init():
    plt.ion()
    fig, ax = plt.subplots()
    x, y = [], []
    plot_scatter = ax.scatter(x, y)
    plot_line, = ax.plot(x, y)
    plt.xlim(0, 10)
    plt.ylim(0, 1)
    plt.axhline(y=0.5, color="grey", linestyle="--", linewidth=1)
    plt.draw()
    plt.show()
    return plot_scatter, plot_line, fig


def plot_window_func(plot_scatter, plot_line, k, num_genes, snp_to_gene_dist, I, sampled_G, window_mu, window_s, fig):
    num_points = 1000
    if k == 0:  # we set the x axis the first time
        x_values = np.empty(0)
        for l in range(num_genes):
            x_values = np.concatenate((x_values, snp_to_gene_dist[l]))
        plt.xlim(min(x_values), max(x_values))
    regression_X = np.empty(0)
    regression_y = np.empty(0)
    for l in range(num_genes):
        # if the gene is causal
        if sampled_G[-1, l, 0]:
            gene_X = snp_to_gene_dist[l]
            gene_y = I[-1, :, l]
            regression_X = np.concatenate((regression_X, gene_X))
            regression_y = np.concatenate((regression_y, gene_y))
    regression_X = regression_X.reshape(-1, 1)
    if regression_X.size != 0:
        x_axis = np.linspace(min(regression_X), max(regression_X), num_points)
        logit = 1 / (1 + np.exp(-(x_axis - window_mu[-1]) / window_s[-1])).reshape(num_points, )
        x = regression_X.reshape(regression_X.shape[0], ).tolist()
        y = regression_y.tolist()
        plot_scatter.set_offsets(np.c_[x, y])
        plot_line.set_data(x_axis, logit)
        fig.canvas.draw_idle()
        plt.title("window of chain 0")
        plt.pause(0.01)
        # plt.waitforbuttonpress()


def update_window_func(num_chains, num_genes, sampled_G, snp_to_gene_dist, I, window_mu, window_s, options):
    for chain in range(num_chains):
        regression_X = np.empty(0)
        regression_y = np.empty(0)
        for l in range(num_genes):
            # if the gene is causal is added to the regression
            if sampled_G[chain, l, 0]:
                gene_X = snp_to_gene_dist[l]
                gene_y = I[chain, :, l]
                # gene_X = np.concatenate((snp_to_gene_dist[l], np.zeros(num_snps)))
                # gene_y = np.concatenate((I[chain, :, l], np.ones(num_snps)))
                regression_X = np.concatenate((regression_X, gene_X))
                regression_y = np.concatenate((regression_y, gene_y))
        if regression_X.size != 0:  # If there is at least one causal gene the window is updated
            regression_X = regression_X.reshape(-1, 1)
            ll_center, s = logistic_regression(data=regression_X, points=regression_y)
            window_mu, window_s = window_restriction(ll_center, s, options, window_mu, window_s, chain)
        else:
            log("No causal genes in this iteration, window was not updated", options)
    return window_mu, window_s


def window_restriction(ll_center, s, options, window_mu, window_s, chain):
    if 0 < ll_center < 1000000:
        window_mu[chain] = ll_center
    else:
        log("Reached invalid windows values: {}, decay: {}, using last iteration values"
            .format(ll_center, s), options, DEBUG)
    if -100000 < s < -15000:
        window_s[chain] = s
    else:
        log("Reached invalid decay values: {}, decay: {}, using last iteration values"
            .format(ll_center, s), options, DEBUG)
    return window_mu, window_s


def update_window_mc_func(num_chains, num_genes, num_snps, snp_to_gene_dist, I, window_mu, window_s, options):
    regression_X = np.repeat(snp_to_gene_dist.T.reshape(num_snps * num_genes).reshape(-1, 1),
                             num_chains, axis=1).T
    regression_y = I.reshape((num_chains, num_snps * num_genes))
    regression_X = np.concatenate((regression_X, np.zeros((num_chains, num_snps * num_genes))), axis=1)
    regression_y = np.concatenate((regression_y, np.ones((num_chains, num_snps * num_genes))), axis=1)
    center, decay = logistic_regression_mc(x=regression_X, points=regression_y)
    for chain in range(num_chains):
        ll_center = center[chain]
        s = decay[chain]
        window_mu, window_s = window_restriction(ll_center, s, options, window_mu, window_s, chain)


def log(message, options, level=TRACE, end_char='\n'):
    log_fh = sys.stdout
    if level <= options.debug_level:
        log_fh.write("%s%s" % (message, end_char))
        log_fh.flush()
