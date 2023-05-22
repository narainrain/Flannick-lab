import math

import numpy as np
from numpy.random import RandomState
from scipy.stats import random_correlation
import scipy as sc
import pickle as pl
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('Qt5Agg')


def plot_snps(data, true_betas, marginal_betas, p_values, inf_beta, inf_p_values, gibbs_betas, gibbs_p_values):
    # sorting snps to plot by genomic position
    sort_snps = dict(sorted(data.snp_to_chr_pos.items(), key=lambda item: item[1]))
    snp_loc_to_plot = np.zeros(data.num_snps)
    to_plot = np.zeros(data.num_snps)
    marginal_betas_to_plot = np.zeros(data.num_snps)
    p_values_to_plot = np.zeros(data.num_snps)
    inf_betas_to_plot = np.zeros(data.num_snps)
    inf_p_values_to_plot = np.zeros(data.num_snps)
    gibbs_betas_to_plot = np.zeros(data.num_snps)
    gibbs_p_values_to_plot = np.zeros(data.num_snps)
    j = 0
    for snp in sort_snps:
        if snp in data.snps_in_region_index:
            snp_loc_to_plot[j] = data.snp_to_chr_pos[snp][1]
            to_plot[j] = true_betas[data.snps_in_region_index[snp]]
            marginal_betas_to_plot[j] = marginal_betas[data.snps_in_region_index[snp]]
            p_values_to_plot[j] = p_values[data.snps_in_region_index[snp]]
            inf_betas_to_plot[j] = inf_beta[data.snps_in_region_index[snp]]
            inf_p_values_to_plot[j] = inf_p_values[data.snps_in_region_index[snp]]
            inf_betas_to_plot[j] = inf_beta[data.snps_in_region_index[snp]]
            inf_p_values_to_plot[j] = inf_p_values[data.snps_in_region_index[snp]]
            gibbs_betas_to_plot[j] = gibbs_betas[data.snps_in_region_index[snp]]
            gibbs_p_values_to_plot[j] = gibbs_p_values[data.snps_in_region_index[snp]]
            j += 1
    plt.figure(1)
    plt.subplot(2, 2, 1)
    plt.scatter(snp_loc_to_plot, to_plot, s=2)
    plt.title("True betas", fontsize=8)
    plt.subplot(2, 2, 3)
    plt.scatter(snp_loc_to_plot, marginal_betas_to_plot, s=2)
    plt.title("Marginal betas", fontsize=8)
    plt.subplot(2, 2, 4)
    plt.scatter(snp_loc_to_plot, -np.log10(p_values_to_plot), s=2)
    plt.title("p-values", fontsize=8)
    plt.figure(2)
    plt.subplot(2, 2, 1)
    plt.scatter(snp_loc_to_plot, inf_betas_to_plot, s=2)
    plt.title("Inf betas", fontsize=8)
    plt.subplot(2, 2, 2)
    plt.scatter(snp_loc_to_plot, -np.log10(inf_p_values_to_plot), s=2)
    plt.title("Inf p-values", fontsize=8)
    plt.subplot(2, 2, 3)
    plt.scatter(snp_loc_to_plot, gibbs_betas_to_plot, s=2)
    plt.title("Gibbs betas", fontsize=8)
    plt.subplot(2, 2, 4)
    plt.scatter(snp_loc_to_plot, -np.log10(gibbs_p_values_to_plot), s=2)
    plt.title("Gibbs p-values", fontsize=8)
    plt.show()


# creation of correlation ===================================
file = open("../data/out/corr", 'rb')
corr = pl.load(file)
# corr = np.asarray(random_correlation.rvs((.5, .8, 1.2, 1.5, .5, .8, 1.2, 1.5, .8, 1.2), random_state=1))
M = corr.shape[0]
# creation of true betas =====================================
true_betas = np.zeros(M)
with open("../data/out/eg.data", 'rb') as f:
    data = pl.load(f)
p = 1/M
M_c = math.ceil(M * p)
half = int(math.floor(M_c/2))
causal_snp = data.snps_in_region_index["rs2803619"]
if (M_c % 2) != 0:
    true_betas[causal_snp - half:causal_snp + half + 1] = -0.03
else:
    true_betas[causal_snp - half:causal_snp + half] = -0.03
non_zero = np.nonzero(true_betas)
h_2 = (true_betas ** 2).sum() / M_c
# Parameters =====================================
n = 95454
sigma_2 = h_2 / (M*p)
sigma_2_e = 1 - h_2

real_marginals = True
if real_marginals:
    marginal_betas = pl.load(open("../data/out/marginal_beta", 'rb'))[0]
    snp_ses = pl.load(open("../data/out/snp_ses", 'rb'))
    # creation of z-scores  ======================================
    z_scores = marginal_betas / snp_ses  # * np.sqrt(n/sigma_2_e)
else:
    # creation of marginal betas =================================
    cov = corr * sigma_2_e / n
    mean = corr.dot(true_betas)
    prng = RandomState(200788)
    marginal_betas = prng.multivariate_normal(mean, cov, 1, )[0]
    # creation of z-scores  ======================================
    z_scores = marginal_betas * np.sqrt(n/sigma_2_e)

# creation of p-values =======================================
p_values = np.zeros(M)
for j in range(M):
    p_values[j] = sc.stats.norm.sf(abs(z_scores[j])) * 2
# ============================================================
# LDPRED ===============================================================
# ======== infinitesimal model =========================================
factor = np.diag(np.ones(M) * (M / (n * h_2)))
# calculation of infinitesimal betas
# normalization of betas
norm_beta = z_scores * np.sqrt(sigma_2_e/n)  # on real data not all SNPs have the same std.
inf_beta = np.linalg.inv(factor + corr).dot(norm_beta)
inf_z_scores = inf_beta * np.sqrt(n/sigma_2_e)
inf_p_values = np.zeros(M)
for j in range(M):
    inf_p_values[j] = sc.stats.norm.sf(abs(inf_z_scores[j])) * 2
# ======== Gibbs =======================================================
burn_in = 100
max_iter = 500
sampled_beta_vec = np.zeros((M, max_iter))
sampled_h_vec = np.zeros(max_iter)
sampled_p_vec = np.zeros(max_iter)
sampled_beta = norm_beta.copy()  # initial value
for gibbs_iter in range(max_iter):
    print("iter: ", gibbs_iter)
    for j in range(M):
        shrinkage = (n * sigma_2) / ((n * sigma_2) + 1)
        # calculation of beta tilde j (joint beta) ========================
        corr_no_j = np.delete(corr[j], j, 0)
        sampled_beta_no_j = np.delete(sampled_beta, j, 0)
        snp_residual_beta = norm_beta[j] - np.dot(sampled_beta_no_j, corr_no_j)
        # Compute pj and sample cj  =======================================
        snp_prob = (1 / (1 + (((1 - p) / p) *
                              np.sqrt(1 + (n * sigma_2)) *
                              np.exp((-1 / 2) * (n * np.power(snp_residual_beta, 2) * shrinkage)))
                         ))
        aux = np.random.random(1)[0]
        if aux <= snp_prob:
            mean_beta = snp_residual_beta * shrinkage
            variance = (1 / n) * shrinkage
            sampled_j = np.random.normal(mean_beta, np.sqrt(variance))
        else:
            sampled_j = 0
        sampled_beta[j] = sampled_j
        sampled_beta_vec[j, gibbs_iter] = sampled_j
    # UPDATE ==============================
    # M_c = np.count_nonzero(sampled_beta)
    # top_mc = sampled_beta[np.argsort(sampled_beta)[-M_c:]]
    # h_2 = (top_mc ** 2).sum() / M_c
    # # p = M_c / M
    # sigma_2 = h_2 / (M * p)
    # =====================================
    sampled_h_vec[gibbs_iter] = np.power(sampled_beta, 2).sum()
    sampled_p_vec[gibbs_iter] = np.count_nonzero(sampled_beta) / M
gibbs_betas = sampled_beta_vec[:, burn_in:max_iter].mean(axis=1)
gibbs_h = sampled_h_vec[burn_in:max_iter].mean()
gibbs_p = sampled_p_vec[burn_in:max_iter].mean()
# ======================================================================
gibbs_z_scores = gibbs_betas * np.sqrt(n/sigma_2_e)
gibbs_p_values = np.zeros(M)
for j in range(M):
    gibbs_p_values[j] = sc.stats.norm.sf(abs(gibbs_z_scores[j])) * 2
plot_snps(data, true_betas, marginal_betas, p_values, inf_beta, inf_p_values, gibbs_betas, gibbs_p_values)
np.set_printoptions(formatter={'float': lambda x: "{0:0.4f}".format(x)})
print("True betas     :", true_betas)
print("Gibbs betas    :", gibbs_betas)
print("Inf betas      :", inf_beta)
print("Marginal betas :", marginal_betas)
print("Gibbs h_2: ", gibbs_h)
print("gibbs p: ", gibbs_p)
# ======================================================================



