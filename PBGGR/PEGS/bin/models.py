import pegs_utils as pg
import sys
import os
import scipy as sc
import numpy as np
from sklearn.linear_model import LogisticRegression
import matplotlib.pyplot as plt
from itertools import cycle
import math
import warnings
import matplotlib
import pickle as pl


class PEGS:
    """class description"""
    options = None
    bad_region = False
    # file handlers
    posterior_betas_fh = None
    # constant variables
    NONE = 0
    INFO = 1
    DEBUG = 2
    TRACE = 3
    # variables with the collected data
    user_range = {}
    extended_range = {}
    component_rsid = {}
    component_snp_betas = {}
    component_correlation = {}
    component_snp_in_comp_index = {}
    genes_in_region = {}
    genes_extension = {}
    genes_position = {}
    component_to_snp = {}
    snp_to_component = {}
    snp_to_chr_pos = {}
    snp_to_z = {}
    snp_to_n = {}
    snp_to_se = {}
    snp_to_beta = {}
    snp_to_p = {}
    snp_to_effect_allele = {}
    snp_to_other_allele = {}
    snps_in_region = None
    snps_in_region_index = {}
    genes_extension_index = {}
    ld_dict = {}
    stored_data = None
    tabix_files_created = False

    def __init__(self, options):
        self.options = options
        self.set_range_restrictions()
        self.rename_files()
        self.create_tabix_files()
        # Loading data ==============================
        if self.options.tabix:
            if self.tabix_files_created:
                self.get_data()
        else:
            self.get_data()
        self.delete_tmp_files()

    def gibbs(self, burn_in=20, gibbs_iter=100, plot=False, plot_window_regression=False, regression_method="BFGS"):
        # ======== initialization ======================================================================================
        exit_status = "success"
        total_iter = burn_in + gibbs_iter
        epsilon = 0.0001
        pi = self.options.gene_prior
        num_snps = len(self.snps_in_region)
        num_genes = len(self.genes_extension)
        snp_rsid, snp_betas, snp_ses, snp_zs, snp_ns = self.extract_region_snps()
        norm_beta = snp_betas  # / se)  # * math.sqrt(n)  used when the marginal betas aren't normalized
        # initializing parameters with whole genome values
        M = self.options.total_num_SNPs
        p = self.options.p
        n = self.options.sample_size
        h_norm = self.options.heritability
        shrinkage = 1 / (1 + ((M * p) / (n * h_norm)))
        p_gene = p
        epsilon_gene = epsilon
        ldpred_h_norm = None
        # PRIORS DEFINITIONS ===========================================================================================
        # prior for window
        prior_window_center = self.options.causal_window
        prior_window_decay = -2000
        window_beta0_strength = 1000
        window_beta1_strength = 30000000000  # inverse of the variance
        # initial values for center and decay in case fail of convergence of optimization method
        mean = prior_window_center
        s = prior_window_decay
        # Prior for sigma
        prior_sigma = h_norm / (M * p)
        prior_alpha = 10000  # confidence on the prior sigma
        prior_beta = prior_sigma * prior_alpha
        posterior_sigma = prior_sigma  # initial value this will be updated
        # Priors for p
        prior_p_alpha = 1000
        prior_p_beta = (1 / p) * prior_p_alpha
        # Priors for p_in
        prior_p_in_alpha = 1
        prior_p_in_beta = (1 / (p * 10)) * prior_p_in_alpha  # the prior says is twice more
        # Priors for p_out
        prior_p_out_alpha = 1
        prior_p_out_beta = (1 / epsilon) * prior_p_out_alpha
        # ======== variables for SNP to gene ===========================================================================
        snp_to_gene_dist = self.get_snp_to_gene_dist()
        I = self.get_gene_to_snp_link()
        window_mu = np.ones(num_genes) * self.options.causal_window
        window_s = np.zeros(num_genes)
        window_mu_sum = np.zeros(num_genes)
        window_s_sum = np.zeros(num_genes)
        # ======== variables for plotting ==============================================================================
        window_mu_vec_print = np.zeros((total_iter, num_genes))
        window_s_vec_print = np.zeros((total_iter, num_genes))
        print_beta = np.zeros(total_iter)
        max_beta_idx = np.argmax(abs(norm_beta))
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
        # ======== infinitesimal model =================================================================================
        factor = np.diag(np.ones(num_snps) * (M / (n * h_norm)))
        region_correlation = self.get_region_correlation()
        inf_beta = np.matmul(np.linalg.inv(factor + region_correlation), norm_beta)
        # variables to update inside the gibbs loop
        pos_expected_beta = np.zeros(num_snps)
        pos_expected_beta_sum = np.zeros(num_snps)
        sampled_beta = np.copy(norm_beta)
        sampled_c = np.zeros(num_snps)
        pos_prob_c = np.zeros(num_snps)
        pos_prob_c_sum = np.zeros(num_snps)
        pos_prob_pi = np.zeros(num_genes)
        pos_prob_pi_sum = np.zeros(num_genes)
        # G is initialized so all genes are causal
        sampled_G = np.ones(num_genes)  # np.random.binomial(1, 0.05, num_genes)# is better to initialize all genes to 1
        # ======== aux variables =======================================================================================
        plot_scatter = None
        plot_line = None
        fig = None
        if plot_window_regression:
            plt.ion()
            fig, ax = plt.subplots()
            x, y = [], []
            plot_scatter = ax.scatter(x, y)
            plot_line, = ax.plot(x, y)
            plt.xlim(0, 10)
            plt.ylim(0, 1)
            plt.draw()
            plt.show()
        # ======= gibbs iterations =====================================================================================
        for k in range(burn_in + gibbs_iter):
            self.log("iter: {}".format(k), self.INFO)
            # ==UPDATE== SNP causal status, true betas
            for component in self.component_to_snp:
                snp_rsid, comp_norm_beta, comp_sample_beta, correlation = \
                    self.get_component_data(component, k, sampled_beta)
                comp_num_snps = len(comp_norm_beta)
                comp_residual_beta = np.zeros(comp_num_snps)
                comp_pos_prob_c = np.zeros(comp_num_snps)
                for j in range(comp_num_snps):
                    snp = snp_rsid[j]
                    # calculation of beta tilde j (joint beta) ========================
                    corr_no_j = np.delete(correlation[j], j, 0)
                    sampled_beta_no_j = np.delete(comp_sample_beta, j, 0)
                    comp_residual_beta[j] = comp_norm_beta[j] - np.dot(sampled_beta_no_j.T, corr_no_j)
                    # Compute pj and sample cj  ==================================
                    if pg.S2CG(j, I, sampled_G):
                        comp_pos_prob_c[j] = 1 / (1 + (((1 - p) / p) * math.sqrt(1 + (n * posterior_sigma))
                                                       * math.exp((-1 / 2) * ((n * comp_residual_beta[j] ** 2) /
                                                                              (1 + (1 / (n * posterior_sigma)))))))
                    else:
                        comp_pos_prob_c[j] = epsilon
                    pos_prob_c[self.snps_in_region_index[snp]] = comp_pos_prob_c[j]
                    sampled_c[self.snps_in_region_index[snp]] = np.random.binomial(1, comp_pos_prob_c[j], 1)[0]
                    # Sample bj  ======================================
                    mean = comp_residual_beta[j] * shrinkage
                    if sampled_c[self.snps_in_region_index[snp]]:
                        variance = (1 / n) * shrinkage
                        comp_sample_beta[j] = np.random.normal(mean, math.sqrt(variance), 1)
                        sampled_beta[self.snps_in_region_index[snp]] = comp_sample_beta[j]
                    else:
                        comp_sample_beta[j] = 0.0
                        sampled_beta[self.snps_in_region_index[snp]] = 0.0
                    pos_expected_beta[self.snps_in_region_index[snp]] = comp_pos_prob_c[j] * mean
            # ==UPDATE== gene causal status
            update_gene = True
            if update_gene:
                for l in range(num_genes):
                    sampled_G[l] = 1  # setting the current gene to one
                    multiplication_1 = 1
                    for j in range(num_snps):
                        if pg.S2CG(j, I, sampled_G):
                            if sampled_c[j] == 1:
                                multiplication_1 *= p_gene
                            else:
                                multiplication_1 *= (1 - p_gene)
                        else:
                            if sampled_c[j] == 1:
                                multiplication_1 *= epsilon_gene
                            else:
                                multiplication_1 *= (1 - epsilon_gene)
                    sampled_G[l] = 0  # setting the current gene to zero
                    multiplication_0 = 1
                    for j in range(num_snps):
                        if pg.S2CG(j, I, sampled_G):
                            if sampled_c[j] == 1:
                                multiplication_0 *= p_gene
                            else:
                                multiplication_0 *= (1 - p_gene)
                        else:
                            if sampled_c[j] == 1:
                                multiplication_0 *= epsilon_gene
                            else:
                                multiplication_0 *= (1 - epsilon_gene)
                    pi_denominator = ((multiplication_1 * pi) + (multiplication_0 * (1 - pi)))
                    if pi_denominator != 0:
                        pos_prob_pi[l] = (multiplication_1 * pi) / pi_denominator
                    else:
                        pos_prob_pi[l] = 0
                    sampled_G[l] = np.random.binomial(1, pos_prob_pi[l], 1)[0]
            # ==UPDATE== p, heritability, p_in, p_out
            update_p = True
            if update_p:
                # update gene probability
                L_c = np.count_nonzero(sampled_G)
                pi = np.random.beta(1 + L_c, 19 + num_genes - L_c)  # 1 and 19 set a week prior of 0.05
                # update of fraction of causal SNPs in the region
                M_c = np.count_nonzero(sampled_beta)
                p = np.random.beta(prior_p_alpha + M_c, prior_p_beta + num_snps - M_c)
                # update the number of SNPs in the region
                M = num_snps
                # calculation of fraction of causal SNPs p_in inside of causal genes
                M_c_in = 0
                M_c_out = 0
                snps_inside_causal_genes = 0
                for j in range(num_snps):
                    if pg.S2CG(j, I, sampled_G):  # inside the genes
                        snps_inside_causal_genes += 1
                        if sampled_beta[j] != 0:
                            M_c_in += 1
                    else:  # outside genes
                        if sampled_beta[j] != 0:
                            M_c_out += 1
                snps_outside_causal_genes = num_snps - snps_inside_causal_genes
                p_in = np.random.beta(prior_p_in_alpha + M_c_in, prior_p_in_beta + snps_inside_causal_genes - M_c_in)
                p_out = np.random.beta(prior_p_out_alpha + M_c_out,
                                       prior_p_out_beta + snps_outside_causal_genes - M_c_out)
                p_gene = p_in
                epsilon_gene = p_out
                # update of heritability
                ldpred_h_norm = np.matmul(np.matmul(sampled_beta.T, region_correlation), sampled_beta)
                pos_alpha = (prior_alpha + (M / 2))
                pos_beta = (prior_beta + (np.matmul(sampled_beta.T, sampled_beta) / 2))
                posterior_sigma = sc.stats.invgamma.rvs(a=pos_alpha, scale=pos_beta, size=1)[0]
                # heritability is no longer used in gibbs it's used to compare
                h_norm = posterior_sigma * M_c
                # =======================================================================
            # ==UPDATE== window and S2G link
            update_window = True
            if update_window:
                # updating window mean and s
                regression_X = np.empty(0)
                regression_y = np.empty(0)
                for l in range(num_genes):
                    # gene_X = np.concatenate((gene_to_snp_distance[gene],
                    # np.negative(gene_to_snp_distance[gene]))).reshape(-1, 1)
                    # gene_y = np.concatenate((I[gene], np.ones(snps_in_component)))
                    # gene_X = snp_to_gene_dist[gene]
                    # gene_y = I[gene]
                    gene_X = np.concatenate((snp_to_gene_dist[l], np.zeros(num_snps)))
                    gene_y = np.concatenate((I[l, :], np.ones(num_snps)))
                    regression_X = np.concatenate((regression_X, gene_X))
                    regression_y = np.concatenate((regression_y, gene_y))
                regression_X = regression_X.reshape(-1, 1)
                with warnings.catch_warnings():
                    warnings.filterwarnings('error')
                    try:
                        if self.options.window_model == "one_window":
                            clf = LogisticRegression(random_state=0, C=100).fit(regression_X, regression_y)
                            mean = - clf.intercept_ / clf.coef_
                            s = 1 / clf.coef_
                        if self.options.window_model == "bayesian_one_window":
                            mean, s = pg.bayesian_logistic_regression(data=regression_X, points=regression_y,
                                                                      prior_center=prior_window_center,
                                                                      prior_decay=prior_window_decay,
                                                                      strength_beta0=window_beta0_strength,
                                                                      strength_beta1=window_beta1_strength,
                                                                      regression_method=regression_method)
                        # print(mean)
                        # print(s)
                        if not (-15000 < s < -1500):
                            self.log("s reached the minimum possible s ={}, using last iteration values".format(s),
                                     self.DEBUG)
                            mean = window_mu[0]
                            s = window_s[0]
                        window_mu = float(mean) * np.ones(num_genes)  # put the same window in all genes
                        window_s = float(s) * np.ones(num_genes)
                    except Warning:
                        self.log("Window not updated, failed to converge, using previous value")
                if plot_window_regression:  # used to check the evolution of the window for a gene
                    num_points = 1000
                    x_axis = np.linspace(min(regression_X), max(regression_X), num_points)
                    logit = 1 / (1 + np.exp(-(x_axis - mean) / s)).reshape(num_points, )
                    x = regression_X.reshape(regression_X.shape[0], ).tolist()
                    y = regression_y.tolist()
                    plt.xlim(min(regression_X), max(regression_X))
                    plot_scatter.set_offsets(np.c_[x, y])
                    plot_line.set_data(x_axis, logit)
                    fig.canvas.draw_idle()
                    plt.pause(0.01)
                    # plt.waitforbuttonpress()
                # update I
                for l in range(num_genes):
                    for j in range(num_snps):
                        dis_j = snp_to_gene_dist[l][j]
                        gene_mu = window_mu[l]
                        gene_s = window_s[l]
                        exponent = np.float64(-(dis_j - gene_mu) / gene_s)
                        # exponent = np.float64(-(dis_j - window_mu[l]) / window_s[l])
                        if abs(exponent) > 100:
                            prior_p_I_l_j_1 = 0
                        else:
                            prior_p_I_l_j_1 = 1 / (1 + np.exp(exponent))
                        I[l][j] = 1
                        c_likelihood_1 = 1
                        if pg.S2CG(j, I, sampled_G):
                            if sampled_c[j] == 1:
                                c_likelihood_1 *= p_gene
                            else:
                                c_likelihood_1 *= (1 - p_gene)
                        else:
                            if sampled_c[j] == 1:
                                c_likelihood_1 *= epsilon_gene
                            else:
                                c_likelihood_1 *= (1 - epsilon_gene)
                        I[l][j] = 0
                        c_likelihood_0 = 1
                        if pg.S2CG(j, I, sampled_G):
                            if sampled_c[j] == 1:
                                c_likelihood_0 *= p_gene
                            else:
                                c_likelihood_0 *= (1 - p_gene)
                        else:
                            if sampled_c[j] == 1:
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
                        I[l][j] = np.random.binomial(1, posterior_p_I_l_j_1, 1)[0]
            # ======= Sums after burn-in for average calculations =============
            if k >= burn_in:  # used to calculate the results by eliminated the first iterations
                pos_expected_beta_sum += pos_expected_beta
                pos_prob_c_sum += pos_prob_c
                pos_prob_pi_sum += pos_prob_pi
                window_mu_sum += window_mu
                window_s_sum += window_s
                # setting the averages  =================
                pi_avg_print[k - burn_in] = pos_prob_pi_sum / (k - burn_in + 1)
                # shrinkage = shrinkage_sum / (k - burn_in + 1)
            # ======= Update of the traces ===================================
            pi_vec_print[k] = pos_prob_pi
            print_beta[k] = pos_expected_beta[max_beta_idx]
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
                self.log("chain collapsed", self.INFO)
                exit_status = "collapsed"
                break
        # === end of gibbs =============================================================================================
        # === calculation of posterior betas and gene prob  ============================================================
        posterior_beta = pos_expected_beta_sum / gibbs_iter
        posterior_gene_prob = pos_prob_pi_sum / gibbs_iter
        # ==== plotting results ========================================================================================
        if exit_status == "success" and plot:
            plt.figure(1)
            plt.subplot(2, 2, 1)
            plt.plot(print_beta, linewidth=1)
            plt.title("=TRACE= convergence of strongest SNP", fontsize=8)
            plt.subplot(2, 2, 2)
            plt.scatter(inf_beta, posterior_beta, s=2)
            plt.title("infinitesimal vs Gibbs results", fontsize=8)
            # sorting snps to plot by genomic position
            sort_snps = dict(sorted(self.snp_to_chr_pos.items(), key=lambda item: item[1]))
            to_plot = np.zeros(num_snps)
            snp_loc_to_plot = np.zeros(num_snps)
            j = 0
            for snp in sort_snps:
                if snp in self.snps_in_region_index:
                    snp_loc_to_plot[j] = self.snp_to_chr_pos[snp][1]
                    to_plot[j] = posterior_beta[self.snps_in_region_index[snp]]
                    j += 1
            plt.subplot(2, 2, 3)
            plt.scatter(snp_loc_to_plot, to_plot, s=2)
            plt.title("=AVG minus burn-in= gibbs results", fontsize=8)
            # plot of gene probabilities
            ax = plt.subplot(2, 2, 4)
            plt.title("=TRACE= gene posteriors", fontsize=8)
            NUM_COLORS = num_genes
            cm = plt.get_cmap('gist_rainbow')
            line_type = ["-", "--", "-.", ":"]
            linecycler = cycle(line_type)
            ax.set_prop_cycle('color', [cm(1. * i / NUM_COLORS) for i in range(NUM_COLORS)])
            for gene in self.genes_extension_index:
                l = self.genes_extension_index[gene]
                ax.plot(pi_vec_print[:, l], next(linecycler), label=gene)
            # plt.legend(loc='upper center', bbox_to_anchor=(0.5, 1.1), ncol=3, fancybox=True, shadow=True)
            plt.tight_layout(pad=0.5)
            # plotting other probabilities
            plt.figure(2)
            plt.subplot(3, 1, 1)
            plt.plot(p_vec_print)
            plt.title("=TRACE= Prior: {:.4f}. "
                      "fraction of causal SNPs p in the region, "
                      "mean {:.4f}".format(self.options.p, p_vec_print.mean()),
                      fontsize=8)
            plt.subplot(3, 1, 2)
            plt.plot(p_in_vec_print)
            plt.title("=TRACE= Prior: {:.4f}. "
                      "fraction of causal SNPs p_in given the gene is causal, "
                      "mean {:.4f}".format(self.options.p, p_in_vec_print.mean()),
                      fontsize=8)
            plt.subplot(3, 1, 3)
            plt.plot(p_out_vec_print)
            plt.title("=TRACE= Prior: {:.4f}. "
                      "fraction of causal SNPs p_out given the gene is not causal, "
                      "mean {:.4f}".format(epsilon, p_out_vec_print.mean()),
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
            plt.title("=TRACE= shrinkage , final: {:.4f}, used: {:.4f}".format(shrinkage_vec_print[total_iter - 1],
                                                                               shrinkage), fontsize=8)
            plt.subplot(3, 1, 2)
            plt.plot(posterior_sigma_vec_print)
            plt.title("=TRACE= Prior sigma: {:.8f}. "
                      "Posterior sigma mean {:.8f}".format(prior_sigma, posterior_sigma_vec_print.mean()),
                      fontsize=8)
            plt.subplot(3, 1, 3)
            plt.plot(window_mu_vec_print)
            plt.title("=TRACE= Prior window: {:.8f}. "
                      "Posterior window mean {:.8f}".format(self.options.causal_window, window_mu_vec_print.mean()),
                      fontsize=8)
            plt.tight_layout(pad=0.5)
            plt.figure(5)
            plt.plot(pi_avg_print)
            plt.title("=AVG= Gene posterior prob")
            # self.plot_genes_pos(window_mu_sum / gibbs_iter)
            plt.show()
        # saving betas results =========================================================================================
        if self.options.betas_out_file is not None:
            try:
                posterior_betas_fh = open(self.options.betas_out_file, 'w')
                line = "MarkerName\tchr\tpos\tref\talt\tposterior_beta\tposterior_prob\tzscore\tabs_pos\n"
                posterior_betas_fh.write(line)
                for snp in self.snps_in_region_index:
                    idx = self.snps_in_region_index[snp]
                    chr = self.snp_to_chr_pos[snp][0]
                    pos = self.snp_to_chr_pos[snp][1]
                    line = "{}\t{}\t{}\tA\tG\t{}\t{}\t{}\t{}\n".format(snp, chr, pos, posterior_beta[idx],
                                                                       pos_prob_c_sum[idx] / gibbs_iter,
                                                                       snp_zs[idx], abs(posterior_beta[idx]))
                    posterior_betas_fh.write(line)
                posterior_betas_fh.close()
            except ValueError:
                pg.bail("Failed to open out file")
        # saving genes results =========================================================================================
        if self.options.genes_out_file is not None:
            try:
                posterior_genes_fh = open(self.options.genes_out_file, 'w')
                line = "gene\tposterior_gene_prob\twindow\tdecay\n"
                posterior_genes_fh.write(line)
                for gene in self.genes_extension_index:
                    l = self.genes_extension_index[gene]
                    line = "{}\t{}\t{}\t{}\n".format(gene, posterior_gene_prob[l],
                                                     window_mu_sum[l] / gibbs_iter,
                                                     window_s_sum[l] / gibbs_iter)
                    posterior_genes_fh.write(line)
                posterior_genes_fh.close()
            except ValueError:
                pg.bail("Failed to open out file")
        if plot_window_regression:
            plt.waitforbuttonpress()
        return posterior_beta, exit_status

    def gibbs_multi_chain(self, num_chains=20, burn_in=50, max_iter=500, plot_window_regression=False,
                          store_data=True, multi_chain_regression=False, convergence_method="weighted",
                          convergence_thr=1.01, infinite=False, local_test=False, file_name=None,
                          prior_alpha_sigma=1000, prior_p_alpha=100, prior_p_in_alpha=1000, prior_p_out_alpha=100,
                          prior_pi_alpha=100):
        # ======== initialization ======================================================================================
        exit_status = "success"
        total_iter = burn_in + max_iter
        epsilon = 0.0001
        pi = np.ones(num_chains) * self.options.gene_prior
        num_snps = len(self.snps_in_region)
        num_genes = len(self.genes_extension)
        snp_rsid, snp_betas, snp_ses, snp_zs, snp_ns = self.extract_region_snps()
        converged = False
        # initializing parameters with whole genome values
        M = self.options.total_num_SNPs
        p = np.ones(num_chains) * self.options.p
        n = self.options.sample_size
        # ======== infinitesimal model =================================================================================
        factor = np.diag(np.ones(num_snps) * (M / (n * self.options.heritability)))
        region_correlation = self.get_region_correlation()
        region_correlation = self.matrix_stabilization(region_correlation)
        # correction of the allele using the correlation matrix
        if not local_test:
            self.log("Allele correction 2 correlation based, infinitesimal...", self.DEBUG)
            snp_zs, snp_betas = pg.allele_correction_2(snp_zs, snp_betas, snp_ns, region_correlation,
                                                       self.options.detect_flips)
        else:
            matplotlib.use('Qt5Agg')  # loading library to visualize results locally
        # normalization of betas
        norm_beta = (snp_betas / snp_ses) / np.sqrt(self.options.sample_size)
        # calculation of infinitesimal betas
        inf_beta = np.matmul(np.linalg.inv(factor + region_correlation), norm_beta)
        # ==============================================================================================================
        h_norm = np.ones(num_chains) * self.options.heritability
        p_in = np.ones(num_chains) * 0.005
        p_out = np.ones(num_chains) * 0.0001
        snp_to_gene_dist = self.get_snp_to_gene_dist()
        # PRIORS DEFINITIONS ===========================================================================================
        # WINDOW ================================
        initial_window_decay = -15000
        # Prior for sigma ======================================
        prior_sigma = self.options.heritability / (self.options.total_num_SNPs * self.options.p)
        prior_beta = prior_sigma * prior_alpha_sigma
        posterior_sigma = np.ones(num_chains) * prior_sigma  # initial value this will be updated
        # initial shrinkage changes after first iteration
        shrinkage = (n * posterior_sigma) / ((n * posterior_sigma) + 1)
        # Priors for p ========================================
        prior_p_beta = (1 / self.options.p) * prior_p_alpha
        # Priors for p_in =======================================
        prior_p_in_beta = (1 / 0.005) * prior_p_in_alpha  # the prior says is twice more
        # Priors for p_out =======================================
        prior_p_out_beta = (1 / 0.0001) * prior_p_out_alpha
        # Prior for pi
        prior_pi_beta = (1 / self.options.gene_prior) * prior_p_alpha  # this will create a week prior of 0.05
        # ======== variables for plotting ==============================================================================
        window_mu_vec_print = np.zeros((num_chains, total_iter))
        window_s_vec_print = np.zeros((num_chains, total_iter))
        print_beta = np.zeros(total_iter)
        max_beta_pos_prob = np.zeros(total_iter)
        max_beta_idx = np.argmax(abs(norm_beta))
        pi_vec_print = np.zeros((total_iter, num_genes))
        pi_avg_print = np.zeros((max_iter, num_genes))
        p_vec_print = np.zeros((num_chains, total_iter))
        p_in_vec_print = np.zeros((num_chains, total_iter))
        p_out_vec_print = np.zeros((num_chains, total_iter))
        heritability_vec_print = np.zeros((num_chains, total_iter))
        posterior_sigma_vec_print = np.zeros((num_chains, total_iter))
        ldpred_heritability_vec = np.zeros((num_chains, total_iter))
        pi_pos_vec = np.zeros((num_chains, total_iter))
        shrinkage_vec_print = np.zeros((num_chains, total_iter))
        ldpred_h_norm = np.zeros(num_chains)
        trace_sampled_G = np.zeros((num_chains, num_genes, total_iter))
        # variables to update inside the gibbs loop ====================================================================
        pos_expected_beta = np.zeros((num_chains, num_snps))
        pos_expected_beta_sum = np.zeros((num_chains, num_snps))
        sampled_c = np.zeros((num_chains, num_snps))
        pos_prob_c = np.zeros((num_chains, num_snps))
        pos_prob_c_sum = np.zeros((num_chains, num_snps))
        pos_prob_pi = np.zeros((num_chains, num_genes))
        pos_prob_pi_sum = np.zeros((num_chains, num_genes))
        pos_prob_pi_trace = np.zeros((num_chains, num_genes, total_iter))
        between = np.zeros((num_genes, total_iter))
        within = np.zeros((num_genes, total_iter))
        ratio_var = np.zeros((num_genes, total_iter))
        # G is initialized so all genes are causal
        sampled_G = np.ones((num_chains, num_genes, 1), dtype=np.bool_)  # np.random.binomial(1, 0.05, num_genes)
        # variables for SNP to gene
        I = np.ones((num_chains, num_snps, num_genes), dtype=np.bool_)
        sampled_beta = np.zeros((num_chains, num_snps))
        # aux1 = self.get_gene_to_snp_link()
        for chain in range(num_chains):
            # I[chain] = aux1.T
            sampled_beta[chain] = np.copy(inf_beta)
        # np.ones(num_chains) * self.options.causal_window
        window_mu = (np.random.randint(100, 400, num_chains) * 1000).astype(np.float64)
        window_s = np.ones(num_chains) * initial_window_decay
        window_mu_sum = np.zeros(num_chains)
        window_s_sum = np.zeros(num_chains)
        # ======== aux variables =======================================================================================
        plot_scatter = None
        plot_line = None
        fig = None
        if plot_window_regression:
            plot_scatter, plot_line, fig = pg.plot_window_init()
        # ======= UPDATE CONTROL =======================================================================================
        update_gene = True
        update_p = True
        update_window = True
        update_links = True
        # ======= gibbs iterations =====================================================================================
        k = -1
        global_iter = -1
        while True:
            k += 1
            global_iter += 1
            self.log("iter: {}".format(global_iter), self.INFO)
            indicator_all = np.matmul(I, sampled_G)
            # ==UPDATE== SNP causal status, true betas
            for component in self.component_to_snp:
                snp_rsid, comp_norm_beta, comp_sample_beta, correlation = \
                    self.get_component_data_for_multi_chain(component, global_iter, sampled_beta, local_test)
                for j in range(len(snp_rsid)):
                    snp = snp_rsid[j]
                    # calculation of beta tilde j (joint beta) ========================
                    corr_no_j = np.delete(correlation[j], j, 0)
                    sampled_beta_no_j = np.delete(comp_sample_beta, j, 1)
                    snp_residual_beta = comp_norm_beta[:, j] - np.dot(sampled_beta_no_j, corr_no_j)
                    # Compute pj and sample cj  =======================================
                    mean_beta = snp_residual_beta * shrinkage
                    snp_prob = (1 / (1 + (((1 - p) / p) *
                                          np.sqrt(1 + (n * posterior_sigma)) *
                                          np.exp((-1 / 2) * (n * np.power(snp_residual_beta, 2) * shrinkage)))
                                     ))
                    indicator = indicator_all[:, j].reshape((num_chains,))  # pg.S2CG_mc(j, I, sampled_G)
                    snp_pos_prob_c = (indicator * snp_prob) + ((1 - indicator) * epsilon)
                    pos_prob_c[:, self.snps_in_region_index[snp]] = snp_pos_prob_c
                    snp_sampled_c = np.random.binomial(1, snp_pos_prob_c, num_chains)
                    sampled_c[:, self.snps_in_region_index[snp]] = snp_sampled_c
                    # Sample bj  ======================================================

                    variance = (1 / n) * shrinkage
                    samples = np.zeros(num_chains)
                    for chain in range(num_chains):
                        samples[chain] = np.random.normal(mean_beta[chain], math.sqrt(variance[chain]))
                    sampled_j = (snp_sampled_c * samples) + ((1 - snp_sampled_c) * np.zeros(num_chains))
                    sampled_beta[:, self.snps_in_region_index[snp]] = sampled_j
                    comp_sample_beta[:, j] = sampled_j
                    pos_expected_beta[:, self.snps_in_region_index[snp]] = sampled_j  # snp_pos_prob_c * mean_beta
            # ==UPDATE== gene causal status
            if update_gene:
                sampled_G = pg.update_gene(num_chains, num_snps, self.genes_extension, self.genes_extension_index,
                                           I, sampled_G, sampled_c, p_in, p_out, pi, pos_prob_pi)
            # ==UPDATE== p, heritability, p_in, p_out, pi
            if update_p:
                # update the number of SNPs in the region
                # M = num_snps  # M is not used we use num_spns
                # update the number of causal genes and snps
                L_c = np.count_nonzero(sampled_G, axis=1)
                M_c = np.count_nonzero(sampled_beta, axis=1)
                # calculation of fraction of causal SNPs p_in inside of causal genes
                M_c_in = np.zeros(num_chains)
                M_c_out = np.zeros(num_chains)
                snps_inside_causal_genes = np.zeros(num_chains)
                indicator_all = np.matmul(I, sampled_G)
                for j in range(num_snps):
                    indicator = indicator_all[:, j].reshape((num_chains,))
                    for chain in range(num_chains):
                        if indicator[chain]:  # inside the genes
                            snps_inside_causal_genes[chain] += 1
                            if sampled_beta[chain, j] != 0:
                                M_c_in[chain] += 1
                        else:  # outside genes
                            if sampled_beta[chain, j] != 0:
                                M_c_out[chain] += 1
                snps_outside_causal_genes = num_snps - snps_inside_causal_genes
                for chain in range(num_chains):
                    # update gene probability
                    pi[chain] = np.random.beta(prior_pi_alpha + L_c[chain], prior_pi_beta + num_genes - L_c[chain])
                    # update of fraction of causal SNPs in the region
                    p[chain] = np.random.beta(prior_p_alpha + M_c[chain], prior_p_beta + num_snps - M_c[chain])
                    # update p_in
                    p_in[chain] = np.random.beta(prior_p_in_alpha + M_c_in[chain],
                                                 prior_p_in_beta + snps_inside_causal_genes[chain] - M_c_in[chain])
                    p_out[chain] = np.random.beta(prior_p_out_alpha + M_c_out[chain],
                                                  prior_p_out_beta + snps_outside_causal_genes[chain] - M_c_out[chain])
                    # update of heritability
                    ldpred_h_norm[chain] = np.abs(np.matmul(np.matmul(sampled_beta[chain].T, region_correlation),
                                                            sampled_beta[chain]))
                    if M_c[chain] != 0:
                        pos_alpha = (prior_alpha_sigma + (M_c[chain] / 2))
                        pos_beta = (prior_beta + ((ldpred_h_norm[chain]) / 2))
                        posterior_sigma[chain] = sc.stats.invgamma.rvs(a=pos_alpha, scale=pos_beta, size=1)[0]
                shrinkage = (n * posterior_sigma) / ((n * posterior_sigma) + 1)
                # =======================================================================
            # ==UPDATE==  I, S2G link
            if update_links:
                I = pg.update_link_func(num_chains, num_genes, num_snps, I, sampled_G, sampled_c, p_in, p_out,
                                        snp_to_gene_dist, window_mu, window_s)
            # ==UPDATE== window
            if update_window:
                if multi_chain_regression:
                    window_mu, window_s = pg.update_window_mc_func(num_chains, num_genes, num_snps, snp_to_gene_dist, I,
                                                                   window_mu, window_s, self.options)
                else:
                    window_mu, window_s = pg.update_window_func(num_chains, num_genes, sampled_G, snp_to_gene_dist, I,
                                                                window_mu, window_s, self.options)
            # ==PLOT WINDOW==
            if plot_window_regression:  # plotting last chain
                pg.plot_window_func(plot_scatter, plot_line, k, num_genes, snp_to_gene_dist, I, sampled_G, window_mu,
                                    window_s, fig)
            # ======= Sums after burn-in for average calculations =============
            if k >= burn_in:  # used to calculate the results by eliminated the first iterations
                pos_expected_beta_sum += pos_expected_beta
                pos_prob_c_sum += pos_prob_c
                pos_prob_pi_sum += pos_prob_pi
                window_mu_sum += window_mu
                window_s_sum += window_s
                # setting the averages all chains for gene probability ======================
                pi_avg_print[k - burn_in] = pos_prob_pi_sum.mean(axis=0) / (k - burn_in + 1)
            # ======= Update of the traces of first chain ===================================
            pi_vec_print[k] = pos_prob_pi[0]
            print_beta[k] = pos_expected_beta[0, max_beta_idx]
            max_beta_pos_prob[k] = pos_prob_c[0, max_beta_idx]
            # ======= Update of the traces for all chains ===================================
            pos_prob_pi_trace[:, :, k] = pos_prob_pi
            p_vec_print[:, k] = p
            heritability_vec_print[:, k] = h_norm
            posterior_sigma_vec_print[:, k] = posterior_sigma
            ldpred_heritability_vec[:, k] = ldpred_h_norm
            p_in_vec_print[:, k] = p_in
            p_out_vec_print[:, k] = p_out
            pi_pos_vec[:, k] = pi
            shrinkage_vec_print[:, k] = shrinkage
            window_mu_vec_print[:, k] = window_mu
            window_s_vec_print[:, k] = window_s
            trace_sampled_G[:, :, k] = sampled_G.reshape((num_chains, num_genes,))
            # ====== Convergence measure =====================================================
            if k > 0 and num_chains > 1:
                overall_mean = np.mean(pos_prob_pi_trace[:, :, 0: (k + 1)], axis=2).mean(axis=0)
                between_var = ((np.mean(pos_prob_pi_trace[:, :, 0: (k + 1)], axis=2) - overall_mean) ** 2).sum(axis=0) \
                              * ((k + 1) / (num_chains - 1))
                within_var = np.var(pos_prob_pi_trace[:, :, 0: (k + 1)], axis=2).mean(axis=0)
                unbiased = ((((k + 1) - 1) / (k + 1)) * within_var) + (((num_chains + 1) / (num_chains * (k + 1)))
                                                                       * between_var)
                ratio = np.sqrt(unbiased / within_var)
                ratio_w = (np.asarray([0 if val < 0 else val for val in (ratio - 1)]) * pos_prob_pi.mean(axis=0)) + 1
                between[:, k] = between_var
                within[:, k] = within_var
                if convergence_method == "weighted":
                    ratio_var[:, k] = ratio_w
                    if max(ratio_w) < convergence_thr:
                        converged = True
                else:
                    ratio_var[:, k] = ratio
                    if max(ratio) < convergence_thr:
                        converged = True
            # ====== Saving results  =========================================================
            if (((k + 1) % (max_iter + burn_in) == 0) and (k != 0)) or (converged and k > (3 * burn_in)):
                gibbs_iter = k - burn_in + 1
                inf_suffix = ""
                if file_name is not None:
                    inf_suffix += file_name
                if infinite:
                    inf_suffix += "_iter_{}".format(global_iter + 1)
                # === calculation of posterior =========================================================================
                posterior_beta = pos_expected_beta_sum / gibbs_iter
                posterior_gene_prob = pos_prob_pi_sum / gibbs_iter
                posterior_prob_c = pos_prob_c_sum / gibbs_iter
                posterior_window = window_mu_sum / gibbs_iter
                posterior_decay = window_s_sum / gibbs_iter
                posterior_gene_freq = trace_sampled_G[:, :, burn_in:burn_in + gibbs_iter].sum(axis=2) / gibbs_iter
                # === aggregation of multiple chains
                posterior_beta = np.mean(posterior_beta, axis=0)
                posterior_gene_prob = np.mean(posterior_gene_prob, axis=0)
                posterior_prob_c = np.mean(posterior_prob_c, axis=0)
                posterior_window = np.mean(posterior_window, axis=0)
                posterior_decay = np.mean(posterior_decay, axis=0)
                posterior_gene_freq = np.mean(posterior_gene_freq, axis=0)
                # storing data for latter plot or fast warmup for gibbs  ===============================================
                gd = GibbsData()
                gd.print_beta = print_beta
                gd.max_beta_pos_prob = max_beta_pos_prob
                gd.inf_beta = inf_beta
                gd.posterior_beta = posterior_beta
                gd.snp_to_chr_pos = self.snp_to_chr_pos
                gd.num_snps = num_snps
                gd.num_genes = num_genes
                gd.snps_in_region_index = self.snps_in_region_index
                gd.genes_extension_index = self.genes_extension_index
                gd.pi_vec_print = pi_vec_print
                gd.p_vec_print = p_vec_print
                gd.epsilon = epsilon
                gd.options = self.options
                gd.p_in_vec_print = p_in_vec_print
                gd.p_out_vec_print = p_out_vec_print
                gd.heritability_vec_print = heritability_vec_print
                gd.ldpred_heritability_vec = ldpred_heritability_vec
                gd.pi_pos_vec = pi_pos_vec
                gd.shrinkage_vec_print = shrinkage_vec_print
                gd.posterior_sigma_vec_print = posterior_sigma_vec_print
                gd.window_mu_vec_print = window_mu_vec_print
                gd.window_s_vec_print = window_s_vec_print
                gd.pi_avg_print = pi_avg_print
                gd.shrinkage = shrinkage
                gd.pos_prob_pi_trace = pos_prob_pi_trace
                gd.within = within
                gd.between = between
                gd.ratio_var = ratio_var
                gd.gibbs_iter = gibbs_iter
                gd.burn_in = burn_in
                gd.prior_window_decay = initial_window_decay
                gd.trace_sampled_G = trace_sampled_G[0, :, :]  # we only save first chain
                gd.marginal_betas = snp_betas
                gd.marginal_se = snp_ses
                #  === save data for latter use
                self.stored_data = gd
                if store_data:
                    with open(self.options.data_out_file + inf_suffix, 'wb') as f:
                        pl.dump(gd, f)
                    f.close()
                # saving betas results =================================================================================
                if self.options.betas_out_file is not None:
                    try:
                        posterior_betas_fh = open(self.options.betas_out_file + inf_suffix, 'w')
                        line = "MarkerName\tchr\tpos\tref\talt\tposterior_beta\tposterior_prob\tzscore\tabs_pos\n"
                        posterior_betas_fh.write(line)
                        for snp in self.snps_in_region_index:
                            idx = self.snps_in_region_index[snp]
                            chr = self.snp_to_chr_pos[snp][0]
                            pos = self.snp_to_chr_pos[snp][1]
                            line = "{}\t{}\t{}\tA\tG\t{}\t{}\t{}\t{}\n".format(snp, chr, pos, posterior_beta[idx],
                                                                               posterior_prob_c[idx],
                                                                               snp_zs[idx], abs(posterior_beta[idx]))
                            posterior_betas_fh.write(line)
                        posterior_betas_fh.close()
                    except ValueError:
                        pg.bail("Failed to open out file")
                # saving genes results =================================================================================
                if self.options.genes_out_file is not None:
                    try:
                        posterior_genes_fh = open(self.options.genes_out_file + inf_suffix, 'w')
                        line = "gene\tposterior_gene_prob\twindow\tdecay\tgene_freq\n"
                        posterior_genes_fh.write(line)
                        for gene in self.genes_in_region:
                            l = self.genes_extension_index[gene]
                            line = "{}\t{}\t{}\t{}\t{}\n".format(gene, posterior_gene_prob[l],
                                                                 posterior_window, posterior_decay,
                                                                 posterior_gene_freq[l])
                            posterior_genes_fh.write(line)
                        posterior_genes_fh.close()
                    except ValueError:
                        pg.bail("Failed to open out file")
                # saving extended genes results ========================================================================
                if self.options.genes_out_file is not None:
                    try:
                        posterior_genes_fh = open(self.options.genes_out_file + "ext" + inf_suffix, 'w')
                        line = "gene\tposterior_gene_prob\twindow\tdecay\tgene_freq\n"
                        posterior_genes_fh.write(line)
                        for gene in self.genes_extension:
                            l = self.genes_extension_index[gene]
                            line = "{}\t{}\t{}\t{}\t{}\n".format(gene, posterior_gene_prob[l],
                                                                 posterior_window, posterior_decay,
                                                                 posterior_gene_freq[l])
                            posterior_genes_fh.write(line)
                        posterior_genes_fh.close()
                    except ValueError:
                        pg.bail("Failed to open out file")
                # reset of variables for new block of gibbs iterations
                k = -1  # reset of k so the process can continue
                pos_expected_beta_sum = np.zeros((num_chains, num_snps))
                pos_prob_pi_sum = np.zeros((num_chains, num_genes))
                pos_prob_c_sum = np.zeros((num_chains, num_snps))
                window_mu_sum = np.zeros(num_chains)
                window_s_sum = np.zeros(num_chains)
                if not infinite:
                    break
        # === end of gibbs =============================================================================================
        if plot_window_regression:
            plt.waitforbuttonpress()
        return posterior_beta, exit_status

    def plot_results(self):
        pg.plot_data(data=self.stored_data)

    # =============== HANDLING FUNCTIONS ==========================

    def rename_files(self):
        chrom = str(list(self.user_range.keys())[0])
        # start = str(self.user_range[chrom][0])
        # end = str(self.user_range[chrom][1])
        # unique_id = chrom + "_" + start + "_" + end
        if self.options.genes_out_file is None:
            genes_out_file = self.options.out_folder + ".genes"  # + unique_id + ".genes"
            self.options.genes_out_file = genes_out_file
        if self.options.betas_out_file is None:
            betas_out_file = self.options.out_folder + ".betas"  # + unique_id + ".betas"
            self.options.betas_out_file = betas_out_file
        if self.options.data_out_file is None:
            data_out_file = self.options.out_folder + ".data"  # + unique_id + ".data"
            self.options.data_out_file = data_out_file
        if self.options.ld_folder is not None:
            self.options.ld_file = self.options.ld_folder + chrom + ".gz"
        if self.options.dentist_folder is not None:
            self.options.dentist_file = self.options.dentist_folder + chrom + ".den"

    def delete_tmp_files(self):
        if self.options.tabix is not None:
            os.remove(self.options.sumstats_file)
            os.remove(self.options.ld_file)

    def create_tabix_files(self):
        if self.options.tabix is not None:
            chrom = str(list(self.extended_range.keys())[0])
            start = str(self.extended_range[chrom][0])
            end = str(self.extended_range[chrom][1])
            my_range = chrom + ":" + start + "-" + end
            tabix = self.options.tabix
            sumstats_tabix = self.options.tmp_folder + "pegs_gwas_" + chrom + "_" + start + "_" + end
            ld_tabix = self.options.tmp_folder + "pegs_ld_" + chrom + "_" + start + "_" + end
            sumstats_command = tabix + " -h " + self.options.sumstats_file + " " + my_range + \
                               " > " + sumstats_tabix + ";"
            ld_command = tabix + " -h " + self.options.ld_file + " " + my_range + " > " + ld_tabix + ";"
            # Remove files if already exist
            if os.path.exists(sumstats_tabix):
                os.remove(sumstats_tabix)
            if os.path.exists(ld_tabix):
                os.remove(ld_tabix)
            # Create new files, read at the end of the command force popen to wait until it finish
            os.popen(sumstats_command).read()
            os.popen(ld_command).read()
            self.options.sumstats_file = sumstats_tabix
            self.options.ld_file = ld_tabix
            self.tabix_files_created = True

    """
    def plot_genes_pos(self, window):
        open_data = []
        close_data = []
        high_data = []
        low_data = []
        gene_names = []
        sort_genes = dict(sorted(self.genes_position.items(), key=lambda item: item[1]))
        for gene in sort_genes:
            if gene in self.genes_extension_index:
                l = self.genes_extension_index[gene]
                gene_names.append(gene)
                f_gene_chrom, f_gene_start, f_gene_end = self.genes_position[gene]
                open_data.append(f_gene_start)
                close_data.append(f_gene_end)
                high_data.append(f_gene_end + window[l])
                low_data.append(f_gene_start - window[l])
        fig = go.Figure(data=[go.Candlestick(x=gene_names,
                                             open=open_data, high=high_data,
                                             low=low_data, close=close_data)])

        fig.show()
    """

    def log(self, message, level=TRACE, end_char='\n'):
        log_fh = sys.stdout
        if level <= self.options.debug_level:
            log_fh.write("%s%s" % (message, end_char))
            log_fh.flush()

    # ============================================================

    def filter_snps(self, snp_to_p_dentist):
        filtered_snps = []
        for snp in self.snp_to_z:
            if snp in snp_to_p_dentist:
                filtered_snps.append(snp)
        return filtered_snps

    def set_range_restrictions(self):
        restrict_chrom = self.options.chrom
        restrict_start = self.options.start
        restrict_end = self.options.end
        if self.options.range:
            if restrict_chrom is None or restrict_start is None or restrict_end is None:
                restrict_range = self.options.range.split(':')
                if len(restrict_range) != 2:
                    pg.bail("Error: Couldn't parse range %s (must be of form chr:start-end)" % restrict_range)
                restrict_chrom = pg.clean_chrom(restrict_range[0])
                restrict_range = restrict_range[1].split('-')
                if len(restrict_range) != 2:
                    pg.bail("Error: Couldn't parse range %s (must be of form chr:start-end)" % restrict_range)
                restrict_start = int(restrict_range[0])
                restrict_end = int(restrict_range[1])
            else:
                self.log("Ignoring --range because --chrom, --start, and --end have all been specified")
        if self.options.cluster:
            instance_id = self.options.job_id
            instance_start = ((instance_id - 1) * self.options.block_size)
            instance_end = (instance_id * self.options.block_size)
            user_range = {restrict_chrom: (instance_start, instance_end)}
            extended_range = \
                {restrict_chrom: (instance_start - self.options.extension, instance_end + self.options.extension)}
        else:
            user_range = {restrict_chrom: (restrict_start, restrict_end)}
            extended_range = \
                {restrict_chrom: (restrict_start - self.options.extension, restrict_end + self.options.extension)}
        self.user_range = user_range
        self.extended_range = extended_range
        self.log("User region: {}:{}-{}".format(restrict_chrom, user_range[restrict_chrom][0],
                                                user_range[restrict_chrom][1]), self.DEBUG)
        self.log("Extended region: {}:{}-{}".format(restrict_chrom, extended_range[restrict_chrom][0],
                                                    extended_range[restrict_chrom][1]), self.DEBUG)

    def read_gene_file(self, chrom_data_ranges):
        self.log("Reading gene file...", self.INFO)
        with pg.open_gz(self.options.gene_file) as gene_fh:
            gene_header = gene_fh.readline().strip().split()
            gene_id_col = pg.get_col(gene_header, self.options.gene_id_col)
            gene_chrom_col = pg.get_col(gene_header, self.options.gene_chrom_col)
            gene_start_col = pg.get_col(gene_header, self.options.gene_start_col)
            gene_end_col = pg.get_col(gene_header, self.options.gene_end_col)
            for line in gene_fh:
                cols = line.strip().split()
                try:
                    gene = pg.get_value(cols, gene_id_col)
                    chrom = pg.clean_chrom(pg.get_value(cols, gene_chrom_col))
                    start = int(pg.get_value(cols, gene_start_col))
                    end = int(pg.get_value(cols, gene_end_col))
                except IndexError:
                    self.log("Skipping line %s: not enough columns" % line.strip())
                if chrom in chrom_data_ranges:
                    if (chrom_data_ranges[chrom][0] - self.options.extension) < start <= (chrom_data_ranges[chrom][1] +
                                                                                          self.options.extension):
                        self.genes_position[gene] = (chrom, start, end)
                        self.genes_extension[gene] = 0
                    if chrom_data_ranges[chrom][0] < start <= chrom_data_ranges[chrom][1]:
                        self.genes_in_region[gene] = 0
        self.log("Genes in the user region: {}. "
                 "Genes in the extended region: {}".format(len(self.genes_in_region),
                                                           len(self.genes_extension)), self.INFO)
        if len(self.genes_in_region) == 0:
            self.bad_region = True
        gene_fh.close()

    def read_sumstats(self):
        self.log("Reading sumstats...", self.DEBUG)
        with pg.open_gz(self.options.sumstats_file) as sumstats_fh:
            sumstats_header = sumstats_fh.readline().strip().split()
            sumstats_id_col = pg.get_col(sumstats_header, self.options.sumstats_id_col)
            sumstats_beta_col = pg.get_col(sumstats_header, self.options.sumstats_beta_col, require_match=False)
            sumstats_se_col = pg.get_col(sumstats_header, self.options.sumstats_se_col, require_match=False)
            sumstats_n_col = pg.get_col(sumstats_header, self.options.sumstats_n_col, require_match=False)
            sumstats_effect_allele_col = pg.get_col(sumstats_header, self.options.sumstats_effect_allele_col,
                                                    require_match=False)
            sumstats_other_allele_col = \
                pg.get_col(sumstats_header, self.options.sumstats_other_allele_col, require_match=False)
            if sumstats_beta_col is None or sumstats_se_col is None:
                pg.bail("Require --sumstats-beta-col and --sumstats-se-col")
            for line in sumstats_fh:
                cols = line.strip().split()
                z = None
                n = None
                try:
                    snp = pg.get_value(cols, sumstats_id_col)
                    n = float(pg.get_value(cols, sumstats_n_col))
                    beta = float(pg.get_value(cols, sumstats_beta_col))
                    se = float(pg.get_value(cols, sumstats_se_col))
                    z = beta / se
                    p = sc.stats.norm.sf(abs(z)) * 2
                    if sumstats_effect_allele_col is not None and sumstats_other_allele_col is not None:
                        self.snp_to_effect_allele[snp] = pg.get_value(cols, sumstats_effect_allele_col).upper()
                        self.snp_to_other_allele[snp] = pg.get_value(cols, sumstats_other_allele_col).upper()
                except IndexError:
                    self.log("Skipping line %s: not enough columns" % line.strip())

                self.snp_to_beta[snp] = beta
                self.snp_to_se[snp] = se
                self.snp_to_z[snp] = z
                self.snp_to_n[snp] = n
                self.snp_to_p[snp] = p
        self.log("SNPs in sumstats: %s" % len(self.snp_to_z), self.DEBUG)
        sumstats_fh.close()

    def allele_alignment(self):
        if len(self.snp_to_effect_allele) > 0 and self.options.ld_allele_file:
            self.log("Allele correction 1 using allele file...", self.INFO)
            with pg.open_gz(self.options.ld_allele_file) as ld_allele_fh:
                ld_allele_header = ld_allele_fh.readline().strip().split()
                ld_allele_id_col = pg.get_col(ld_allele_header, self.options.ld_allele_id_col)
                ld_allele_effect_allele_col = pg.get_col(ld_allele_header, self.options.ld_allele_effect_allele_col)
                ld_allele_other_allele_col = pg.get_col(ld_allele_header, self.options.ld_allele_other_allele_col)
                for line in ld_allele_fh:
                    cols = line.strip().split()
                    snp = pg.get_value(cols, ld_allele_id_col)
                    effect_allele = pg.get_value(cols, ld_allele_effect_allele_col).upper()
                    other_allele = pg.get_value(cols, ld_allele_other_allele_col).upper()
                    if snp in self.snp_to_z and snp in self.snp_to_effect_allele and snp in self.snp_to_other_allele:
                        if self.snp_to_effect_allele[snp] == effect_allele and \
                                self.snp_to_other_allele[snp] == other_allele:
                            pass
                        elif self.snp_to_effect_allele[snp] == other_allele and \
                                self.snp_to_other_allele[snp] == effect_allele:
                            self.snp_to_z[snp] = -self.snp_to_z[snp]
                        else:
                            other_allele_c = pg.complement_allele(other_allele)
                            effect_allele_c = pg.complement_allele(effect_allele)
                            if other_allele_c is not None and effect_allele_c is not None:
                                if self.snp_to_effect_allele[snp] == effect_allele_c \
                                        and self.snp_to_other_allele[snp] == other_allele_c:
                                    pass
                                elif self.snp_to_effect_allele[snp] == other_allele_c \
                                        and self.snp_to_other_allele[snp] == effect_allele_c:
                                    self.snp_to_z[snp] = -self.snp_to_z[snp]
                                else:
                                    self.log("Could not match alleles in LD file (%s, %s) "
                                             "to alleles in sumstats file (%s, %s) for SNP %s" % (
                                                 effect_allele, other_allele,
                                                 self.snp_to_effect_allele[snp],
                                                 self.snp_to_other_allele[snp],
                                                 snp))
                                    del self.snp_to_z[snp]
                                    del self.snp_to_n[snp]
                            else:
                                self.log("Could not convert alleles in LD file (%s, %s) for SNP %s"
                                         % (effect_allele, other_allele, snp))
                                del self.snp_to_z[snp]
                                del self.snp_to_n[snp]
            ld_allele_fh.close()
        else:
            self.log("Allele correction skipped, allele file was not provided.", self.INFO)

    def read_dentist(self):
        snp_to_p_dentist = {}
        self.log("Reading dentist, using threshold: {}...".format(self.options.dentist_thr), self.DEBUG)
        with pg.open_gz(self.options.dentist_file) as dentist_fh:
            sumstats_id_col = 0
            sumstats_p_col = 2
            for line in dentist_fh:
                cols = line.strip().split()
                try:
                    snp = pg.get_value(cols, sumstats_id_col)
                    p_dentist = float(pg.get_value(cols, sumstats_p_col))
                except IndexError:
                    self.log("Skipping line %s: not enough columns" % line.strip())
                if p_dentist <= self.options.dentist_thr:
                    snp_to_p_dentist[snp] = p_dentist
        dentist_fh.close()
        return snp_to_p_dentist

    def dentist_filter(self):
        if self.options.dentist_file is not None:
            snp_to_p_dentist = self.read_dentist()  # SNP with values below the parameter are eliminated
            self.snps_in_region = self.filter_snps(snp_to_p_dentist)
            self.log("SNPs with z-score and dentist p-value %s" % len(self.snps_in_region), self.DEBUG)
        else:
            self.snps_in_region = list(self.snp_to_z.keys())
            self.log("Dentist file not found, SNPs with z-score %s" % len(self.snps_in_region), self.DEBUG)

    def read_ld(self, extended_range):
        have_negative = False
        self.log("Reading LD, reading and setting SNP correlation values...", self.DEBUG)
        # in this block we get the chrom and pos of SNPs
        with pg.open_gz(self.options.ld_file) as ld_fh:
            ld_header = ld_fh.readline().strip().split()
            ld_id1_col = pg.get_col(ld_header, self.options.ld_id1_col)
            ld_chrom1_col = pg.get_col(ld_header, self.options.ld_chrom1_col)
            ld_pos1_col = pg.get_col(ld_header, self.options.ld_pos1_col)
            ld_id2_col = pg.get_col(ld_header, self.options.ld_id2_col)
            ld_chrom2_col = pg.get_col(ld_header, self.options.ld_chrom2_col)
            ld_pos2_col = pg.get_col(ld_header, self.options.ld_pos2_col)
            ld_r_col = pg.get_col(ld_header, self.options.ld_r_col)
            for line in ld_fh:
                try:
                    cols = line.strip().split()
                    value = float(pg.get_value(cols, ld_r_col))
                    if value < 0:
                        have_negative = True
                    if abs(value) < self.options.min_ld_threshold:
                        continue
                    snp_1 = pg.get_value(cols, ld_id1_col)
                    snp_2 = pg.get_value(cols, ld_id2_col)

                    if snp_1 not in self.snps_in_region or snp_2 not in self.snps_in_region:
                        continue
                    snp_1_chr = pg.clean_chrom(pg.get_value(cols, ld_chrom1_col))
                    snp_1_pos = int(pg.get_value(cols, ld_pos1_col))
                    snp_2_chr = pg.clean_chrom(pg.get_value(cols, ld_chrom2_col))
                    snp_2_pos = int(pg.get_value(cols, ld_pos2_col))

                    if snp_1_chr not in extended_range or snp_2_chr not in extended_range:
                        continue
                    if snp_1_pos < extended_range[snp_1_chr][0] or snp_2_pos < extended_range[snp_1_chr][0]:
                        continue
                    if snp_1_pos > extended_range[snp_1_chr][1] or snp_2_pos > extended_range[snp_1_chr][1]:
                        continue

                    self.snp_to_chr_pos[snp_1] = (snp_1_chr, snp_1_pos)
                    self.snp_to_chr_pos[snp_2] = (snp_2_chr, snp_2_pos)
                    self.ld_dict[snp_1, snp_2] = value
                except IndexError:
                    self.log("Skipping line %s: not enough columns" % line.strip())
                    continue
        if not have_negative and not self.options.all_positive:
            pg.bail("All r values are positive; did you pass in r2 values? If not, specify --all-positive to proceed")
        # Erasing snps in the region without correlation information
        self.snps_in_region = [snp for snp in self.snps_in_region if snp in self.snp_to_chr_pos]
        self.log("SNPs with ld, z-score, and dentist in the selected extended region: %s" % len(self.snps_in_region),
                 self.DEBUG)
        ld_fh.close()

    def create_components(self):
        component = 0
        self.component_to_snp[component] = set()
        max_component_size = 2
        self.log("Creation of ld blocks...", self.DEBUG)
        for entry in dict(sorted(self.ld_dict.items(), key=lambda x: x[1], reverse=True)):
            snp_1 = entry[0]
            snp_2 = entry[1]
            if snp_1 not in self.snp_to_component and snp_2 not in self.snp_to_component:
                # component += 1
                self.snp_to_component[snp_1] = component
                self.snp_to_component[snp_2] = component
                # self.component_to_snp[component] = set()
                self.component_to_snp[component].add(snp_1)
                self.component_to_snp[component].add(snp_2)
            elif snp_1 in self.snp_to_component and snp_2 not in self.snp_to_component:
                if len(self.component_to_snp[self.snp_to_component[snp_1]]) < self.options.max_component_size:
                    self.snp_to_component[snp_2] = self.snp_to_component[snp_1]
                    self.component_to_snp[self.snp_to_component[snp_1]].add(snp_2)
                else:
                    component += 1
                    self.snp_to_component[snp_2] = component
                    self.component_to_snp[component] = set()
                    self.component_to_snp[component].add(snp_2)
            elif snp_2 in self.snp_to_component and snp_1 not in self.snp_to_component:
                if len(self.component_to_snp[self.snp_to_component[snp_2]]) < self.options.max_component_size:
                    self.snp_to_component[snp_1] = self.snp_to_component[snp_2]
                    self.component_to_snp[self.snp_to_component[snp_2]].add(snp_1)
                else:
                    component += 1
                    self.snp_to_component[snp_1] = component
                    self.component_to_snp[component] = set()
                    self.component_to_snp[component].add(snp_1)
            elif snp_2 in self.snp_to_component and snp_1 in self.snp_to_component \
                    and self.snp_to_component[snp_1] != self.snp_to_component[snp_2]:
                if len(self.component_to_snp[self.snp_to_component[snp_1]]) + \
                        len(self.component_to_snp[self.snp_to_component[snp_2]]) <= \
                        self.options.max_component_size:
                    component_1 = self.snp_to_component[snp_1]
                    component_2 = self.snp_to_component[snp_2]
                    for snp in self.component_to_snp[component_2]:
                        self.snp_to_component[snp] = component_1
                    self.component_to_snp[component_1] = \
                        self.component_to_snp[component_1].union(self.component_to_snp[component_2])
                    self.component_to_snp.pop(component_2)

            if len(self.component_to_snp[self.snp_to_component[snp_1]]) > max_component_size:
                max_component_size = len(self.component_to_snp[self.snp_to_component[snp_1]])
        self.log("size of biggest component: {}".format(max_component_size), self.DEBUG)

    def get_data(self):
        # Get chromosomes and position of the genes to be analyzed
        self.read_gene_file(self.user_range)
        # Setting values from the GWAS study all SNPs
        self.read_sumstats()
        # Read Dentist file and intersect with sumstats
        self.dentist_filter()
        # Sample size filter
        if len(self.snps_in_region) != 0:
            max_n = max(self.snp_to_n.values())
            to_delete = [snp for snp in self.snps_in_region if self.snp_to_n[snp] < self.options.n_threshold * max_n]
            before_snps = len(self.snps_in_region)
            for snp in to_delete:
                self.snps_in_region.remove(snp)
            self.log("Deleted %d SNPs due to sample size filtering at %.3g of max; %d remaining" % (
                before_snps - len(self.snps_in_region), self.options.n_threshold, len(self.snps_in_region)), self.DEBUG)
        # check if there are SNPs in region ====
        if len(self.snps_in_region) == 0:
            self.bad_region = True
            return
        # ======================================
        # fill all the dictionaries with the SNPs inside the extended region from the LD file
        self.read_ld(self.extended_range)
        # check if there are SNPs in region ====
        if len(self.snps_in_region) == 0:
            self.bad_region = True
            return
        # ======================================
        # Allele correction using the alleles defined in different file
        self.allele_alignment()
        # Using LD information to create components
        self.create_components()
        # Create a index for each snp in the region
        self.create_region_index()

    def add_genes_to_component(self, genes, snp_to_index, snp_to_chr_pos):
        genes_in_component = set()
        for gene in genes:
            causal_chrom, causal_start, causal_end = genes[gene]
            causal_start -= self.options.causal_window
            causal_end += self.options.causal_window
            for snp in snp_to_index:
                snp_chrom, snp_pos = snp_to_chr_pos[snp]
                if snp_chrom == causal_chrom and causal_start < snp_pos < causal_end:
                    genes_in_component.add(gene)
                    break
        return genes_in_component

    def get_single_snp_link(self, f_gene, f_snp):
        f_snp_chrom, f_snp_pos = self.snp_to_chr_pos[f_snp]
        f_causal_chrom, f_causal_start, f_causal_end = self.genes_position[f_gene]
        f_causal_start -= self.options.causal_window
        f_causal_end += self.options.causal_window
        if f_snp_chrom == f_causal_chrom and f_causal_end > f_snp_pos > f_causal_start:
            f_snp_link = True
        else:
            f_snp_link = False
        return f_snp_link

    def get_gene_to_snp_link(self):
        gene_to_snp_link = np.zeros((len(self.genes_extension), len(self.snps_in_region)))
        for gene in self.genes_extension:
            self.log("SNPs in gene %s: " % gene, self.TRACE)
            for snp in self.snps_in_region:
                snp_link = self.get_single_snp_link(gene, snp)
                gene_to_snp_link[self.genes_extension_index[gene]][self.snps_in_region_index[snp]] = snp_link
        return gene_to_snp_link

    def matrix_stabilization(self, matrix):
        # stabilization achieved by multiply the diagonal by a small factor to avoid singular matrix
        matrix_size = matrix.shape[0]  # this is a square matrix
        null_cov_diag = np.diag(matrix)
        stab_matrix = np.diag(np.ones(matrix_size)) * (np.mean(np.diag(matrix)) * self.options.stab_term_frac)
        null_cov_matrix = matrix + stab_matrix

        # eigen value decomposition to assure the resulting correlation matrix is positive semi-definite
        try:
            l, Q = np.linalg.eigh(null_cov_matrix)
        except np.linalg.linalg.LinAlgError:
            self.log("Exception raised during eigen decomposition")
            raise

        l[l < 0] = 0
        null_cov_matrix = np.dot(Q, np.dot(np.diag(l), Q.transpose())) + stab_matrix

        # transformation to restore the values changed in the stabilization step
        null_cov_norm_matrix = np.diag(np.sqrt(null_cov_diag / np.diag(null_cov_matrix)))
        null_cov_matrix = np.dot(null_cov_norm_matrix, np.dot(null_cov_matrix, null_cov_norm_matrix))
        return null_cov_matrix

    def get_component_data(self, component, k, sampled_beta):
        if k == 0:
            # Assigment of the SNPs to use in the analysis
            snps_in_component = self.component_to_snp[component]
            # Create an index for each snp and for each component
            snp_in_comp_index = pg.create_component_index(snps_in_component)
            # set the z-scores and sample size n of the SNP in the component
            snp_rsid, snp_betas, snp_ses, snp_zs, snp_ns = self.extract_component_snps(snps_in_component,
                                                                                       snp_in_comp_index)
            # set the correlation values of the SNPs in the component
            correlation = self.get_component_correlation(snps_in_component, snp_in_comp_index)
            # get the correlation matrix and de corr matrix of the SNP with z-score i.e. the observable
            correlation = self.matrix_stabilization(correlation)
            # correction of the allele using the correlation matrix
            snp_zs, snp_betas = pg.allele_correction_2(snp_zs, snp_betas, snp_ns, snp_rsid, correlation)

            self.component_rsid[component] = snp_rsid
            self.component_snp_betas[component] = snp_betas
            self.component_correlation[component] = correlation
            self.component_snp_in_comp_index[component] = snp_in_comp_index
            comp_sample_beta = np.copy(snp_betas)
        else:
            snp_rsid = self.component_rsid[component]
            snp_betas = self.component_snp_betas[component]
            correlation = self.component_correlation[component]
            snp_in_comp_index = self.component_snp_in_comp_index[component]
            comp_sample_beta = np.copy(snp_betas)
            for snp in snp_rsid:
                comp_sample_beta[snp_in_comp_index[snp]] = sampled_beta[self.snps_in_region_index[snp]]
        return snp_rsid, snp_betas, comp_sample_beta, correlation

    def extract_region_snps(self):
        num = len(self.snps_in_region)
        snp_betas = np.zeros(num)
        snp_ses = np.zeros(num)
        snp_zs = np.zeros(num)
        snp_ns = np.zeros(num)
        snp_rsid = np.array([None] * num)
        for snp in self.snps_in_region:  # SNPs with LD
            if snp in self.snp_to_z:  # SNPs with z-score
                snp_betas[self.snps_in_region_index[snp]] = self.snp_to_beta[snp]
                snp_ses[self.snps_in_region_index[snp]] = self.snp_to_se[snp]
                snp_zs[self.snps_in_region_index[snp]] = self.snp_to_z[snp]
                snp_ns[self.snps_in_region_index[snp]] = self.snp_to_n[snp]
                snp_rsid[self.snps_in_region_index[snp]] = snp
        return snp_rsid, snp_betas, snp_ses, snp_zs, snp_ns

    def get_region_correlation(self):
        cor = np.identity(len(self.snps_in_region), dtype=np.float64)
        index_list = list(self.snps_in_region)
        for i in range(len(index_list)):
            for j in range(i):
                snp_1 = index_list[i]
                snp_2 = index_list[j]
                if (snp_1, snp_2) in self.ld_dict:
                    cor[self.snps_in_region_index[snp_1], self.snps_in_region_index[snp_2]] = self.ld_dict[snp_1, snp_2]
                    cor[self.snps_in_region_index[snp_2], self.snps_in_region_index[snp_1]] = self.ld_dict[snp_1, snp_2]
                elif (snp_2, snp_1) in self.ld_dict:
                    cor[self.snps_in_region_index[snp_1], self.snps_in_region_index[snp_2]] = self.ld_dict[snp_2, snp_1]
                    cor[self.snps_in_region_index[snp_2], self.snps_in_region_index[snp_1]] = self.ld_dict[snp_2, snp_1]
        return cor

    def create_region_index(self):
        index = 0
        for snp in sorted(self.snps_in_region):
            self.snps_in_region_index[snp] = index
            index += 1
        gene_index = 0
        for gene in sorted(self.genes_extension):
            self.genes_extension_index[gene] = gene_index
            gene_index += 1

    def extract_component_snps(self, snps_in_component, snp_to_index):
        num = len(snps_in_component)
        snp_betas = np.zeros(num)
        snp_ses = np.zeros(num)
        snp_zs = np.zeros(num)
        snp_ns = np.zeros(num)
        snp_rsid = np.array([None] * num)
        for snp in snps_in_component:
            snp_betas[snp_to_index[snp]] = self.snp_to_beta[snp]
            snp_ses[snp_to_index[snp]] = self.snp_to_se[snp]
            snp_zs[snp_to_index[snp]] = self.snp_to_z[snp]
            snp_ns[snp_to_index[snp]] = self.snp_to_n[snp]
            snp_rsid[snp_to_index[snp]] = snp
        return snp_rsid, snp_betas, snp_ses, snp_zs, snp_ns

    def get_component_correlation(self, snps_in_component, snp_to_index):
        cor = np.identity(len(snps_in_component), dtype=np.float64)
        index_list = list(snps_in_component)
        for i in range(len(index_list)):
            for j in range(i):
                snp_1 = index_list[i]
                snp_2 = index_list[j]
                if (snp_1, snp_2) in self.ld_dict:
                    cor[snp_to_index[snp_1], snp_to_index[snp_2]] = self.ld_dict[snp_1, snp_2]
                    cor[snp_to_index[snp_2], snp_to_index[snp_1]] = self.ld_dict[snp_1, snp_2]
                elif (snp_2, snp_1) in self.ld_dict:
                    cor[snp_to_index[snp_1], snp_to_index[snp_2]] = self.ld_dict[snp_2, snp_1]
                    cor[snp_to_index[snp_2], snp_to_index[snp_1]] = self.ld_dict[snp_2, snp_1]
        return cor

    def get_snp_to_gene_dist(self):
        gene_to_snp_distance = np.zeros((len(self.genes_extension), len(self.snps_in_region)))
        for gene in self.genes_extension:
            self.log("SNPs in gene %s: " % gene, self.TRACE)
            for snp in self.snps_in_region:
                snp_dist = self.get_snp_dist(gene, snp)
                gene_to_snp_distance[self.genes_extension_index[gene]][self.snps_in_region_index[snp]] = snp_dist
        return gene_to_snp_distance

    def get_snp_dist(self, f_gene, f_snp):
        f_snp_chrom, f_snp_pos = self.snp_to_chr_pos[f_snp]
        f_gene_chrom, f_gene_start, f_gene_end = self.genes_position[f_gene]
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

    def get_component_data_for_multi_chain(self, component, k, sampled_beta, local_test):
        # Assigment of the SNPs to use in the analysis
        num_chains = sampled_beta.shape[0]
        snps_in_component = self.component_to_snp[component]
        num_snps_in_component = len(snps_in_component)
        if k == 0:
            # Create an index for each snp and for each component
            snp_in_comp_index = pg.create_component_index(snps_in_component)
            # set the z-scores and sample size n of the SNP in the component
            snp_rsid, snp_betas, snp_ses, snp_zs, snp_ns = self.extract_component_snps(snps_in_component,
                                                                                       snp_in_comp_index)
            # set the correlation values of the SNPs in the component
            correlation = self.get_component_correlation(snps_in_component, snp_in_comp_index)
            # get the correlation matrix and de corr matrix of the SNP with z-score i.e. the observable
            correlation = self.matrix_stabilization(correlation)
            # correction of the allele using the correlation matrix
            if not local_test:
                self.log("Allele correction 2 correlation based, components...", self.DEBUG)
                snp_zs, snp_betas = pg.allele_correction_2(snp_zs, snp_betas, snp_ns, correlation,
                                                           self.options.detect_flips)
            # normalization of betas
            snp_betas = (snp_betas / snp_ses) / np.sqrt(self.options.sample_size)

            self.component_rsid[component] = snp_rsid
            self.component_correlation[component] = correlation
            self.component_snp_in_comp_index[component] = snp_in_comp_index
            comp_sample_beta = np.zeros((num_chains, num_snps_in_component))
            for snp in snp_rsid:
                comp_sample_beta[:, snp_in_comp_index[snp]] = sampled_beta[:, self.snps_in_region_index[snp]]
            component_norm_betas = np.zeros((num_chains, num_snps_in_component))
            for chain in range(num_chains):
                component_norm_betas[chain] = np.copy(snp_betas)
            self.component_snp_betas[component] = component_norm_betas
        else:
            snp_rsid = self.component_rsid[component]
            component_norm_betas = self.component_snp_betas[component]
            correlation = self.component_correlation[component]
            snp_in_comp_index = self.component_snp_in_comp_index[component]
            comp_sample_beta = np.zeros((num_chains, num_snps_in_component))
            for snp in snp_rsid:
                comp_sample_beta[:, snp_in_comp_index[snp]] = sampled_beta[:, self.snps_in_region_index[snp]]
        return snp_rsid, component_norm_betas, comp_sample_beta, correlation


class GibbsData:
    print_beta = None
    inf_beta = None
    posterior_beta = None
    snp_to_chr_pos = None
    num_snps = None
    num_genes = None
    snps_in_region_index = None
    genes_extension_index = None
    pi_vec_print = None
    p_vec_print = None
    options = None
    p_in_vec_print = None
    p_out_vec_print = None
    heritability_vec_print = None
    ldpred_heritability_vec = None
    pi_pos_vec = None
    shrinkage_vec_print = None
    posterior_sigma_vec_print = None
    window_mu_vec_print = None
    pi_avg_print = None
    epsilon = None
    shrinkage = None
    pos_prob_pi_trace = None
    within = None
    between = None
    ratio_var = None
    gibbs_iter = None
    burn_in = None
    prior_shrinkage = None
    max_beta_pos_prob = None
    window_s_vec_print = None
    prior_window_decay = None
    trace_sampled_G = None
    marginal_betas = None
    marginal_se = None
