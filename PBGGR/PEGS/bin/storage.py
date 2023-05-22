


def ldpred(beta, n, correlation, p, h_norm, M, snp_to_index, snp_zs, snp_to_chr_pos, burn_in=10, gibbs_iter=500):
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
    :param burn_in number of iterations dump from the mean
    :param gibbs_iter number of samples used for gibbs calculations
    """
    # ======== initialization ===================================
    snps_in_component = beta.shape[0]
    norm_beta = beta  # / se)  # * math.sqrt(n)  used when the marginal betas aren't normalized
    shrinkage = 1 / (1 + ((M * p) / (n * h_norm)))
    # ======== variables for plotting ===================
    print_beta = np.zeros(gibbs_iter)
    max_beta_idx = np.argmax(abs(beta))
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

    for k in range(burn_in + gibbs_iter):
        for j in range(snps_in_component):
            # calculation of beta tilde j (joint beta) ========================
            corr_no_j = np.delete(correlation[j], j, 0)
            sampled_beta_no_j = np.delete(sampled_beta, j, 0)
            joint_beta[j] = norm_beta[j] - np.dot(sampled_beta_no_j.T, corr_no_j)
            # Compute pj according to (6)   ==================================
            p_vec[j] = 1 / (1 + (((1 - p) / p) * math.sqrt(1 + ((n * h_norm) / (M * p))) * math.exp(
                ((-1 / 2) * n * joint_beta[j] ** 2) /
                (1 + ((M * p) / (n * h_norm))))))
            # Sample bj according to (7) ======================================
            aux = np.random.uniform(0, 1, 1)[0]
            mean = joint_beta[j] * shrinkage
            if 0 <= aux <= p_vec[j]:
                variance = (1 / n) * shrinkage
                sampled_beta[j] = np.random.normal(mean, math.sqrt(variance), 1)
            else:
                sampled_beta[j] = 0.0
            omega_vec[j] = p_vec[j] * mean
        # end of cycle of variants inside the component
        if k >= burn_in:
            omega_cap_vec = omega_cap_vec + omega_vec
            # var to print SNP with the highest effect
            print_beta[k - burn_in] = (omega_cap_vec[max_beta_idx] / (k - burn_in + 1))
    # end of gibbs cycle
    # === calculation of posterior betas   =================
    omega_cap_vec = omega_cap_vec / gibbs_iter
    posterior_beta = omega_cap_vec
    # ==== plotting results ================================
    if options.debug_level >= 4:
        plt.figure(1)
        plt.plot(print_beta)
        plt.title("convergence of strongest SNP")
        plt.figure(2)
        plt.scatter(inf_beta, posterior_beta)
        plt.title("infinitesimal vs Gibbs results")
        # sorting snps to plot by genomic position
        sort_snps = dict(sorted(snp_to_chr_pos.items(), key=lambda item: item[1]))
        to_plot = np.zeros(snps_in_component)
        j = 0
        for snp in sort_snps:
            if snp in snp_to_index:
                to_plot[j] = posterior_beta[snp_to_index[snp]]
                j += 1
        plt.figure(3)
        plt.scatter(np.linspace(0, 1, snps_in_component), to_plot)
        plt.title("gibbs results")
        plt.show()
    # saving data ======================================================
    for snp in snp_to_index:
        idx = snp_to_index[snp]
        chr = snp_to_chr_pos[snp][0]
        pos = snp_to_chr_pos[snp][1]
        line = "{}\t{}\t{}\tA\tG\t{}\t{}\t{}\t{}\n".format(snp, chr, pos, posterior_beta[idx], p_vec[idx],
                                                           snp_zs[idx], abs(posterior_beta[idx]))
        file_ldpred.write(line)
    return posterior_beta


def ldpred_gene(beta, n, correlation, p, h_norm, M, snp_to_index, snp_zs, snp_to_chr_pos, gene_to_snp_probs, burn_in=10,
                gibbs_iter=500):
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
    "param gene_to_snp_probs used to plot in genomic order
    :param burn_in number of iterations dump from the mean
    :param gibbs_iter number of samples used for gibbs calculations
    """
    # ======== initialization ===================================
    epsilon = 0.0001
    pi = options.gene_prior
    snps_in_component = beta.shape[0]
    num_genes = len(gene_to_snp_probs)
    norm_beta = beta  # / se)  # * math.sqrt(n)  used when the marginal betas aren't normalized
    shrinkage = 1 / (1 + ((M * p) / (n * h_norm)))
    p_gene = p
    epsilon_gene = epsilon
    # priors definition =================================
    # variables for heritability update
    prior_sigma = h_norm / (M * p)
    prior_alpha = 10000  # confidence on the prior sigma
    prior_beta = prior_sigma * prior_alpha
    posterior_sigma = prior_sigma
    # variables for p update
    prior_strength = 10
    prior_p_alpha = prior_strength
    prior_p_beta = (1/p) * prior_strength
    # ======== variables for plotting ===================
    print_beta = np.zeros(gibbs_iter)
    max_beta_idx = np.argmax(abs(beta))
    pi_vec_print = np.zeros((gibbs_iter, num_genes))
    p_vec_print = np.zeros(gibbs_iter)
    p_in_vec_print = np.zeros(gibbs_iter)
    p_out_vec_print = np.zeros(gibbs_iter)
    heritability_vec_print = np.zeros(gibbs_iter)
    posterior_sigma_vec_print = np.zeros(gibbs_iter)
    ldpred_heritability_vec = np.zeros(gibbs_iter)
    pi_pos_vec = np.zeros(gibbs_iter)
    shrinkage_vec_print = np.zeros(gibbs_iter)
    # ======== printing priors information =====================
    log("prior sigma = {}".format(prior_sigma))
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
    shrinkage_sum = 0
    # G is initialized so all genes are causal
    sampled_G = {}
    for gene in gene_to_snp_probs:
        sampled_G[gene] = 1
    # ======== aux variables ===========================
    gene_list = list(gene_to_snp_probs)
    # ======= gibbs iterations =========================
    for k in range(burn_in + gibbs_iter):
        for j in range(snps_in_component):
            # calculation of beta tilde j (joint beta) ========================
            corr_no_j = np.delete(correlation[j], j, 0)
            sampled_beta_no_j = np.delete(sampled_beta, j, 0)
            joint_beta[j] = norm_beta[j] - np.dot(sampled_beta_no_j.T, corr_no_j)
            # Compute pj and sample cj  ==================================
            if indicator(j, gene_to_snp_probs, sampled_G):
                # p_vec[j] = 1 / (1 + (((1 - p) / p) * math.sqrt(1 + ((n * h_norm) / (M * p))) * math.exp(
                #     ((-1 / 2) * n * joint_beta[j] ** 2) /
                #     (1 + ((M * p) / (n * h_norm))))))
                p_vec[j] = 1 / (1 + (((1 - p) / p) * math.sqrt(1 + (n * posterior_sigma))
                                     * math.exp((-1/2) * ((n * joint_beta[j]**2)/(1 + (1/(n*posterior_sigma)))))))
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
        # end of cycle of variants inside the component
        # compute pi_l and sample gl
        for l in range(num_genes):
            gene = gene_list[l]
            iter_sampled_g = sampled_G.copy()
            iter_sampled_g[gene] = 1  # setting the current gene to one
            multiplication_1 = 1
            for j in range(snps_in_component):
                if indicator(j, gene_to_snp_probs, iter_sampled_g):
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
                if indicator(j, gene_to_snp_probs, iter_sampled_g):
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
        # update of p, heritability, p_gene, epsilon gene
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
                if indicator(j, gene_to_snp_probs, sampled_G):  # inside the genes
                    snps_inside_causal_genes += 1
                    if sampled_beta[j] != 0:
                        M_c_in += 1
                else:  # outside genes
                    if sampled_beta[j] != 0:
                        M_c_out += 1
            snps_outside_causal_genes = snps_in_component - snps_inside_causal_genes
            p_in = np.random.beta(1+M_c_in, 1 + snps_inside_causal_genes - M_c_in)
            p_out = np.random.beta(1+M_c_out, 1 + snps_outside_causal_genes - M_c_out)
            # print("{},{}".format(p_in, p_out))
            p_gene = p_in
            epsilon_gene = p_out
            # epsilon_gene = 0.0001  # p_out
            # print(epsilon_gene)
            # update of heritability
            ldpred_h_norm = np.matmul(np.matmul(sampled_beta.T, correlation), sampled_beta)
            pos_alpha = (prior_alpha + (M/2))
            pos_beta = (prior_beta + (np.matmul(sampled_beta.T, sampled_beta) / 2))
            posterior_sigma = sc.stats.invgamma.rvs(a=pos_alpha, scale=pos_beta, size=1)[0]
            h_norm = posterior_sigma * M * p_in  # heritability is no longer used in gibbs it's used to compare
            # =======================================================================
        # sums for the average calculations
        if k >= burn_in:
            omega_cap_vec = omega_cap_vec + omega_vec
            pi_vec_sum = pi_vec_sum + pi_vec
            # variables to print
            pi_vec_print[k - burn_in] = pi_vec_sum / (k - burn_in + 1)
            print_beta[k - burn_in] = (omega_cap_vec[max_beta_idx] / (k - burn_in + 1))
            p_vec_print[k - burn_in] = p
            heritability_vec_print[k - burn_in] = h_norm
            posterior_sigma_vec_print[k - burn_in] = posterior_sigma
            ldpred_heritability_vec[k - burn_in] = ldpred_h_norm
            p_in_vec_print[k - burn_in] = p_gene
            p_out_vec_print[k - burn_in] = epsilon_gene
            pi_pos_vec[k - burn_in] = pi
            shrinkage_sum += (n * posterior_sigma) / ((n * posterior_sigma) + 1)
            shrinkage_vec_print[k - burn_in] = shrinkage_sum / (k - burn_in + 1)
            shrinkage = shrinkage_sum / (k - burn_in + 1)
            # p = 0.0000000000000000061442

    # end of gibbs cycle
    # === calculation of posterior betas   =================
    omega_cap_vec = omega_cap_vec / gibbs_iter
    posterior_beta = omega_cap_vec
    # ==== plotting results ================================
    if options.debug_level >= 2:
        plt.figure(1)
        plt.subplot(2, 2, 1)
        plt.plot(print_beta, linewidth=1)
        plt.title("convergence of strongest SNP", fontsize=8)
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
        plt.title("gibbs results", fontsize=8)
        # plot of gene probabilities
        ax = plt.subplot(2, 2, 4)
        plt.title("gene posteriors", fontsize=8)
        avg_pi_vec = pi_vec_sum / gibbs_iter
        gene_file = open('../data/ldpred_results/gene.results', 'w')
        NUM_COLORS = num_genes
        cm = plt.get_cmap('gist_rainbow')
        line_type = ["-", "--", "-.", ":"]
        linecycler = cycle(line_type)
        ax.set_prop_cycle('color', [cm(1. * i / NUM_COLORS) for i in range(NUM_COLORS)])
        for l in range(num_genes):
            ax.plot(pi_vec_print[:, l], next(linecycler), label=gene_list[l])
            gene_file.write("{}\t{}\n".format(gene_list[l], avg_pi_vec[l]))
        # plt.legend(loc='upper center', bbox_to_anchor=(0.5, 1.1), ncol=3, fancybox=True, shadow=True)
        gene_file.close()
        plt.tight_layout(pad=0.5)
        # plotting other probabilities
        plt.figure(2)
        plt.subplot(3, 1, 1)
        plt.plot(p_vec_print)
        plt.title("fraction of causal SNPs p in the region, mean {:.4f}".format(p_vec_print.mean()), fontsize=8)
        plt.subplot(3, 1, 2)
        plt.plot(p_in_vec_print)
        plt.title("fraction of causal SNPs p_in given the gene is causal, mean {:.4f}".format(p_in_vec_print.mean()), fontsize=8)
        plt.subplot(3, 1, 3)
        plt.plot(p_out_vec_print)
        plt.title("fraction of causal SNPs p_out given the gene is not causal, mean {:.4f}".format(p_out_vec_print.mean()), fontsize=8)
        plt.tight_layout(pad=0.5)
        plt.figure(3)
        plt.subplot(3, 1, 1)
        plt.plot(heritability_vec_print)
        plt.title("My heritability, mean {:.6f}".format(heritability_vec_print.mean()), fontsize=8)
        plt.subplot(3, 1, 2)
        plt.plot(ldpred_heritability_vec)
        plt.title("LDpred heritability, mean {:.6f}".format(ldpred_heritability_vec.mean()), fontsize=8)
        plt.subplot(3, 1, 3)
        plt.plot(pi_pos_vec)
        plt.title("Posterior gene probability, mean {:.6f}".format(pi_pos_vec.mean()), fontsize=8)
        plt.tight_layout(pad=0.5)
        plt.figure(4)
        plt.subplot(3, 1, 1)
        plt.plot(shrinkage_vec_print)
        plt.title("Mean shrinkage on each iteration", fontsize=8)
        plt.subplot(3, 1, 2)
        plt.plot(posterior_sigma_vec_print)
        plt.title("Posterior sigma, mean {:.6f}".format(posterior_sigma_vec_print.mean()), fontsize=8)
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
    return posterior_beta


def ldpred_gene_window(beta, n, correlation, p, h_norm, M, snp_to_index, snp_zs, snp_to_chr_pos, gene_to_snp_probs,
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
    # priors definition =================================
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
            for l in range(num_genes):
                gene = gene_list[l]
                regression_X = np.concatenate((gene_to_snp_distance[gene], np.zeros(snps_in_component))).reshape(-1, 1)
                regression_y = np.concatenate((I[gene], np.ones(snps_in_component)))
                # regression_X = np.concatenate((gene_to_snp_distance[gene], np.negative(gene_to_snp_distance[gene]))).reshape(-1, 1)
                # regression_y = np.concatenate((I[gene], np.ones(snps_in_component)))
                # regression_X = gene_to_snp_distance[gene].reshape(-1, 1)
                # regression_y = I[gene]
                with warnings.catch_warnings():
                    warnings.filterwarnings('error')
                    try:
                        clf = LogisticRegression(random_state=0).fit(regression_X, regression_y)
                        mean = - clf.intercept_ / clf.coef_
                        s = 1 / clf.coef_
                        window_mu[l] = mean
                        window_s[l] = s
                    except Warning:
                        log("Window not updated for gene {}, failed to converge, using previous value".format(gene))
                        mean = window_mu[l]
                        s = window_s[l]
                if gene == "FEN1" and False:  # used to check the evolution of the window for a gene
                    # debugging plot
                    # log("mean: {}, s: {}".format(mean, s))
                    num_points = 1000
                    x_axis = np.linspace(min(regression_X), max(regression_X), num_points)
                    logit = 1 / (1 + np.exp(-(x_axis - mean) / s)).reshape(num_points, )
                    plt.figure(1)
                    plt.scatter(regression_X, regression_y)
                    plt.plot(x_axis, logit)
                    plt.title("Gene: {}, iter {}".format(gene, k))
                    plt.show()

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
                    if link_denominator == 0:
                        posterior_p_I_l_j_1 = 0
                    else:
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
        limit = 1e-20  # 1e-07
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
    return posterior_beta, exit_status


