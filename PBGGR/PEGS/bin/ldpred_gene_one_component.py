from help_functions import *
import time


def BGAM(options):
    start_time = time.time()
    # Variables, initialization and validation -------------------------------------------------------------------------
    open_files()
    initialize_global_var()
    arg_validation()
    # ------------------------------------------------------------------------------------------------------------------

    # Reading and setting data -----------------------------------------------------------------------------------------
    if options.simulate:
        gene_to_log_bf, component_to_snp, gene_extend_to_chrom_start_stop, snp_to_chr_pos, snp_to_z, snp_to_n, ld_dict \
            = simulate_data()
    elif options.create_data:
        gene_to_log_bf = create_data()
        return gene_to_log_bf
    else:
        gene_to_log_bf, component_to_snp, gene_extend_to_chrom_start_stop, snp_to_chr_pos, snp_to_z, snp_to_n, ld_dict, \
        snp_to_se, snp_to_beta = get_data()
    # ------------------------------------------------------------------------------------------------------------------
    # gene_to_log_bf = {'BEST1': 0, 'FTH1': 0}
    # variables for the Calculation of Bayes Factor --------------------------------------------------------------------
    gene_to_alt_num_snps = dict(gene_to_log_bf)  # used to count the number of SNPs per gene
    gene_to_null_num_snps = dict(gene_to_log_bf)  # used to count the number of SNPs per gene
    gene_to_alt_h = dict(gene_to_log_bf)
    gene_to_null_h = dict(gene_to_log_bf)
    # ------------------------------------------------------------------------------------------------------------------
    # [207] best1
    # [256] APOC1
    for component in [207]:
        # Assigment of the SNPs to use in the analysis  ----------------------------------------------------------------
        snp_to_analyze = component_to_snp[component]
        # --------------------------------------------------------------------------------------------------------------
        # Create an index for each snp and for each component
        snp_to_index = create_index(snp_to_analyze)

        # add genes to the component if they have at least one SNP inside the causal region ----------------------------
        genes_in_component = add_genes_to_component(gene_extend_to_chrom_start_stop, snp_to_index, snp_to_chr_pos)
        if len(genes_in_component) == 0:
            log("Component without genes: %s" % snp_to_analyze, TRACE)
            continue
        # --------------------------------------------------------------------------------------------------------------
        # print(component)
        # print(genes_in_component)
        # continue
        # only calculate component if the genes in gene to log bf are in the component ---------------------------------
        if len(set(gene_to_log_bf.keys()).intersection(genes_in_component)) == 0:
            log("Component without genes for BF: %s" % snp_to_analyze, TRACE)
            continue
        # --------------------------------------------------------------------------------------------------------------

        # set the z-scores and sample size n of the SNP in the component -----------------------------------------------
        snp_zs, snp_ns, component_snps_obs, has_z, snp_betas, snp_ses = extract_snps(snp_to_index, snp_to_z, snp_to_n,
                                                                                     snp_to_beta, snp_to_se)
        print_collected_statistics(snp_to_index, component_snps_obs, snp_to_chr_pos, snp_zs, snp_ns)
        if len(snp_zs) == 0:
            continue
        # --------------------------------------------------------------------------------------------------------------

        # For each SNP in the component calculate the link to gene -----------------------------------------------------
        gene_to_snp_probs = gene_to_snp(genes_in_component, snp_to_index, snp_to_chr_pos,
                                        gene_extend_to_chrom_start_stop, has_z)
        gene_to_snp_distance = gene_to_snp_distance_function(genes_in_component, snp_to_index, snp_to_chr_pos,
                                                             gene_extend_to_chrom_start_stop, has_z)
        # --------------------------------------------------------------------------------------------------------------

        # set the correlation values of the SNPs in the component ------------------------------------------------------
        correlation = set_correlation(snp_to_index, ld_dict)
        # get the correlation matrix and de corr matrix of the SNP with z-score i.e. the observable
        cor_matrix_obs = correlation[has_z, :][:, has_z]
        cor_matrix_obs = matrix_stabilization(cor_matrix_obs)
        # correction of the allele using the correlation matrix
        allele_correction_2(snp_zs, snp_ns, component_snps_obs, cor_matrix_obs)
        # --------------------------------------------------------------------------------------------------------------
        while True:
            omega, status = ldpred_gene_window_prior(snp_betas, options.sample_size, cor_matrix_obs, options.p,
                                       options.heritability, options.total_num_SNPs, snp_to_index, snp_zs,
                                       snp_to_chr_pos, gene_to_snp_probs, gene_to_snp_distance, gene_extend_to_chrom_start_stop)
            if status == "success":
                break
        # end of component loop
    # print_results(gene_to_log_bf, gene_to_null_num_snps, gene_to_alt_num_snps)
    # print_model_setup(DEBUG)
    close_files()
    log("--- %s seconds ---" % (time.time() - start_time), DEBUG)
    return gene_to_log_bf


if __name__ == "__main__":
    options = arg_settings()
    BGAM(options)
