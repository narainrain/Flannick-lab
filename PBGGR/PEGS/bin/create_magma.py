from gene_bf import *

if __name__ == "__main__":
    options = arg_settings()
    priors = [0.05]
    genes = ["greedy"]
    snps = ["Full_SNP"]  # 'greedy', 'Full_SNP', 'All_one'
    aggregations = [0]  # 1 is additive, 0 is max
    windows = [100]
    for prior in priors:
        for gene in genes:
            for snp in snps:
                for aggregation in aggregations:
                    for window in windows:
                        # change options
                        options.gene_prior = prior
                        options.model = gene
                        options.SNP_model = snp
                        options.additive = aggregation
                        options.causal_window = window
                        print(prior, gene, snp, aggregation, window)
                        BGAM(options)
