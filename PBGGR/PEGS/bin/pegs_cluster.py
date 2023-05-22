from models import PEGS
import pegs_utils as pg
import time

# Instance creation
start_time = time.time()
cluster_options = pg.arg_settings()
my_model = PEGS(cluster_options)
if not my_model.bad_region:
    my_model.gibbs_multi_chain(num_chains=cluster_options.num_chains,
                               burn_in=cluster_options.burn_in,
                               max_iter=cluster_options.max_iter,
                               convergence_thr=cluster_options.convergence_thr)
else:
    if my_model.options.genes_out_file is not None:
        print("Bad region check range")
        try:
            posterior_genes_fh = open(my_model.options.genes_out_file, 'w')
            line = "gene\tposterior_gene_prob\twindow\tdecay\n"
            posterior_genes_fh.write(line)
            posterior_genes_fh.close()
        except ValueError:
            pg.bail("Failed to open out file")
print("--- %s seconds ---" % (time.time() - start_time))
