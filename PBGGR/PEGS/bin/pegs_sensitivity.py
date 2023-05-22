from models import PEGS
import pegs_utils as pg
import time
import numpy as np

# sensitivity configurations
num_conf = 10
sigma_strength = np.linspace(1, 10000, num_conf)
p_strength = np.logspace(-2, 3, num_conf)
p_in_strength = np.logspace(-2, 3, num_conf)
combination = np.array(np.meshgrid(sigma_strength, p_strength, p_in_strength)).T.reshape(-1, 3)
# Instance creation
start_time = time.time()
sensitivity_options = pg.arg_settings()
my_model = PEGS(sensitivity_options)
idx = sensitivity_options.job_id - 1
if not my_model.bad_region:
    my_model.gibbs_multi_chain(num_chains=20,
                               burn_in=100,
                               max_iter=3000,
                               convergence_thr=1.000001,
                               prior_alpha_sigma=combination[idx][0],
                               prior_p_alpha=combination[idx][1],
                               prior_p_in_alpha=combination[idx][2]
                               )
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
