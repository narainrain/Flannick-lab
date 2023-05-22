import math

from models import PEGS
import pegs_utils as pg
import numpy as np
import pickle as pl
import matplotlib
import matplotlib.pyplot as plt
matplotlib.use('Qt5Agg')


test_options = pg.arg_settings()
# ==========================================================================================================
from_file = False
# =========================================
gene = "APOB"  # "NPAP1" or "GPAM" or "APOB"
base_file = "/Users/ellamas/PycharmProjects/PEGS/PEGS/data/sigma_sensitivity_variation_of_h/test.data" + gene + "_"
# ==========================================================================================================
if gene == "NPAP1":
    test_options.range = "15:24000000-25000000"
if gene == "GPAM":
    test_options.range = "10:113000000-114000000"
if gene == "APOB":
    test_options.range = "2:21000000-22000000"

num_samples = 10
# Instance creation
my_model = PEGS(test_options)
gene_index = my_model.genes_extension_index
num_genes = len(gene_index)
gene_posterior_probabilities = np.zeros((num_samples, num_genes))
posterior_betas = np.zeros(num_samples)
# p_vec = np.linspace(0.00001, 0.01, num_samples)
h_vec = np.logspace(-3, 0, 10)

for i in range(num_samples):
    if from_file:
        file = base_file + str(i)
        with open(file, 'rb') as f:
            data = pl.load(f)
        gene_posterior_probabilities[i, :] = data.pi_avg_print[data.gibbs_iter - 1]
        posterior_betas[i] = max(abs(data.posterior_beta))
    else:
        # my_model.options.p = p_vec[i]
        my_model.options.heritability = h_vec[i]
        my_model.gibbs_multi_chain(num_chains=20, burn_in=100, max_iter=1000, local_test=True,
                                   convergence_thr=1.0001, file_name=gene + "_" + str(i))
        gene_posterior_probabilities[i, :] = my_model.stored_data.pi_avg_print[my_model.stored_data.gibbs_iter-1]
        posterior_betas[i] = max(abs(my_model.stored_data.posterior_beta))

plt.figure(0)
plt.title("max posterior beta")
plt.xscale("log")
plt.plot(h_vec, posterior_betas)
plt.savefig("../data/out/"+gene+"_betas.pdf", bbox_inches='tight')
size_of_plot = math.ceil(math.sqrt(num_genes))
plt.figure(1)
plt.tight_layout()
for gene_loop in gene_index:
    l = gene_index[gene_loop]
    ax = plt.subplot(size_of_plot, size_of_plot, l+1)
    ax.set_ylim([0, 1])
    ax.set_xscale("log")
    plt.title("Gene: " + str(gene_loop), fontsize=8)
    plt.plot(h_vec, gene_posterior_probabilities[:, l])

plt.savefig("../data/out/"+gene+"_genes.pdf", bbox_inches='tight')
plt.show()

