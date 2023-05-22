from optparse import OptionParser
import numpy as np
import pickle as pl
import matplotlib.pyplot as plt

# arguments
parser = OptionParser("usage: %prog [options]")
parser.add_option("", "--input-folder", default=None)
parser.add_option("", "--out-file", default=None)
(options, _) = parser.parse_args()

# ======================================================
base_name = options.input_folder + "sensitivity_idx_{}.data"
num_files = 1000
gene_posterior_probabilities = []
posterior_betas = []
flag = True
gene_index = None
num_genes = None
# sorted priors ==========================================================
num_conf = 10
sigma_strength = np.logspace(-2, 3, num_conf)
p_strength = np.logspace(-2, 3, num_conf)
p_in_strength = np.logspace(-2, 3, num_conf)
combination = np.array(np.meshgrid(sigma_strength, p_strength, p_in_strength)).T.reshape(-1, 3)
sort_comb = np.argsort(combination.sum(axis=1))
# =========================================================================
for i in sort_comb:
    file = base_name.format(i)
    try:
        with open(file, 'rb') as f:
            data = pl.load(f)
        gene_posterior_probabilities.append(data.pi_avg_print[data.gibbs_iter - 1])
        posterior_betas.append(max(abs(data.posterior_beta)))
        if flag:
            gene_index = data.genes_extension_index
            num_genes = data.num_genes
            flag = False
    except IOError:
        print("file does not exist")

gene_posterior_probabilities_to_plot = np.asarray(gene_posterior_probabilities)
posterior_betas_to_plot = np.asarray(posterior_betas)


column_values = list(gene_index.keys())

fig, ax = plt.subplots()
im = ax.imshow(gene_posterior_probabilities_to_plot, cmap='Reds')
ax.set_xticks(np.arange(len(column_values)), labels=column_values)
plt.setp(ax.get_xticklabels(), rotation=45, ha="right",
         rotation_mode="anchor")
plt.tick_params(left = False, right = False , labelleft = False ,
                labelbottom = True, bottom = True)
fig.colorbar(im, ax=ax, location='right', anchor=(0, 0.3), shrink=0.7)
ax.set_title("Sensitivity")
fig.tight_layout()
plt.savefig(options.out_file)
