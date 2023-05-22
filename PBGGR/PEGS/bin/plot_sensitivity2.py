from optparse import OptionParser
import numpy as np
import pickle as pl
import matplotlib.pyplot as plt
from bokeh.plotting import figure, output_file, save
from bokeh.plotting import ColumnDataSource

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
sigma_real = np.linspace(1, 10000, num_conf)
p_strength = np.logspace(-2, 3, num_conf)
p_in_strength = np.logspace(-2, 3, num_conf)
combination = np.array(np.meshgrid(sigma_strength, p_strength, p_in_strength)).T.reshape(-1, 3)
combination2 = np.array(np.meshgrid(sigma_real, p_strength, p_in_strength)).T.reshape(-1, 3)
sort_comb = np.argsort(combination.sum(axis=1))
idx_filter_list = []
sigma_strength_list = []
p_strength_list = []
p_in_strength_list = []
# =========================================================================
for i in sort_comb:
    file = base_name.format(i + 1)
    try:
        with open(file, 'rb') as f:
            data = pl.load(f)
        gene_posterior_probabilities.append(data.pi_avg_print[data.gibbs_iter - 1])
        posterior_betas.append(max(abs(data.posterior_beta)))
        idx_filter_list.append(combination[i].sum())
        sigma_strength_list.append(combination2[i][0])
        p_strength_list.append(combination2[i][1])
        p_in_strength_list.append(combination2[i][2])
        if flag:
            gene_index = data.genes_extension_index
            num_genes = data.num_genes
            flag = False
    except IOError:
        print("file does not exist")

gene_posterior_probabilities_to_plot = np.asarray(gene_posterior_probabilities).sum(axis=1)
samples = gene_posterior_probabilities_to_plot.size
posterior_betas_to_plot = np.asarray(posterior_betas)
column_values = list(gene_index.keys())
idx_filter = np.asarray(idx_filter_list)
# ======================================================
TOOLTIPS = [
    ("(x,y)", "(@x, @y)"),
    ("Sigma strength:", "@sigma"),
    ("p strength:", "@p"),
    ("p_in strength:", "@p_in"),
    ("Gene prob:", "@gene_prob"),
    ("Gene names:", "@gene_name")
]
source = ColumnDataSource(data=dict(
    x=list(idx_filter),
    y=list(gene_posterior_probabilities_to_plot),
    sigma=sigma_strength_list,
    p=p_in_strength_list,
    p_in=p_in_strength_list,
    gene_prob=gene_posterior_probabilities,
    gene_name=[column_values] * samples
))
# ======================================================
# set output to static HTML file
output_file(filename=options.out_file)
# create a new plot with a specific size
tools = ['box_select', 'box_zoom', 'pan', 'tap', 'wheel_zoom', 'undo', 'redo', 'reset', 'save', 'hover']
p = figure(sizing_mode="stretch_width", max_width=2000, height=500,
           x_axis_label='Combined strength', y_axis_label='Sum of probabilities',
           title=options.out_file, tooltips=TOOLTIPS, tools=tools)
p.line(idx_filter, gene_posterior_probabilities_to_plot, line_width=1)
p.circle('x', 'y', source=source, size=5)
# save the results to a file
save(p)
