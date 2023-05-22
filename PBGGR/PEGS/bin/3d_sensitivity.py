from optparse import OptionParser
import pickle as pl
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.ticker import LinearLocator
import numpy as np
import matplotlib
matplotlib.use('Qt5Agg')

# arguments
parser = OptionParser("usage: %prog [options]")
parser.add_option("", "--input-folder", default=None)
parser.add_option("", "--out-file", default=None)
(options, _) = parser.parse_args()

# ======================================================
gene_to_plot = 3
for p_val in range(10):
    # ==== sigma to plot ====
    # p_val = 0
    # ======================================================
    base_name = options.input_folder + "sensitivity_idx_{}.data"
    gene_posterior_probabilities = np.zeros((10, 10))
    X = np.zeros((10, 10))
    Y = np.zeros((10, 10))
    posterior_betas = []
    flag = True
    gene_index = None
    num_genes = None
    # sorted priors ==========================================================
    num_conf = 10
    sigma_strength = np.linspace(1, 2000, num_conf)
    p_strength = np.logspace(-2, 3, num_conf)
    p_in_strength = np.logspace(-2, 3, num_conf)
    combination = np.array(np.meshgrid(sigma_strength, p_strength, p_in_strength)).T.reshape(-1, 3)
    my_filter = list(np.nonzero(combination[:, 1] == p_strength[p_val])[0])
    # =========================================================================
    j = 0
    k = 0
    for i in my_filter:
        file = base_name.format(i + 1)
        try:
            with open(file, 'rb') as f:
                data = pl.load(f)
            gene_posterior_probabilities[k, j] = data.pi_avg_print[data.gibbs_iter - 1][gene_to_plot]  # .sum()
            X[k, j] = combination[i][0]
            Y[k, j] = np.log10(combination[i][2])
            k += 1
            if k == 10:
                k = 0
                j += 1
            if flag:
                gene_index = data.genes_extension_index
                num_genes = data.num_genes
                flag = False
        except IOError:
            print("file does not exist")

    column_values = list(gene_index.keys())

    # ======================================================
    # X1 = np.linspace(1, 10000, num_conf)
    # Y1 = np.log10(np.logspace(-2, 3, num_conf))
    # X1, Y1 = np.meshgrid(X1, Y1)
    fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
    plt.title("P strength: " + str(p_strength[p_val]) + " Gene: " + str(column_values[gene_to_plot]))
    # Plot the surface.
    surf = ax.plot_surface(X, Y, gene_posterior_probabilities, cmap=cm.coolwarm,
                           linewidth=0, antialiased=False)

    # Customize the z axis.
    # ax.set_zlim(-1.01, 1.01)
    ax.zaxis.set_major_locator(LinearLocator(10))
    # A StrMethodFormatter is used automatically
    ax.zaxis.set_major_formatter('{x:.02f}')

    # Add a color bar which maps values to colors.
    fig.colorbar(surf, shrink=0.5, aspect=5)
    ax.set_xlabel('sigma')
    ax.set_ylabel('p_in')
    ax.set_zlabel('sum prob', rotation=60)
plt.show()
