import pickle as pl
import matplotlib.pyplot as plt
import matplotlib
import glob
from matplotlib.widgets import Slider, Button
import numpy
import numpy as np
import math

matplotlib.use('Qt5Agg')


def bar_plot(name, gene_prob_by_file):
    n = len(gene_prob_by_file[0])
    values = gene_prob_by_file.mean(axis=1)
    # Figure Size
    fig, ax = plt.subplots()

    # Horizontal Bar Plot
    ax.barh(name, values)

    # Remove axes splines
    for s in ['top', 'bottom', 'left', 'right']:
        ax.spines[s].set_visible(False)

    # Remove x, y Ticks
    ax.xaxis.set_ticks_position('none')
    ax.yaxis.set_ticks_position('none')

    # Add padding between axes and labels
    ax.xaxis.set_tick_params(pad=5)
    ax.yaxis.set_tick_params(pad=10)

    # Add x, y gridlines
    ax.grid(color='grey',
            linestyle='-.', linewidth=0.5,
            alpha=0.2)

    # Show top values
    ax.invert_yaxis()

    # Add annotation to bars
    for i in ax.patches:
        plt.text(i.get_width(), i.get_y() + 0.5,
                 str(round((i.get_width()), 2)),
                 fontsize=10, fontweight='bold',
                 color='grey')

    # Add Plot Title
    ax.set_title('Gene probability from {} files and {} samples'.format(len(gene_prob_by_file[0]), len(window_list)),
                 loc='left', )
    # Create 1 axes for 1 sliders red
    axred = plt.axes([0.25, 0.05, 0.65, 0.02])

    # Create a slider from 0.0 to 1.0 in axes axred
    # with 1.0 as initial value.
    red = Slider(axred, 'Percentage of Files', 0.0, 1.0, 1.0)

    # Create function to be called when slider value is changed

    def update(val):
        r = red.val
        num_samples = int(math.ceil(n * red.val))
        data = gene_prob_by_file[:, 0:num_samples].mean(axis=1)
        ax.barh(name, data)
        plt.title("Files: {}".format(num_samples))

    # Call update function when slider value is changed
    red.on_changed(update)

    # Create axes for reset button and create button
    resetax = plt.axes([0.8, 0.0, 0.1, 0.04])
    button = Button(resetax, 'Reset', color='gold', hovercolor='skyblue')

    # Create a function resetSlider to set slider to
    # initial values when Reset button is clicked

    def resetSlider(event):
        red.reset()

    # Call resetSlider function when clicked on reset button
    button.on_clicked(resetSlider)
    plt.show()


def plot_hist(samples):
    # Create a subplot
    fig, ax = plt.subplots()
    plt.subplots_adjust(bottom=0.35)
    n = len(samples)
    plt.title("window distribution")
    my_text = fig.text(0.02, 0.1, "")
    r = 1.0
    # Create and plot a bar chart

    plt.hist(samples, bins=200)

    # Create 3 axes for 3 sliders red,green and blue
    axred = plt.axes([0.25, 0.2, 0.65, 0.03])
    axred2 = plt.axes([0.25, 0.15, 0.65, 0.03])

    # Create a slider from 0.0 to 1.0 in axes axred
    # with 0.6 as initial value.
    red = Slider(axred, 'Remove samples end', 0.0, 1.0, 1.0)
    remove = Slider(axred2, 'Remove samples start', 0.0, 1.0, 0.0)

    # Create function to be called when slider value is changed

    def update(val):
        r = red.val
        rem = remove.val
        num_samples = int(math.ceil(n * red.val))
        num_samples_rem = int(math.ceil(n * rem))
        data = samples[num_samples_rem:num_samples]
        ax.hist(data, bins=200)
        my_text.set_text("Total samples: {}. Removed start: {}. "
                         "Removed end: {}. Reminding: {}".format(n, num_samples_rem, n - num_samples,
                                                                 num_samples-num_samples_rem))

    # Call update function when slider value is changed
    red.on_changed(update)
    remove.on_changed(update)

    # Create axes for reset button and create button
    resetax = plt.axes([0.8, 0.025, 0.1, 0.04])
    button = Button(resetax, 'Reset', color='gold', hovercolor='skyblue')

    # Create a function resetSlider to set slider to
    # initial values when Reset button is clicked

    def resetSlider(event):
        red.reset()
        remove.reset()

    # Call resetSlider function when clicked on reset button
    button.on_clicked(resetSlider)
    # plt.title("Samples: {}".format(n), loc='left')
    # Display graph
    plt.show()


# ===== Parameters ==========================================
folder = '/Users/ellamas/PycharmProjects/PEGS/PEGS/data/out/'
plot_mod = "sorted"  # sorted, unsorted
# ===== Parameters ==========================================
names = glob.glob(folder + "test.data*")
if plot_mod == "sorted":
    iter_list = []
    for name in names:
        iter_list.append(int(name.split('/')[-1].split('_')[2]))
    sorted_idx = np.argsort(iter_list)
else:
    sorted_idx = range(len(names))
window_list = []
gene_list = None
gene_names = None
first = True
for idx in sorted_idx:
    file = names[idx]
    with open(file, 'rb') as f:
        data = pl.load(f)
    stop = data.gibbs_iter + data.burn_in
    window_list += list(data.window_mu_vec_print[:, 0:stop].flatten(order='F'))
    if first:
        gene_list = data.pi_avg_print[0:stop-data.burn_in, :].T
        gene_names = numpy.empty(data.num_genes, dtype=object)
        for gene in data.genes_extension_index:
            gene_names[data.genes_extension_index[gene]] = gene
        first = False
    else:
        gene_list = np.concatenate((gene_list, data.pi_avg_print[0:stop-data.burn_in, :].T), axis=1)
plot_hist(window_list)
bar_plot(gene_names, gene_list)
