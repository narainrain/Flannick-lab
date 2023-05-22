from bokeh.plotting import figure, output_file, save
from optparse import OptionParser
import numpy as np
import pandas as pd

# arguments
parser = OptionParser("usage: %prog [options]")
parser.add_option('--input-file', action='append')
parser.add_option("", "--out-file", default=None)
(options, _) = parser.parse_args()
# ======================================================
# prepare some data
data = []
for file in options.input_file:
    data = data + pd.read_csv(file, sep='\t', usecols=['Weight'])['Weight'].to_list()
# set output to static HTML file
output_file(filename=options.out_file)

# create a new plot with a specific size
p = figure(sizing_mode="stretch_width", max_width=2000, height=500, x_axis_label='Sample size', y_axis_label='points',
           title=options.out_file)

# data = np.random.normal(0, 1, 1000)
arr_hist, edges = np.histogram(data, bins=200)
p.quad(top=arr_hist, bottom=0, left=edges[:-1], right=edges[1:],
       fill_color="navy", line_color="white", alpha=0.5)
# p.line(x, np.log10(y), legend_label="Enrichment", line_width=1)

# save the results to a file
save(p)
