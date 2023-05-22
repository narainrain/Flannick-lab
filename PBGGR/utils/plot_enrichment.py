from bokeh.plotting import figure, output_file, save
from optparse import OptionParser
import numpy as np

# arguments
parser = OptionParser("usage: %prog [options]")
parser.add_option("", "--input-file", default=None)
parser.add_option("", "--out-file", default=None)
(options, _) = parser.parse_args()

# ======================================================
# prepare some data
x, y = np.loadtxt(options.input_file, delimiter='\t', skiprows=1, usecols=(0, 1), unpack=True)

# set output to static HTML file
output_file(filename=options.out_file)

# create a new plot with a specific size
p = figure(sizing_mode="stretch_width", max_width=2000, height=500, x_axis_label='Rank', y_axis_label='log of p', title=options.out_file)

# add a circle renderer
p.line(x, np.log10(y), legend_label="Enrichment", line_width=1)
# circle = p.circle(x, y, fill_color="red", size=15)

# save the results to a file
save(p)
