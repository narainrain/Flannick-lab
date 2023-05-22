from bokeh.plotting import figure, output_file, save
from optparse import OptionParser
import numpy as np
from bokeh.models import Span


# arguments
parser = OptionParser("usage: %prog [options]")
parser.add_option("", "--magma-input-file", default=None)
parser.add_option("", "--bf-input-file", default=None)
parser.add_option("", "--out-file", default=None)
(options, _) = parser.parse_args()

# ======================================================
# prepare some data
xm, ym = np.loadtxt(options.magma_input_file, delimiter='\t', skiprows=1, usecols=(0, 1), unpack=True)
xb, yb = np.loadtxt(options.bf_input_file, delimiter='\t', skiprows=1, usecols=(0, 1), unpack=True)

# set output to static HTML file
output_file(filename=options.out_file)

# create a new plot with a specific size
p = figure(sizing_mode="stretch_width", max_width=2000, height=500, x_axis_label='Rank', y_axis_label='log of p', title=options.out_file)

# add a circle renderer
p.line(xm, np.log10(ym), legend_label="Magma", color="red", line_width=1)
p.line(xb, np.log10(yb), legend_label="Bayes Factor", color="blue", line_width=1)

sum_m = np.sum(-np.log10(ym))
magma_enr = np.sum(-(np.log10(ym)/sum_m)*np.log10(ym))

sum_xm = np.sum(xm)
magma_enr2 = np.sum(((max(xm)-xm+1)/sum_xm)*np.log10(ym))

sum_b = np.sum(-np.log10(yb))
bf_enr = np.sum(-(np.log10(yb)/sum_b)*np.log10(yb))

sum_xb = np.sum(xb)
bf_enr2 = np.sum(((max(xb)-xb+1)/sum_xb)*np.log10(yb))

# Vertical line
vline = Span(location=1000, dimension='height', line_color='grey', line_width=2)

hlinem = Span(location=magma_enr, dimension='width', line_color='red', line_width=2)
hlineb = Span(location=bf_enr, dimension='width', line_color='blue', line_width=2)

hlinem2 = Span(location=magma_enr2, dimension='width', line_dash='dashed', line_color='red', line_width=2)
hlineb2 = Span(location=bf_enr2, dimension='width', line_dash='dashed', line_color='blue', line_width=2)
p.renderers.extend([vline, hlinem, hlineb, hlinem2, hlineb2])
# circle = p.circle(x, y, fill_color="red", size=15)

# save the results to a file
save(p)
