from bokeh.plotting import figure, output_file, save
from optparse import OptionParser
import numpy as np
from bokeh.io import show
from bokeh.plotting import figure
from bokeh.models import CustomJS, CheckboxGroup, Div
from bokeh.layouts import Row, column
from bokeh.palettes import Category20_20
from random import choice


# arguments
parser = OptionParser("usage: %prog [options]")
parser.add_option('--input-magma', action='append')
parser.add_option('--input-bf', action='append')
parser.add_option("", "--out-file", default=None)
(options, _) = parser.parse_args()

# ======================================================
num_files_magma = len(options.input_magma)
num_files_bf = len(options.input_bf)

bf_enr_vec = []
magma_enr_vec = []
bf_enr_vec2 = []
magma_enr_vec2 = []
for file_bf in options.input_bf:
    x, y = np.loadtxt(file_bf, delimiter='\t', skiprows=1, usecols=(0, 1), unpack=True)
    sum_m = np.sum(-np.log10(y))
    enr = np.sum(-(np.log10(y) / sum_m) * np.log10(y))
    bf_enr_vec.append(enr)
    sum_x = np.sum(x)
    bf_enr2 = np.sum(((max(x)-x+1)/sum_x)*np.log10(y))
    bf_enr_vec2.append(bf_enr2)

for file_magma in options.input_magma:
    x, y = np.loadtxt(file_magma, delimiter='\t', skiprows=1, usecols=(0, 1), unpack=True)
    sum_m = np.sum(-np.log10(y))
    enr = np.sum(-(np.log10(y) / sum_m) * np.log10(y))
    magma_enr_vec.append(enr)
    sum_x = np.sum(x)
    magma_enr2 = np.sum(((max(x)-x+1)/sum_x)*np.log10(y))
    magma_enr_vec2.append(magma_enr2)

x1 = np.asarray(bf_enr_vec)
x2 = np.asarray(magma_enr_vec)
x3 = np.asarray(bf_enr_vec2)
x4 = np.asarray(magma_enr_vec2)
x_label = ['Bayes Factor\nweight\non p-value', 'Magma\nweight\non p-value', 'Bayes Factor\nweight\non rank', 'Magma\nweight\non rank']

q1 = np.asarray([np.quantile(x1, 0.25, axis=0), np.quantile(x2, 0.25, axis=0), np.quantile(x3, 0.25, axis=0), np.quantile(x4, 0.25, axis=0)])
q2 = np.asarray([np.quantile(x1, 0.5, axis=0), np.quantile(x2, 0.5, axis=0), np.quantile(x3, 0.5, axis=0),  np.quantile(x4, 0.5, axis=0)])
q3 = np.asarray([np.quantile(x1, 0.75, axis=0), np.quantile(x2, 0.75, axis=0), np.quantile(x3, 0.75, axis=0), np.quantile(x4, 0.75, axis=0)])
min =np.asarray([np.quantile(x1, 0, axis=0), np.quantile(x2, 0, axis=0), np.quantile(x3, 0, axis=0),    np.quantile(x4, 0, axis=0)])
max =np.asarray([np.quantile(x1, 1, axis=0), np.quantile(x2, 1, axis=0), np.quantile(x3, 1, axis=0),    np.quantile(x4, 1, axis=0)])

p = figure(tools="", background_fill_color="#efefef", x_range=x_label, toolbar_location=None)

# stems
p.segment(x_label, max, x_label, q3, line_color="black")
p.segment(x_label, min, x_label, q1, line_color="black")

# boxes
p.vbar(x_label, 0.7, q2, q3, fill_color="#E08E79", line_color="black")
p.vbar(x_label, 0.7, q1, q2, fill_color="#3B8686", line_color="black")

# whiskers (almost-0 height rects simpler than segments)
p.rect(x_label, min, 0.2, 0.01, line_color="black")
p.rect(x_label, max, 0.2, 0.01, line_color="black")

p.xgrid.grid_line_color = None
p.ygrid.grid_line_color = "white"
p.grid.grid_line_width = 2
p.xaxis.major_label_text_font_size = "16px"

output_file(filename=options.out_file)
save(p)
