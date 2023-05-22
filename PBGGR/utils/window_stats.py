from optparse import OptionParser
import numpy as np
from bokeh.plotting import figure, output_file, save
from bokeh.models import Label
from bokeh.layouts import Row, column
from bokeh.models.widgets import Div

# arguments
parser = OptionParser("usage: %prog [options]")
parser.add_option('--input-file', action='append')
parser.add_option("", "--out-file", default=None)
(options, _) = parser.parse_args()

# ======================================================
window_vec = np.zeros(len(options.input_file))
pip_avg_vec = np.zeros(len(options.input_file))
i = 0
for file in options.input_file:
    window_vec[i] = int(file.split('_')[0].split('/')[-1].split('w')[1])
    pip_avg_vec[i] = np.average(np.loadtxt(file, delimiter=' ', usecols=(1), unpack=True))
    i += 1

coeff = np.polyfit(window_vec, pip_avg_vec, 1)
coeff2 = np.polyfit(window_vec, pip_avg_vec, 2)
x = np.linspace(min(window_vec), max(window_vec), 100)
y = coeff[0]*x + coeff[1]
y2 = coeff2[0]*(x**2) + coeff2[1]*x + coeff2[2]

# Max value calculation
max_x = (-coeff2[1]) / (2*coeff2[0])
max_y = coeff2[0]*(max_x**2) + coeff2[1]*max_x + coeff2[2]
# set output to static HTML file
output_file(filename=options.out_file)

# create a new plot with a specific size
p = figure(width=1900, height=500, x_axis_label='window size', y_axis_label='posterior probability', title=options.out_file)

# add a circle renderer
p.circle(window_vec, pip_avg_vec, legend_label="Points", line_width=5)
p.line(x, y, legend_label="Degree 1", line_width=2, color='orange')
p.line(x, y2, legend_label="Degree 2", line_width=2, color='green')
p.circle(max_x, max_y, legend_label="Max", color='red', line_width=5)

div = Div(
    text=""" <b>Max Value</b> """ + str(max_x),
    width=2000, height=100)

# save the results to a file
save(column(p, div))
