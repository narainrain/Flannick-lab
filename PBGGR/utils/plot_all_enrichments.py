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
parser.add_option('--input-file', action='append')
parser.add_option("", "--out-file", default=None)
(options, _) = parser.parse_args()

# ======================================================
num_files = len(options.input_file)
fig = figure(width=1900, height=500, x_axis_label='Rank', y_axis_label='log of p', title='Enrichment of all options')
args = []
code = ""
i = 0
for file in options.input_file:
    x, y = np.loadtxt(file, delimiter='\t', skiprows=1, usecols=(0, 1), unpack=True)
    glyph = fig.line(x, np.log10(y), color=choice(Category20_20), line_width=1)
    args += [('glyph' + str(i), glyph)]
    code += "glyph{}.visible = this.active.includes({});".format(i, i)
    i += 1

code += "console.log('checkbox_button_group: active=' + this.active);"
checkbox = CheckboxGroup(width=1900, height=500, labels=[file for file in options.input_file],
                         active=[1])
callback = CustomJS(args={key: value for key, value in args}, code=code)
checkbox.js_on_click(callback)
output_file(filename=options.out_file)
save(column(fig, checkbox))
