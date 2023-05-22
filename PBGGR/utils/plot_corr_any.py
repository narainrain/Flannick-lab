import math

import pandas as pd
from bokeh.plotting import figure, output_file, save
from optparse import OptionParser
import numpy as np
from bokeh.io import show
from bokeh.plotting import figure, ColumnDataSource
from bokeh.models import CustomJS, CheckboxGroup, Div, TapTool, OpenURL, ColorBar
from bokeh.layouts import Row, column
from bokeh.palettes import Category20_20, Spectral6
from bokeh.transform import linear_cmap
from random import choice
import itertools


# arguments
parser = OptionParser("usage: %prog [options]")
parser.add_option('--input-file', action='append')
parser.add_option("", "--gene-file", default=None)
parser.add_option("", "--out-file", default=None)
parser.add_option("", "--name1", default=0)
parser.add_option("", "--data1", default=1)
parser.add_option("", "--header1", default=0)
parser.add_option("", "--name2", default=0)
parser.add_option("", "--data2", default=1)
parser.add_option("", "--header2", default=1)

(options, _) = parser.parse_args()

# ======================================================
region_size = 10000000
num_files = len(options.input_file)
TOOLTIPS = [
    ("(x,y)", "(@x, @y)"),
    ("Gene:", "@gene"),
    ("Diff:", "@diff"),
    ("Chr:start-end:", "@chrom:@start-@end"),
    ("Region:", "@region")
]
x_name = options.input_file[0]
y_name = options.input_file[1]
tools = ['box_select', 'box_zoom', 'pan', 'tap', 'wheel_zoom', 'undo', 'redo', 'reset', 'save', 'hover']
fig = figure(width=1900, height=500, x_axis_label=x_name, y_axis_label=y_name,
             title='Correlation', tooltips=TOOLTIPS, tools=tools)
args = []
code = ""
i = 0
names = []
names.append(x_name + " -VS- " + y_name)
data1 = pd.read_csv(x_name, sep='\t')
data1.columns = ['gene', 'prob', 'p', 'num']
data1['gene'] = data1['gene'].str.strip()
data2 = pd.read_csv(y_name, sep='\t')
aux = data1.set_index('gene').join(data2.set_index('gene'), how='inner', lsuffix='_A', rsuffix='_B')
aux['gene'] = aux.index
gene_df = pd.read_csv(options.gene_file, sep='\t', header=0)
gene_df.columns = ['gene', 'chr', 'start', 'end', 'dir', 'num']
join = aux.set_index('gene').join(gene_df.set_index('gene'), how='left', lsuffix='_A', rsuffix='_B')
join['gene'] = join.index
join['diff'] = abs(join['prob']-join['posterior_gene_prob'])
join = join.fillna(0)
x1 = join['prob'].values
x2 = join['posterior_gene_prob'].values
genes = join['gene'].values
chrom = join['chr'].values
start = join['start'].values
end = join['end'].values
diff = join['diff'].values
region = np.floor(start/region_size)
source = ColumnDataSource(data=dict(
    x=list(x1),
    y=list(x2),
    gene=list(genes),
    diff=list(diff),
    chrom=list(chrom),
    start=list(start),
    end=list(end),
    region=list(region),
))
# glyph = fig.scatter(x1, x2, color=choice(Category20_20), line_width=5)
mapper = linear_cmap(field_name='diff', palette=Spectral6, low=min(diff), high=max(diff))
color_bar = ColorBar(color_mapper=mapper['transform'], width=8)
fig.add_layout(color_bar, 'right')
glyph = fig.circle('x', 'y', line_color=mapper, color=mapper, fill_alpha=1, size=20, source=source)
args += [('glyph' + str(i), glyph)]
code += "glyph{}.visible = this.active.includes({});".format(i, i)
i += 1
# aux = x_name.split("projects/")[1].split("/")
# file_path = aux[0] + '/' + aux[1] + '/' + aux[2] + '/Magma/plots_data/magma1.@chrom.@region.json'
# locus_path = "https://internal.broadinstitute.org/~ellamas/bgam_out/projects/locuszoom/" \
#              "gene_loc.html?chrom=@chrom&start=@start&end=@end&file="
# aux2 = x_name.split("projects/")[1].split("/")
# file_path2 = aux2[0] + '/' + aux2[1] + '/' + aux2[2] + '/Magma/plots_data/magma1.@chrom.@region.json'
# url = locus_path + file_path + "&file2=" + file_path2
url = 'asdfad'
taptool = fig.select(type=TapTool)
taptool.callback = OpenURL(url=url)

code += "console.log('checkbox_button_group: active=' + this.active);"
# checkbox = CheckboxGroup(width=1900, height=500, labels=[file for file in names],
#                          active=[0])
callback = CustomJS(args={key: value for key, value in args}, code=code)
# checkbox.js_on_click(callback)
output_file(filename=options.out_file)
save(column(fig))
