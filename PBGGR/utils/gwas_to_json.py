import math
from optparse import OptionParser
import pandas as pd
import pickle
import json

# arguments
parser = OptionParser("usage: %prog [options]")
parser.add_option("", "--locus-file", default=None)
parser.add_option("", "--analysis-file", default=None)
parser.add_option("", "--out-folder", default=None)
(options, _) = parser.parse_args()

analysis = pd.read_csv(options.analysis_file, sep=' ')
locus = pd.read_csv(options.locus_file, sep='\t')
ext = 600000
for index, row in analysis.iterrows():
    gene = row['gene']
    chr = row['chr']
    start = int(row['start'])
    end = int(row['end'])
    print(gene)
    my_dict = {}
    my_dict2 = {}
    filter1 = locus.loc[locus['Chr'] == chr]
    filter1 = filter1.astype({'bp': 'int32'})
    filter2 = filter1.loc[(start - ext) <= filter1['bp']]
    filter3 = filter2.loc[filter2['bp'] <= (end + ext)]

    my_dict['variant'] = list( str(chr) + ":" + filter3['bp'].astype(str) + "_" + filter3['A1'].astype(str) + "/" +
        filter3['A2'].astype(str))
    my_dict['position'] = list(filter3['bp'])
    my_dict['log_pvalue'] = list(filter3['p'].apply(lambda x: -math.log10(x) if x != 0 else 0))

    my_dict2['variant'] = list( str(chr) + ":" + filter3['bp'].astype(str) + "_" + filter3['A2'].astype(str) + "/" +
        filter3['A1'].astype(str))
    my_dict2['position'] = list(filter3['bp'])
    my_dict2['log_pvalue'] = list(filter3['p'].apply(lambda x: -math.log10(x) if x != 0 else 0))

    output = options.out_folder + gene + ".json"
    output2 = options.out_folder + "inv_" + gene + ".json"

    filehandler = open(output, 'wt')
    filehandler.write(json.dumps(my_dict))

    filehandler2 = open(output2, 'wt')
    filehandler2.write(json.dumps(my_dict2))
