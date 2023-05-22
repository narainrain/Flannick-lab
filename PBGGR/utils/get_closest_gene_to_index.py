import math

import pandas as pd
from optparse import OptionParser
import numpy as np


def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx


# arguments
parser = OptionParser("usage: %prog [options]")
parser.add_option("", "--locus-file", default=None)
parser.add_option("", "--genes-file", default=None)
parser.add_option("", "--out", default=None)
parser.add_option("", "--region-size", type="int", default=1000000)
parser.add_option("", "--p-value-th", type="int", default=8)
(options, _) = parser.parse_args()

locus = pd.read_csv(options.locus_file, sep='\t')
genes = pd.read_csv(options.genes_file, sep='\t')
genes.columns = ["gene", "chrom", "start", "end", "type", "unknown"]
genes['TSS'] = np.where(genes['type'] == '+', genes['start'], genes['end'])
chrom_to_len = {'1': 250000000,
                '2': 250000000,
                '3': 210000000,
                '4': 200000000,
                '5': 190000000,
                '6': 180000000,
                '7': 160000000,
                '8': 150000000,
                '9': 140000000,
                '10': 140000000,
                '11': 140000000,
                '12': 140000000,
                '13': 120000000,
                '14': 110000000,
                '15': 110000000,
                '16': 100000000,
                '17': 90000000,
                '18': 80000000,
                '19': 70000000,
                '20': 70000000,
                '21': 50000000,
                '22': 50000000
                }
selected_genes = []
for chrom in chrom_to_len:
    iterations = math.ceil(chrom_to_len[chrom] / options.region_size)
    for job in range(iterations):
        start = job * options.region_size
        end = (job+1) * options.region_size
        filter1 = locus.loc[locus['Chr'] == int(chrom)]
        filter1 = filter1.astype({'bp': 'int32'})
        filter2 = filter1.loc[start <= filter1['bp']]
        filter3 = filter2.loc[filter2['bp'] <= end]
        filter4 = filter3.loc[filter3['p'] > 0]
        filter4['log_p'] = filter4['p'].apply(lambda x: -math.log10(x) if x != 0 else 0)
        nonzero = filter4.loc[filter4['log_p'] >= options.p_value_th]
        if not nonzero.empty:
            pos_of_min_p = nonzero.loc[nonzero['log_p'].idxmax()]['bp']
            gene_filter1 = genes.loc[genes['chrom'] == chrom]
            idx_gene = find_nearest(gene_filter1['TSS'].values, pos_of_min_p)
            closest_gene = gene_filter1.iloc[idx_gene]['gene']
            selected_genes.append(closest_gene)

with open(options.out, 'w') as fp:
    for item in selected_genes:
        # write each item on a new line
        fp.write("%s\n" % item)
