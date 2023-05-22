import math
from optparse import OptionParser
import pandas as pd
import pickle
import json

# arguments
parser = OptionParser("usage: %prog [options]")
parser.add_option("", "--locus-file", default=None)
parser.add_option("", "--out-folder", default=None)
(options, _) = parser.parse_args()

region_size = 10000000
locus = pd.read_csv(options.locus_file, sep='\t')
locus = locus.astype({'bp': 'int32'})
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
for chrom in chrom_to_len:
    regions = math.ceil(chrom_to_len[chrom] / region_size)
    for region in range(regions):
        start = region * region_size
        end = (region+1) * region_size
        filter_pd = locus.loc[(start <= locus['bp']) & (locus['bp'] <= end) & (locus['Chr'] == int(chrom))]
        my_dict = {}
        my_dict2 = {}
        my_dict['variant'] = list(filter_pd['Chr'].astype(str) + ":" + filter_pd['bp'].astype(str) + "_" + filter_pd['A1'].astype(str) + "/" +
            filter_pd['A2'].astype(str))
        my_dict['position'] = list(filter_pd['bp'])
        my_dict['log_pvalue'] = list(filter_pd['p'].apply(lambda x: -math.log10(x) if x != 0 else 0))

        my_dict2['variant'] = list(filter_pd['Chr'].astype(str) + ":" + filter_pd['bp'].astype(str) + "_" + filter_pd['A2'].astype(str) + "/" +
            filter_pd['A1'].astype(str))
        my_dict2['position'] = list(filter_pd['bp'])
        my_dict2['log_pvalue'] = list(filter_pd['p'].apply(lambda x: -math.log10(x) if x != 0 else 0))

        output = options.out_folder  + "." + str(chrom) + "." + str(region) + ".json"
        output2 = options.out_folder + "." + str(chrom) + "." + str(region) + ".json.inv"

        filehandler = open(output, 'wt')
        filehandler.write(json.dumps(my_dict))

        filehandler2 = open(output2, 'wt')
        filehandler2.write(json.dumps(my_dict2))

output3 = options.out_folder + ".json"
filehandler3 = open(output3, 'wt')
filehandler3.write("Done")
