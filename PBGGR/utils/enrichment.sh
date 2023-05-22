#!/bin/bash
# variables used in gene enrichment analysis -----------------------------------------------------------------------
file1="/humgen/diabetes2/users/ellamas/LDL/LDL_results/s_body_max_a_1_0kb.txt"
p_value_file="/humgen/diabetes2/users/ellamas/LDL/LDL_exome/min_p_reg.ldl.ukbb.tsv"
out="/humgen/diabetes2/users/ellamas/LDL/LDL_min_p_GE/ge_body_max_a_1_0kb.txt"
# ----------------------------------------------------------------------------------------------------------------------
source /broad/software/scripts/useuse
use Python-3.6
source /home/unix/ellamas/test/bin/activate
head -n 4000 $file1 | cut -f1 | \
python /home/unix/ellamas/project1/bin/exomes_enrichment.py \
--gene-list-file /dev/stdin --p-value-file $p_value_file \
--p-value-file-gene-col 1 --p-value-file-p-col 2 --debug-level 2 \
--gene-stats-file /humgen/diabetes/t2d_exomes/pipeline/enrichment/gene_stats.tsv --gene-stats-file-gene-col 1 \
--gene-stats-file-stats-col 2 --gene-stats-file-stats-col 3 --gene-stats-file-stats-col 4 \
--gene-stats-file-stats-col 5 --num-match 50 \
--progressive-mode > $out
