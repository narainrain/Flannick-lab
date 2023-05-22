#!/bin/bash
while getopts w: flag
do
    case "${flag}" in
        w) window=${OPTARG};;
    esac
done
/home/unix/flannick/links/targeted/lib/magma/magma --annotate window=$window --snp-loc snp.loc --gene-loc gene.loc --out "${window}kb";
# synonyms is the only difference but in the analysis those are automatically detected
/home/unix/flannick/links/targeted/lib/magma/magma --bfile /humgen/diabetes2/users/lthakur/LDSC/ldsc/1000G_EUR_Phase3_plink/g1000_eur --gene-annot "${window}kb.genes.annot" --pval snp.pval N=188577 --out "${window}kb";
# =================== GET GENE NAMES FROM NUMBER ID ===================================
sort -k1 "${window}kb.genes.out" > aux1;
sort -k1 gene.loc > aux2;
join aux2 aux1 > aux3;
cut -d' ' -f6-14 aux3 > aux4;
awk -v OFS='\t' '{ $1=$1; print }' aux4 > aux5;
cat header aux5 > "${window}kb.magma";
rm aux*;
# ======================================================================================
source /broad/software/scripts/useuse;
use .python-3.8.3;
source /humgen/diabetes2/users/ellamas/enviroments/env_3_8_3/bin/activate;
sort -gk9 "${window}kb.magma" | awk 'NR > 1 && NR <= 1000 {print $0}' | cut -f1 -d" " | python /humgen/diabetes2/users/ellamas/lap/projects/targeted/bin/exomes_enrichment.py --gene-list-file /dev/stdin --p-value-file exomes.ldl  --p-value-file-gene-col 1 --p-value-file-p-col 2 --debug-level 2 --gene-stats-file /humgen/diabetes/t2d_exomes/pipeline/enrichment/gene_stats.tsv --gene-stats-file-gene-col 1 --gene-stats-file-stats-col 2 --gene-stats-file-stats-col 3 --gene-stats-file-stats-col 4 --gene-stats-file-stats-col 5 --num-match 50 --progressive-mode > "${window}kb.magma.enr";
python /humgen/diabetes2/users/ellamas/pipeline/source/PEGS/utils/plot_enrichment.py --input-file "${window}kb.magma.enr" --out-file ~ellamas/private_html/tmp/"${window}kb.magma.enr.plot";
