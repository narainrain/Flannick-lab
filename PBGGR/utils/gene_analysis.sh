#!/bin/bash
while getopts m:b:o:n:l:g:e: flag
do
    case "${flag}" in
        m) magma_file=${OPTARG};;
        b) bf_file=${OPTARG};;
        o) out=${OPTARG};;
        n) genes=${OPTARG};;
        l) locus_file=${OPTARG};;
        g) gene_file=${OPTARG};;
        e) exome_file=${OPTARG};;
    esac
done
# variables used in the analysis ---------------------------------------------------------------------------------------
# exome_file="/humgen/diabetes2/users/ellamas/LDL/LDL_exome/min_p" # exome are not sorted
# ----------------------------------------------------------------------------------------------------------------------
genes_plus=$((genes+1))
sort -gk3 $bf_file |  awk -v num_gen="$genes" 'NR <= num_gen {print $1}' > top.bf
sort -gk9 $magma_file | awk -v num_gen="$genes_plus" 'NR > 1 && NR <= num_gen {print $1}' > top.magma
sort top.* | uniq -d > top.both
sort top.{magma,bf,bf} | uniq -u > top.unique_magma
sort top.{magma,magma,bf} | uniq -u > top.unique_bf
sort -gk9 $magma_file | awk '{print $1 "\t" NR "\t" $9}' | fgrep -wf top.unique_bf - | sort -k1 > aux1
sort -gk3 $bf_file | awk '{print $1 "\t" NR "\t" $3}' | fgrep -wf top.unique_bf - | sort -k1 > aux2  # warning fgrep add every prefix match
tail -n +2 $exome_file | sort -gk2 | awk '{print $1 "\t" NR "\t" $2}' | fgrep -wf top.unique_bf - | sort -k1 > aux3
awk '{print $1 "\t" $2 "\t" $3 "\t" $4}' $gene_file | fgrep -wf top.unique_bf - | sort -k1 > aux4
join aux1 aux2 | join - aux3 | join - aux4 > aux5
my_path_aux="${bf_file%/*}/"
my_path=${my_path_aux#*humgen/diabetes2/users/ellamas/pipeline/out/projects/}
path_for_plots="/home/unix/ellamas/private_html/bgam_out/projects/locuszoom/PEGS_plots/${my_path}"
mkdir -p $path_for_plots
awk -v the_path=$my_path '{new_var= "https://internal.broadinstitute.org/~ellamas/bgam_out/projects/locuszoom/?chrom="$8"&start="$9"&end="$10"&file="the_path""$1".json"; print $0 " " new_var}' aux5 > aux6
echo -e "gene magma_pos magma_p bf_pos bf exome_pos exome_p chr start end link" | cat - aux6 > ${out}
rm aux*
rm top*

source /humgen/diabetes2/users/ellamas/enviroments/env_3_8_3/bin/activate
python /humgen/diabetes2/users/ellamas/pipeline/source/PEGS/utils/gwas_to_json.py \
--locus-file $locus_file \
--analysis-file $out \
--out-folder $path_for_plots


