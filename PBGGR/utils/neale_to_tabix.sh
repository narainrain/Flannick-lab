#!/bin/bash
while getopts i:o: flag
do
    case "${flag}" in
        i) file=${OPTARG};;
        o) outfile=${OPTARG};;
    esac
done
cut -f1 $file | cut -d':' -f3 | paste - $file > aux1
cut -f1 $file | cut -d':' -f2 | paste - aux1 > aux2
cut -f1 $file | cut -d':' -f1 | paste - aux2 > aux3
tail -n +2 aux3 | sort -k4 > aux3.sort
join -1 4 -2 1 aux3.sort /humgen/diabetes2/users/ryank/data/dbsnp/All_20180423.rsid.sort > aux4
sort -k2nr -k3g aux4 > aux5
echo -e "#variant chrom pos Allele1 Allele2 minor_AF low_confidence_variant Weight AC ytx Effect StdErr GC.Zscore GC.Pvalue MarkerName" > header
cat header aux5 > aux6
awk -v OFS='\t' '{ $1=$1; print }' aux6 > $outfile
rm aux*
bgzip -c $outfile > "${outfile}.gz"
tabix -s2 -b3 -e3 "${outfile}.gz"
