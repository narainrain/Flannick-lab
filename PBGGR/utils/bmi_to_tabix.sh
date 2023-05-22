join -1 1 -2 2 bmi.sort /humgen/diabetes2/users/ryank/data/dbsnp/All_20180423.rsid.sort > aux1
awk -v OFS='\t' '{ $1=$1; print }' aux1 > aux2

cut -f9 aux2 | cut -d':' -f2 | paste - aux2 > aux3
cut -f9 aux2 | cut -d':' -f1 | paste - aux3 > aux4
sort -k1nr -k2g aux4 > aux5
echo -e "#chrom\tpos\tMarkerName\tAllele1\tAllele2\tFreq1.Hapmap\tEffect\tStdErr\tGC.Pvalue\tWeight\tvariant" > header
cat header aux5 > bmi_study.tbl
bgzip -c bmi_study.tbl > bmi_study.tbl.gz
tabix -s1 -b2 -e2 bmi_study.tbl.gz
