#!/bin/bash
while getopts c: flag
do
    # shellcheck disable=SC2220
    case "${flag}" in
        c) chr=${OPTARG};;
    esac
done
ldd="/humgen/diabetes2/users/lthakur/LDSC/ldsc/1000G_EUR_Phase3_plink/ldl_dentist"
ld="/humgen/diabetes2/users/lthakur/LDSC/ldsc/1000G_EUR_Phase3_plink/1000G.EUR.QC.${chr}_row.ld"
out="/humgen/diabetes2/users/lthakur/LDSC/ldsc/1000G_EUR_Phase3_plink/dentist_ld_ldl_${chr}.txt"
echo "Chromosome: $chr"
awk '{print $3}' "$ld" > s1
awk '{print $6}' "$ld" > s2
cat s1 s2 | sort | uniq -d > ld.txt
cut -f1 $ldd > cuted_ldd
sort ld.txt cuted_ldd cuted_ldd | uniq -u > ld_no
fgrep -v -f ld_no $ld > aux
rm cuted_ldd
rm ld.txt
rm ld_no
rm s1
rm s2
awk '{print $3}' aux > l1
awk '{print $6}' aux > l2
cat l1 l2 | sort | uniq -d > keep
fgrep -f keep $ldd > $out
rm aux
rm l1
rm l2
rm keep
