#!/bin/bash
cut -f2 jointGwasMc_LDL.txt | cut -d':' -f1 | cut -d'r' -f2 | paste - jointGwasMc_LDL.txt > aux1
cut -f3 aux1 | cut -d':' -f2 | paste - aux1 > aux2
awk '{print $2":"$1":"toupper($6)":"toupper($7)"\t"$2"\t"$1"\t"$5"\t"toupper($6)"\t"toupper($7)"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12}' aux2 > aux3
tail -n +2 aux3 | sort -k2nr -k3g > aux4
echo -e "#variant\tchrom\tpos\tMarkerName\tAllele1\tAllele2\tEffect\tStdErr\tWeight\tGC.Pvalue\tminor_AF" > header
cat header aux4 > ldlMetaboChip_study.tbl
rm aux*
bgzip -c ldlMetaboChip_study.tbl > ldlMetaboChip_study.tbl.gz
tabix -s2 -b3 -e3 ldlMetaboChip_study.tbl.gz
