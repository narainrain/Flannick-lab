# mount server folder in local machine
sshfs ellamas@dig-ae-dev-02:/humgen/diabetes2/users/ellamas Broad

#pipeline
source ~flannick/bin/setup.sh
perl /humgen/diabetes2/users/ellamas/lap/trunk/bin/run.pl --meta /humgen/diabetes2/users/ellamas/lap/projects/pegs/meta/pegs.meta --only ldl --only s_1 --only window_1 --init --mkdir
perl /humgen/diabetes2/users/ellamas/lap/trunk/bin/run.pl --meta /humgen/diabetes2/users/ellamas/lap/projects/pegs/meta/pegs.meta --only ldl --only s_1 --only window_1 --check
perl /humgen/diabetes2/users/ellamas/lap/trunk/bin/run.pl --meta /humgen/diabetes2/users/ellamas/lap/projects/pegs/meta/pegs.meta --only ldl --only s_1 --only window_1 --bsub
perl /humgen/diabetes2/users/ellamas/lap/trunk/bin/run.pl --meta /humgen/diabetes2/users/ellamas/lap/projects/pegs/meta/pegs.meta --only ldl --check --update-cmd-key-different .

# gene enrichment
head -n 4000 /humgen/diabetes2/users/ellamas/LDL/LDL_results/s_body_max_p_0kb.txt | cut -f1 | \
python /home/unix/ellamas/project1/bin/exomes_enrichment.py \
--gene-list-file /dev/stdin --p-value-file /humgen/diabetes2/users/ellamas/exome_results/exome_LoF.txt \
--p-value-file-gene-col 1 --p-value-file-p-col 9 --debug-level 2 \
--gene-stats-file /humgen/diabetes/t2d_exomes/pipeline/enrichment/gene_stats.tsv --gene-stats-file-gene-col 1 \
--gene-stats-file-stats-col 2 --gene-stats-file-stats-col 3 --gene-stats-file-stats-col 4 \
--gene-stats-file-stats-col 5 --num-match 50 \
--progressive-mode >> /humgen/diabetes2/users/ellamas/LDL/LDL_LoF_GE/ge_body_max_p_0kb.txt


# our model run in local
python \
gene_bf.py \
--gene-file ../data_chr2/ref_gene.bed \
--sumstats-file ../data_chr2/ldl_study.tbl \
--ld-file ../data_chr2/ld_rows.r.gz \
--range 2:43054411-44233145

# LDSC regression runs in local machine
# clean the GWAS data
/Users/ellamas/ldsc/munge_sumstats.py \
--snp SNP \
--N-col N \
--a1 A1 \
--a2 A2 \
--p p \
--signed-sumstats b,0 \
--sumstats /Users/ellamas/ldsc/data/bmi.txt \
--out /Users/ellamas/ldsc/data/bmi
# Execute h_2
/Users/ellamas/ldsc/ldsc.py \
--h2 /Users/ellamas/ldsc/data/bmi.sumstats.gz \
--ref-ld-chr /Users/ellamas/ldsc/data/eur_w_ld_chr/ \
--w-ld-chr /Users/ellamas/ldsc/data/eur_w_ld_chr/ \
--out /Users/ellamas/ldsc/data/bmi_h2


# Magma
/home/unix/ellamas/magma/magma --annotate window=0 \
--snp-loc /humgen/diabetes2/users/lthakur/LDSC/ldsc/1000G_EUR_Phase3_plink/g1000_eur.bim \
--gene-loc /humgen/diabetes2/users/lthakur/LDSC/ldsc/1000G_EUR_Phase3_plink/geneloc37.tsv \
--out /humgen/diabetes2/users/lthakur/LDSC/ldsc/1000G_EUR_Phase3_plink/results/0kb

/home/unix/ellamas/magma/magma --bfile /humgen/diabetes2/users/lthakur/LDSC/ldsc/1000G_EUR_Phase3_plink/g1000_eur \
--pval /home/unix/ellamas/pheno_data/ldl-allSNPs.ma ncol=N \
--gene-annot /humgen/diabetes2/users/lthakur/LDSC/ldsc/1000G_EUR_Phase3_plink/results/0kb.genes.annot \
--out /humgen/diabetes2/users/lthakur/LDSC/ldsc/1000G_EUR_Phase3_plink/results/0kb

# plink ld calculation
ish -pe smp 2 -binding linear:2 -l h_vmem=16G
/home/unix/ellamas/cojo/plink --r2 inter-chr --memory 14400 --bfile /humgen/diabetes2/users/lthakur/LDSC/ldsc/1000G_EUR_Phase3_plink/1000G.EUR.QC.6 --out /humgen/diabetes2/users/lthakur/LDSC/ldsc/1000G_EUR_Phase3_plink/1000G.EUR.QC.6_row

# tabix
tabix /humgen/diabetes2/users/ryank/data/dbsnp/All_20180423.vcf.gz chr:pos-pos

# dentist
/home/unix/ellamas/dentist/DENTIST_1.1.0.0 --gwas-summary /home/unix/ellamas/pheno_data/ldl-allSNPs.ma --bfile /humgen/diabetes2/users/lthakur/LDSC/ldsc/1000G_EUR_Phase3_plink/1000G.EUR.QC.2 --out /home/unix/ellamas/chr2
#cluster option 1
qsub -pe smp 1 -binding linear:1 -l h_vmem=1G sp1.sh

# cluster option 2
use UGER
qsub -t 1-4 sp2.sh
# cluster option 3
use UGER
qsub -t 1-5 -pe smp 1 -binding linear:1 -l h_vmem=2G -o /home/unix/ellamas/project1/out/dentist/chr1/ sp.sh

# interactive node
ish -pe smp 1 -binding linear:1 -l h_vmem=16G

# our model run in server
python \
/home/unix/ellamas/project1/bin/gene_bf.py \
--gene-file /home/unix/ellamas/project1/data_chr2/ref_gene.bed \
--sumstats-file /home/unix/ellamas/project1/data_chr2/ldl_study.tbl \
--ld-file /home/unix/ellamas/cojo/data_chrom2/chrom2_rows.ld \
--dentist-file /home/unix/ellamas/cojo/data_chrom2/dentist_chr2.txt \
--all-positive \
--debug-level 2 \
--range 2:21000000-21500000



