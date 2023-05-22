source /broad/software/scripts/useuse
use Python-3.6
source /home/unix/ellamas/test/bin/activate
export PATH=$PATH:/humgen/diabetes2/users/ellamas/PEGS/BGAM/bin
python \
/humgen/diabetes2/users/ellamas/PEGS/BGAM/bin/gene_bf.py \
--gene-file /home/unix/ellamas/magma/geneloc37.tsv \
--sumstats-file /humgen/diabetes2/users/ellamas/pheno_data/ldl_study.tbl \
--sumstats-id-col MarkerName \
--sumstats-beta-col Effect \
--sumstats-se-col StdErr \
--sumstats-p-col GC.Pvalue \
--sumstats-n-col Weight \
--gene-id-col 1 \
--gene-chrom-col 2 \
--gene-start-col 3 \
--gene-end-col 4 \
--all-positive \
--dentist-thr 218.70 \
--n-threshold 0.5 \
--p 0.00061442 \
--heritability 0.1212 \
--heritability-causal 0.9 \
--max-component-size 10 \
--SNP-log-diff-ignore 0 \
--causal-window 0 \
--SNP-model Greedy_SNP \
--cluster --block-size 1000000 \
--debug-level 1 \
--ld-file /humgen/diabetes2/users/ellamas/PEGS/BGAM/data/ld/HMGCR_ld_rows.ld \
--dentist-file /humgen/diabetes2/users/ellamas/PEGS/BGAM/data/dentist/chr5.DENTIST.txt \
--chrom 5 \
--job-id 75 \
--model greedy \
--add-cov-uncertainty
