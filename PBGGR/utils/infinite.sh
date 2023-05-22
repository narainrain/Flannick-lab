#! /bin/bash

source /broad/software/scripts/useuse;
use .python-3.9.2;
cd /humgen/diabetes/users/ellamas/PEGS_pipeline/source/PEGS/PEGS/bin;
export PATH=$PATH:/humgen/diabetes/users/ellamas/PEGS_pipeline/source/PEGS/PEGS/bin;
source /humgen/diabetes2/users/ellamas/enviroments/env_3_9_2/bin/activate;
python test.py --gene-file /humgen/diabetes/users/ellamas/PEGS_pipeline/raw/gene_pos/gene_loc_37.tab --gene-id-col 1 --gene-chrom-col 2 --gene-start-col 3 --gene-end-col 4 --sumstats-file /humgen/diabetes/users/ellamas/PEGS_pipeline/raw/pheno/ldl_study.tbl.gz --tabix /broad/software/free/Linux/redhat_7_x86_64/pkgs/tabix/tabix_0.2.6/bin/tabix --tmp-folder ../data/tmp/ --sumstats-id-col MarkerName --sumstats-beta-col Effect --sumstats-se-col StdErr --sumstats-p-col GC.Pvalue --sumstats-n-col Weight --sumstats-z-col GC.Zscore --dentist-folder /humgen/diabetes/users/ellamas/PEGS_pipeline/raw/dentist/ --ld-folder /humgen/diabetes/users/ellamas/PEGS_pipeline/raw/ld/ --ld-chrom1-col "#CHR_A" --all-positive --n-threshold 0.5 --dentist-thr 218.7 --range 5:74000000-75000000 --extension 500000 --max-component-size 1000000 --debug-level 2 --causal-window 100 --gene-prior 0.05 --p 0.0004 --heritability 0.1212 --total_num_SNPs 84700000 --window-model one_window --out-folder ../data/out/test --sample-size 85000;
