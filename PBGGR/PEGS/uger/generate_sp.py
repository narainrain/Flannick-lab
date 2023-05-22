parameters = []
parameters.append("--gene-file /home/unix/ellamas/magma/geneloc37.tsv")
parameters.append("--sumstats-file /home/unix/ellamas/project1/data/ldl_study.tbl")
parameters.append("--sumstats-id-col MarkerName")
parameters.append("--sumstats-beta-col Effect")
parameters.append("--sumstats-se-col StdErr")
parameters.append("--sumstats-p-col GC.Pvalue")
parameters.append("--sumstats-n-col Weight")
parameters.append("--gene-id-col 1")
parameters.append("--gene-chrom-col 2")
parameters.append("--gene-start-col 3")
parameters.append("--gene-end-col 4")
parameters.append("--all-positive")
parameters.append("--dentist-thr 218.70")
parameters.append("--n-threshold 0.95")
parameters.append("--p 0.00061442")
parameters.append("--heritability 0.1212")
parameters.append("--heritability-causal 0.9")
parameters.append("--max-component-size 50")
parameters.append("--SNP-log-diff-ignore 0")
parameters.append("--causal-window 0")
parameters.append("--SNP-model Greedy_SNP")
parameters.append("--cluster --block-size 1000000")
parameters.append("--debug-level 1")

ld_loc = "--ld-file /humgen/diabetes2/users/lthakur/LDSC/ldsc/1000G_EUR_Phase3_plink/1000G.EUR.QC.{}_row.ld"
dentist_loc = "--dentist-file /home/unix/ellamas/dentist/results/chr{}.DENTIST.txt"
scripts_path = "/humgen/diabetes2/users/ellamas/PEGS/BGAM/bin"
script_name = "/humgen/diabetes2/users/ellamas/PEGS/BGAM/bin/gene_bf.py"

for i in range(1, 23):
    file_name = "sp{}.sh".format(i)
    ld = (ld_loc + " \\\n").format(i)
    dentist = (dentist_loc + " \\\n").format(i)
    my_range = "--range {}:1-0 \n".format(i)
    f = open(file_name, "w")
    f.write("source /broad/software/scripts/useuse\n"
            "reuse Python-3.6\n"
            "source /home/unix/ellamas/test/bin/activate\n"
            "export PATH=$PATH:" + scripts_path + "\n"
            "python \\\n" + script_name + " \\\n")
    for parameter in parameters:
        f.write(parameter + " \\\n")
    f.write(ld)
    f.write(dentist)
    f.write(my_range)
    f.close()
