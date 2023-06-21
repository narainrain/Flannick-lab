"""
Simulation Generation
====================================
This class is used to simulate generation of beta values for snps given LD scores, different methods
to determine causal snps: assigning from true data, random generation from distribution(s), random generation
The run method is the only method needed to run the entire process of generating LD scores which
will be retrieved from LD input files. Then after generation of the LD scores, then the true beta scores
will be generated from multiple options: generation from (generate_causal_genes_or_snps) to determine which
model will be used to generate causal snps. Then marginal snps will be calculated from
the true betas and then stored in 2 files: locus files and a gwas file

"""

# Author: Sumit Narain <narain@broadinstitute.org>
# Author: Alex Ellamas <ellamas@broadinstitute.org>
# Author: Jason Flannick <flannick@broadinstitute.org>
# License:

import sys
import numpy as np
import scipy.stats
import gzip
import subprocess
import math


def bail(message):
    sys.stderr.write("%s\n" % message)
    sys.exit(1)


def warn(message):
    sys.stderr.write("Warning: %s\n" % message)


def log(message):
    sys.stderr.write("%s\n" % message)


def open_gz(file, sort_col=None, reverse=True, header=False):
    if sort_col is None:
        if file[-3:] == ".gz":
            fh = gzip.open(file)
            fh.readline()
        else:
            fh = open(file)
            fh.readline()
    else:
        reverse_flag = ""
        if reverse:
            reverse_flag = "r"

        if file[-3:] == ".gz":
            return subprocess.Popen(
                ["zcat %s | tail -n+%d | sort -g%sk%s" % (file, 1 if header else 2, reverse_flag, sort_col)],
                stdout=subprocess.PIPE, shell=True).stdout
        else:
            return subprocess.Popen(
                ["tail -n+%d %s | sort -g%sk%s" % (1 if header else 2, file, reverse_flag, sort_col)],
                stdout=subprocess.PIPE, shell=True).stdout


# This method is used to stablise the LD matrix when doing a cholensky operation.
def matrix_stabilization(matrix):
    # stabilization achieved by multiply the diagonal by a small factor to avoid singular matrix
    matrix_size = matrix.shape[0]  # this is a square matrix
    null_cov_diag = np.diag(matrix)
    stab_term_frac = 1e-4
    stab_matrix = np.diag(np.ones(matrix_size)) * (np.mean(np.diag(matrix)) * stab_term_frac)
    null_cov_matrix = matrix + stab_matrix
    # eigen value decomposition to assure the resulting correlation matrix is positive semi-definite
    try:
        l, Q = np.linalg.eigh(null_cov_matrix)
    except np.linalg.linalg.LinAlgError:
        return null_cov_matrix
    l[l < 0] = 0
    null_cov_matrix = np.dot(Q, np.dot(np.diag(l), Q.transpose())) + stab_matrix
    # transformation to restore the values changed in the stabilization step
    null_cov_norm_matrix = np.diag(np.sqrt(null_cov_diag / np.diag(null_cov_matrix)))
    null_cov_matrix = np.dot(null_cov_norm_matrix, np.dot(null_cov_matrix, null_cov_norm_matrix))
    return null_cov_matrix


class Simulation:

    def __init__(self, options):
        self.options = options

    # Only one function here that runs the entire process of generation of LD to generation of Marginal betas
    def run(self):
        # variables  -----------------------------------------------------------------------------------------------------------
        assume_r = self.options.assume_r
        ld_file = self.options.ld_file
        # ----------------------------------------------------------------------------------------------------------------------

        if ld_file is None:
            bail("Need --ld-file")

        if self.options.h2 < 0 or self.options.h2 > 1:
            bail("--h2 must be between 0 and 1")

        print("Reading ld file... ", ld_file)
        ld_rev_fh = open_gz(ld_file, 7, True)
        ld_for_fh = open_gz(ld_file, 7, False)
        snp_to_component = {}
        component_to_snp = {}
        index = 0
        component = 0
        max_component_size = 100
        chr_to_snp_pos = {}
        snp_to_chr_pos = {}
        snp_to_alt = {}
        snp_to_ref = {}
        chr_pos_to_snp = {}

        line_rev = None
        valid_rev = True
        line_for = None
        valid_for = True
        count = 0

        # Used to calculate the variance later
        total_number_snps = 0
        snp_container = []

        # generation of LD matrix will start first. Obtain data from the LD file and BIM file for SNP information
        while True:
            count += 1

            # read through the file forwards and backwards at the same time, taking the highest abs
            if valid_rev and line_rev is None:
                line_rev = ld_rev_fh.readline()
                cols_rev = line_rev.strip().split()
                if len(cols_rev) != 7:
                    log("Bad line number " + str(count))
                    break
                value_rev = float(cols_rev[6])
                if value_rev > 0:
                    valid_rev = True
                else:
                    valid_rev = False
            if valid_for and line_for is None:
                line_for = ld_for_fh.readline()
                cols_for = line_for.strip().split()
                if len(cols_for) != 7:
                    log("Bad line number " + str(count))
                    break
                value_for = float(cols_for[6])
                if value_for < 0:
                    valid_for = True
                else:
                    valid_for = False

            if valid_rev and (not valid_for or abs(value_rev) >= abs(value_for)):
                line = line_rev
                cols = cols_rev
                value = value_rev
                line_rev = None
            elif valid_for and (not valid_rev or abs(value_rev) < abs(value_for)):
                line = line_for
                cols = cols_for
                value = value_for
                line_for = None
            else:
                break

            if value < 0:
                assume_r = True

            if abs(value) < self.options.ld_threshold:  # if is not, I think it will also continue
                continue
            snp_1 = cols[2]
            snp_2 = cols[5]
            if type(snp_1) == bytes:
                snp_1 = bytes.decode(snp_1)
            if type(snp_2) == bytes:
                snp_2 = bytes.decode(snp_2)

            snp_1_chr = bytes.decode(cols[0])
            snp_1_pos = int(cols[1])
            snp_to_chr_pos[snp_1] = (snp_1_chr, snp_1_pos)
            snp_2_chr = bytes.decode(cols[3])
            snp_2_pos = int(cols[4])
            snp_to_chr_pos[snp_2] = (snp_2_chr, snp_2_pos)

            if snp_1_chr not in chr_to_snp_pos:
                chr_to_snp_pos[snp_1_chr] = set()
            chr_to_snp_pos[snp_1_chr].add(snp_1_pos)
            chr_pos_to_snp[(snp_1_chr, snp_1_pos)] = snp_1
            if snp_2_chr not in chr_to_snp_pos:
                chr_to_snp_pos[snp_2_chr] = set()
            chr_to_snp_pos[snp_2_chr].add(snp_2_pos)
            chr_pos_to_snp[(snp_2_chr, snp_2_pos)] = snp_2

            if snp_1 not in snp_to_component and snp_2 not in snp_to_component:
                component += 1
                snp_to_component[snp_1] = component
                snp_to_component[snp_2] = component
                component_to_snp[component] = set()
                component_to_snp[component].add(snp_1)
                component_to_snp[component].add(snp_2)
            elif snp_1 in snp_to_component and snp_2 not in snp_to_component:
                if len(component_to_snp[snp_to_component[snp_1]]) < self.options.max_component_size:
                    snp_to_component[snp_2] = snp_to_component[snp_1]
                    component_to_snp[snp_to_component[snp_1]].add(snp_2)
                else:
                    component += 1
                    snp_to_component[snp_2] = component
                    component_to_snp[component] = set()
                    component_to_snp[component].add(snp_2)
            elif snp_2 in snp_to_component and snp_1 not in snp_to_component:
                if len(component_to_snp[snp_to_component[snp_2]]) < self.options.max_component_size:
                    snp_to_component[snp_1] = snp_to_component[snp_2]
                    component_to_snp[snp_to_component[snp_2]].add(snp_1)
                else:
                    component += 1
                    snp_to_component[snp_1] = component
                    component_to_snp[component] = set()
                    component_to_snp[component].add(snp_1)
            elif snp_2 in snp_to_component and snp_1 in snp_to_component and snp_to_component[snp_1] != \
                    snp_to_component[
                        snp_2]:
                if len(component_to_snp[snp_to_component[snp_1]]) + len(
                        component_to_snp[snp_to_component[snp_2]]) <= self.options.max_component_size:
                    component_1 = snp_to_component[snp_1]
                    component_2 = snp_to_component[snp_2]
                    for snp in component_to_snp[component_2]:
                        snp_to_component[snp] = component_1
                    component_to_snp[component_1] = component_to_snp[component_1].union(component_to_snp[component_2])
                    component_to_snp.pop(component_2)

            if len(component_to_snp[snp_to_component[snp_1]]) > max_component_size:
                max_component_size = len(component_to_snp[snp_to_component[snp_1]])
        log("Max component size: %s" % max_component_size)

        ld_for_fh.close()
        ld_rev_fh.close()
        if self.options.bim_file:
            bim_fh = open(self.options.bim_file)
            for line in bim_fh:
                total_number_snps = total_number_snps + 1
                cols = line.strip().split()
                snp = cols[1]
                chr = cols[0]
                pos = int(cols[3])

                snp_to_alt[snp] = str(cols[4])
                snp_to_ref[snp] = str(cols[5])

                snp_container.append(snp)

                if snp not in snp_to_chr_pos:
                    snp_to_chr_pos[snp] = (chr, pos)
                    if chr not in chr_to_snp_pos:
                        chr_to_snp_pos[chr] = set()
                    chr_to_snp_pos[chr].add(pos)
                    if (chr, pos) not in chr_pos_to_snp:
                        chr_pos_to_snp[(chr, pos)] = snp

                    if snp not in snp_to_component:
                        component += 1
                        snp_to_component[snp] = component
                        component_to_snp[component] = set()
                        component_to_snp[component].add(snp)

            bim_fh.close()

        for chrom in chr_to_snp_pos:
            chr_to_snp_pos[chrom] = sorted(list(chr_to_snp_pos[chrom]))
        snp_to_index = {}
        component_to_cor = {}
        for component in component_to_snp:
            index = 0
            # print "%s %s" % (component, component_to_snp[component])
            for snp in sorted(component_to_snp[component]):
                snp_to_index[snp] = index
                index += 1
            # print("%s %s" % (component, len(component_to_snp[component])))
            component_to_cor[component] = np.identity(len(component_to_snp[component]), dtype=np.float64)
        ld_fh = open_gz(ld_file, 7)
        for line in ld_fh:
            cols = line.strip().split()
            value = float(cols[6])
            if abs(value) < self.options.ld_threshold:
                continue
            if not assume_r:
                value = math.sqrt(value)
            snp_1 = cols[2]
            snp_2 = cols[5]
            if type(snp_1) == bytes:
                snp_1 = bytes.decode(snp_1)
            if type(snp_2) == bytes:
                snp_2 = bytes.decode(snp_2)

            if snp_to_component[snp_1] == snp_to_component[snp_2]:
                component_to_cor[snp_to_component[snp_1]][snp_to_index[snp_1], snp_to_index[snp_2]] = value
                component_to_cor[snp_to_component[snp_1]][snp_to_index[snp_2], snp_to_index[snp_1]] = value
        ld_fh.close()

        snp_to_true_beta_params = {}

        # Generate the true betas based on variance and proportion fraction
        # This beta generation occurs with a distribution file
        if self.options.beta_dist_file:

            import random
            import bisect

            beta_dist_fh = open(self.options.beta_dist_file)
            for line in beta_dist_fh:
                cols = line.strip().split()
                if len(cols) != 6 and len(cols) != 7:
                    warn(
                        "Ignoring line without six columns for (chrom, start, end, p, mean, var [replicate]):  %s" % line.strip())
                    continue
                chrom = cols[0]
                start = int(cols[1])
                end = int(cols[2])
                p = float(cols[3])
                if p < 0 or p > 1:
                    log("Error: bad value for bernoulli (%s)" % (p))
                    continue
                mean = float(cols[4])
                var = float(cols[5])
                if var < 0:
                    log("Error: bad value for var (%s)" % (var))
                    continue

                if len(cols) == 7:
                    rep = int(cols[6])
                    if rep < 1:
                        log("Error: bad value for rep (%s)" % (rep))
                        continue
                else:
                    rep = 1

                # find all overlapping snps

                if chrom in chr_to_snp_pos:
                    start_ind = bisect.bisect_left(chr_to_snp_pos[chrom], start)
                    for ind in range(start_ind, len(chr_to_snp_pos[chrom])):
                        if chr_to_snp_pos[chrom][ind] < start:
                            log("There is a bug -- %s is less than %s" % (chr_to_snp_pos[chrom][ind], start))
                            continue
                        if chr_to_snp_pos[chrom][ind] > end:
                            break
                        cur_snp = chr_pos_to_snp[(chrom, chr_to_snp_pos[chrom][ind])]
                        if cur_snp not in snp_to_true_beta_params:
                            snp_to_true_beta_params[cur_snp] = []
                        snp_to_true_beta_params[cur_snp].append((p, mean, var, rep - 1))

        snp_to_true_beta = {}
        snp_container_causal_for_genome = []
        count_of_all_causal_snps_genome = 0

        # there should be a generated input file for all the causal snps for an entire genome
        if self.options.input_causal_snps:

            snps_fh = open(self.options.input_causal_snps)

            # skips the header line
            next(snps_fh)

            for line in snps_fh:

                cols = line.strip().split()
                chromosome = cols[0]
                variant = cols[1]
                count_of_all_causal_snps_genome = count_of_all_causal_snps_genome + 1

                if chromosome == self.options.chrom:
                    snp_container_causal_for_genome.append(variant)

        if self.options.beta_file:
            beta_fh = open(self.options.beta_file)
            for line in beta_fh:
                cols = line.strip().split()
                if len(cols) != 2:
                    warn("Ignoring line without two columns for (snp, beta):  %s" % line.strip())
                    continue
                snp = cols[0]
                if snp in snp_to_index:
                    true_beta = float(cols[1])
                    snp_to_true_beta[snp] = true_beta
            beta_fh.close()

        snp_to_beta = {}
        snp_to_p = {}
        locus_fh = open(str(self.options.output_file) + ".locus", 'w')
        gwas_fh = open(str(self.options.output_file) + ".gwas", 'w')

        if self.options.ldsc_format:
            sys.stdout.write("SNP\tReplicate\tA1\tA2\tZ\tN\n")
        else:
            header_locus = "SNP\tChrom\tPos\tRef\tAlt\tReplicate\tEffect\tStdErr\tP-value\n"
            header_gwas = "MarkerName\tChrom\tPos\tRef\tAlt\tWeight\tGC.Zscore\tGC.Pvalue\tOverall\tDirection\tEffect\tStdErr\n"
            locus_fh.write(header_locus)
            gwas_fh.write(header_gwas)

        num_true_causal_snps = {}
        snp_true_beta_container = {}

        # start of the simulation by generating true betas and calculating marginals
        # The first part is to determine how much of this
        if True:

            import random

            portion_of_causal_snps = 0.0
            float(portion_of_causal_snps)

            # This is used in case the p provided is 0
            # the p cannot be 0 because it is used to calculate the heritability
            # This will change the value of p depending on the number of causal snps that are currently stored

            if float(self.options.p) == 0:

                if len(snp_container_causal_for_genome) == 0:
                    portion_of_causal_snps = 0.0000000001
                    print("There were no causal snps in the gene")
                else:
                    portion_of_causal_snps = (
                                float(count_of_all_causal_snps_genome) / float(self.options.total_snp_genome))
                    print("This is the proportion of causal snps: ", portion_of_causal_snps)

            else:
                portion_of_causal_snps = self.options.p

            # heritability is spread depending on the size of the chromosome, so each snp will have an equal h2
            mean = 0
            h2_proportion = self.options.h2 / (self.options.total_snp_genome * portion_of_causal_snps)

            # Bottom heritability is used for single chromosome testing
            # h2_proportion = self.options.h2 / (total_number_snps * portion_of_causal_snps)

            var = h2_proportion

            # Generation of causal snps either by causal snps input file or by randomisation
            if self.options.input_causal_snps:
                snp_container_causal = snp_container_causal_for_genome
            else:
                causal_snps = math.ceil(total_number_snps * portion_of_causal_snps)
                random.shuffle(snp_container)

                if causal_snps == 0:
                    snp_container_causal = []
                else:
                    snp_container_causal = snp_container[0:int(causal_snps - 1)]

            # random sampling per each component
            for component in component_to_cor:

                cor_matrix = component_to_cor[component]

                if np.linalg.det(cor_matrix) == 0:
                    cor_matrix = matrix_stabilization(cor_matrix)

                # Associate each snp with a true beta
                for it in range(0, self.options.num_sim):

                    if it not in num_true_causal_snps:
                        num_true_causal_snps[it] = 0

                    cur_true_beta = np.zeros(len(component_to_snp[component]))
                    cur_geno_var = np.ones(len(component_to_snp[component])) * (1 - self.options.h2)

                    # If data was obtained from beta_dist_file
                    if len(snp_to_true_beta_params) > 0:
                        for snp in component_to_snp[component]:
                            if snp in snp_to_true_beta_params:
                                for true_beta_params in snp_to_true_beta_params[snp]:
                                    if true_beta_params[3] is None or true_beta_params[3] == it:
                                        if random.random() <= true_beta_params[0]:
                                            cur_true_beta[snp_to_index[snp]] = \
                                                np.random.normal(loc=true_beta_params[1],
                                                                 scale=np.sqrt(true_beta_params[2]),
                                                                 size=1)[0]
                                            snp_true_beta_container[snp] = cur_true_beta[snp_to_index[snp]]
                                        else:
                                            snp_true_beta_container[snp] = 0
                                        break

                    # If no beta_dist_file then use the mean,variance from h2 provided and sample size
                    else:
                        for snp in component_to_snp[component]:
                            if snp in snp_container_causal:
                                cur_true_beta[snp_to_index[snp]] = np.random.normal(mean, scale=np.sqrt(var), size=1)[0]
                                snp_true_beta_container[snp] = cur_true_beta[snp_to_index[snp]]
                            else:
                                snp_true_beta_container[snp] = 0

                    # if the true beta was provided in a different file; it will overwrite any snps assigned before
                    if len(snp_to_true_beta) > 0:
                        for snp in component_to_snp[component]:
                            if snp in snp_to_true_beta:
                                cur_true_beta[snp_to_index[snp]] = snp_to_true_beta[snp]
                                snp_true_beta_container[snp] = cur_true_beta[snp_to_index[snp]]
                            else:
                                snp_true_beta_container[snp] = 0

                                # keeps track of all the causal snps per simulation
                    num_true_causal_snps[it] += np.sum(cur_true_beta != 0)

                    try:
                        L = np.linalg.cholesky(cor_matrix)
                    except np.linalg.linalg.LinAlgError:
                        l, Q = np.linalg.eigh(cor_matrix)
                        l[l < 0] = 0
                        L = np.dot(Q, np.diag(np.sqrt(l)))
                        # normalize diagonal to be 1
                        D2 = np.diag(1 / np.sqrt(np.diag(np.dot(L, L.transpose()))))
                        L = np.dot(D2, L)

                    # noise term
                    v = np.random.normal(loc=0, scale=1, size=len(component_to_snp[component]))

                    # creating a distribution of the size parameter
                    size = np.random.laplace(loc=self.options.N, scale=2000)

                    S = ((np.matmul(cor_matrix, cur_true_beta) * np.sqrt(size)) / np.sqrt(cur_geno_var)) + np.matmul(L,
                                                                                                                     v)
                    B = S * np.sqrt(cur_geno_var / size)

                    # Another way to calculate marginal betas
                    # marginal_beta = (np.matmul(cor_matrix, cur_true_beta)) + np.matmul(L, v)
                    # standard_error = np.sqrt(cur_geno_var / size)
                    # Z_score = marginal_beta / standard_error

                    # cor_matrix_2 = np.identity(len(cor_matrix.diagonal()))
                    # Z_score = ((np.matmul(cor_matrix_2, cur_true_beta) * np.sqrt(size)) / np.sqrt(
                    #    cur_geno_var)) + np.matmul(L, v)
                    # marginal_beta = Z_score * np.sqrt(cur_geno_var / size)

                    for snp in component_to_snp[component]:
                        if len(component_to_snp[component]) == 1:
                            beta = B[0]
                        else:
                            beta = B[snp_to_index[snp]]

                        se = np.sqrt(cur_geno_var[snp_to_index[snp]] / size)
                        z_score = float(beta) / float(se)
                        pvalue = str(2 * scipy.stats.norm.sf(abs(z_score)))

                        if self.options.ldsc_format:
                            sys.stdout.write("%s\t%d\t%s\t%s\t%.3g\t%d\n" % (
                                snp, it + 1, "R", "A", beta / se, size))
                        else:
                            chrom = str(snp_to_chr_pos[snp][0])
                            pos = str(snp_to_chr_pos[snp][1])
                            alt = snp_to_alt[snp]
                            ref = snp_to_ref[snp]
                            N = str(size)

                            if snp in snp_true_beta_container:
                                temp_true_beta_holder = snp_true_beta_container[snp]
                                snp_true_beta_container[snp] = [snp, chrom, temp_true_beta_holder, pos, ref, alt]

                            line_gwas = snp + "\t" + chrom + "\t" + pos + "\t" + ref + "\t" + alt + "\t" + N + "\t" + str(
                                z_score) + "\t" + pvalue + "\t" + "+\t+-+-+\t" + str(
                                beta) + "\t" + str(se)
                            line_locus = snp + "\t" + chrom + "\t" + pos + "\t" + ref + "\t" + alt + "\t" + str(
                                it + 1) + "\t" + str(beta) + "\t" + str(se) + "\t" + pvalue
                            gwas_fh.write(line_gwas + "\n")
                            locus_fh.write(line_locus + "\n")

            # Store all the snps and associated betas into a file. Sort them by the absolute value of the betas
            file_output = open(str(self.options.output_file) + ".true_beta", "w")
            file_output.write("SNP\tChrom\tTrue Beta\tPos\tRef\tAlt\n")
            for key, value in sorted(snp_true_beta_container.items(), key=lambda k: abs(k[1][2]), reverse=True):
                file_output.write(
                    '{}\t{}\t{}\t{}\t{}\t{}\n'.format(value[0], value[1], value[2], value[3], value[4], value[5]))
            file_output.close()

        if self.options.num_causal_snps_out is not None:
            out_suma_fh = open(self.options.num_causal_snps_out, 'w')
            out_suma_fh.write("Replicate\tNum_Causal\n")
            for it in range(0, self.options.num_sim):
                out_suma_fh.write("%d\t%d\n" % (it + 1, num_true_causal_snps[it] if it in num_true_causal_snps else 0))
            out_suma_fh.close()

        print('Number of SNPs simulated: ', np.sum([len(component_to_snp[comp]) for comp in component_to_cor]))
        locus_fh.close()
        gwas_fh.close()
