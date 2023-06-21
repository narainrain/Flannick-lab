"""
Simulation Models
====================================
This module implements different models to generate causal snps. The different models are null, infinitesimal,
non-infinitesimal, non-infinitesimal with snp clustering, and non-infinitesimal with snp clustering with causal snps
outside the casual gene region (PEGS)
"""

# Author: Sumit Narain <narain@broadinstitute.org>
#
# License:


class GenerateCausalGenesOrSnps:

    def __init__(self, options):
        self.options = options

    # method used to check if two ranges overlap
    def overlap_with_window(self, current_start, current_end, stored_start, stored_end, increase_window):
        return max(current_start - increase_window, stored_start) < min(current_end + increase_window, stored_end)

    # ----------------------------------------------------------------------------------------------------------------------
    # This is for the null case with no causal genes or snps => true betas should be 0
    def null_case(self):
        if (
                self.options.p_snps_in_causal_gene == 0 or self.options.p_causal_gene == 0) and self.options.p_causal_snps == 0:
            print("Entering Null case: therefore 0 causal snps")

            file_output = open(str(self.options.output_file) + ".causal_snps", "w")
            file_output.write("Chromosome\tSNP\tUnknown\tPosition\tRef\tAlt\n")
            file_output.close()

            file_output = open(str(self.options.output_file) + ".causal_genes", "w")
            file_output.write("Gene\tChromosome\tPosition Start\tPosition End\n")
            file_output.close()

        else:

            raise Exception("Parameter Error: Do not have the correct parameters to run null cases")

    # ----------------------------------------------------------------------------------------------------------------------
    # This is for the PEGs model where there are causal genes and causal snps inside, but there is a chance to generate causal
    # snps outside the causal genes as well
    def pegs(self):
        if self.options.p_causal_gene != 0 and self.options.p_snps_in_causal_gene != 0 and self.options.p_snps_out_causal_gene != 0:
            print(
                "Entered PEGS: Generation of snps within causal genes and generation of causal snps outside the causal gene")

            import random

            window = int(self.options.window)
            proportion_causal = self.options.p_causal_gene
            number_snps_causal_gene = 0
            gene_fh = open(self.options.gene_loc)
            all_genes = []

            # To determine how many genes there are and generate proportion of causal genes
            for line in gene_fh:
                all_genes.append(line)

            random.shuffle(all_genes)
            portion = len(all_genes) * proportion_causal
            causal_genes = all_genes[0:int(portion - 1)]
            chromosome_to_causal_gene_to_position = {}

            # Going through all the causal genes in the chromosome
            # If the window overlaps with another gene then ignore that portion
            # At the end of this will have control_mapping which will have chromosome -> gene -> position ranges
            # This will be used later when parsing through the bim file to see which variants will be sorted and keep count
            # Of all the causal variants, this will just provide us ranges with windows

            file_output = open(str(self.options.output_file) + ".causal_genes", "w")
            file_output.write("Gene\tChromosome\tPosition_Start\tPosition_End\n")

            for line in causal_genes:
                cols = line.strip().split()
                gene = cols[0]
                chromosome = cols[1]
                current_gene_start = int(cols[2])
                current_gene_end = int(cols[3])

                if chromosome in chromosome_to_causal_gene_to_position:

                    overlap = False

                    # Checking the current gene position versus the all stored genes' positions
                    for stored_gene in chromosome_to_causal_gene_to_position[chromosome]:

                        stored_gene_start = chromosome_to_causal_gene_to_position[chromosome][stored_gene][0]
                        stored_gene_end = chromosome_to_causal_gene_to_position[chromosome][stored_gene][1]

                        overlap = self.overlap_with_window(current_gene_start, current_gene_end, stored_gene_start,
                                                           stored_gene_end, window)
                        if overlap:
                            break

                    if overlap:
                        chromosome_to_causal_gene_to_position[chromosome][gene] = [current_gene_start, current_gene_end]
                        file_output.write(
                            '{}\t{}\t{}\t{}\n'.format(gene, chromosome, current_gene_start, current_gene_end))
                    else:
                        start_with_window = current_gene_start - window
                        end_with_window = current_gene_end + window
                        chromosome_to_causal_gene_to_position[chromosome][gene] = [start_with_window, end_with_window]
                        file_output.write(
                            '{}\t{}\t{}\t{}\n'.format(gene, chromosome, start_with_window, end_with_window))

                # If the chromosome is not currently in the array, add a new chromosome and add the gene as well with window
                else:
                    start_with_window = current_gene_start - window
                    end_with_window = current_gene_end + window
                    chromosome_to_causal_gene_to_position[chromosome] = {}
                    chromosome_to_causal_gene_to_position[chromosome][gene] = [start_with_window, end_with_window]
                    file_output.write('{}\t{}\t{}\t{}\n'.format(gene, chromosome, start_with_window, end_with_window))

            file_output.close()

            # Store causal snps to a file
            file_output = open(str(self.options.output_file) + ".causal_snps", "w")
            file_output.write("Chromosome\tSNP\tPosition\n")

            # Use the BIM file to randomly choose causal variants
            if True:
                number_snps_outside_causal_gene = 0
                for i in range(self.options.chrom_start, 23):
                    file_location = self.options.raw_file_location

                    ld_file = file_location + str(i) + ".ld"
                    ld_fh = open(ld_file)

                    # Skip the header when processing variants
                    ld_fh.readline()

                    # Set should not store duplicates
                    variant_already_processed = set()

                    # Store all variants in array from ld_file to make sure it's also in the bim file
                    # If not in bim file then ignore it
                    for line in ld_fh:
                        cols = line.strip().split()
                        chromosome_1 = cols[0]
                        position_1 = cols[1]
                        variant_1 = cols[2]
                        chromosome_2 = cols[3]
                        position_2 = cols[4]
                        variant_2 = cols[5]

                        variant_holder = [[chromosome_1, position_1, variant_1], [chromosome_2, position_2, variant_2]]

                        for variant_info in variant_holder:

                            chromosome = variant_info[0]
                            position = int(variant_info[1])
                            variant = variant_info[2]
                            line = chromosome + "\t" + variant + "\t" + str(position) + "\n"

                            variant_in_causal_gene = False

                            if chromosome in chromosome_to_causal_gene_to_position:
                                if variant not in variant_already_processed:
                                    for gene in chromosome_to_causal_gene_to_position[chromosome]:
                                        if chromosome_to_causal_gene_to_position[chromosome][gene][0] <= \
                                                position <= chromosome_to_causal_gene_to_position[chromosome][gene][1]:

                                            variant_in_causal_gene = True
                                            variant_already_processed.add(variant)

                                            if random.random() <= self.options.p_snps_in_causal_gene:
                                                number_snps_causal_gene = number_snps_causal_gene + 1
                                                file_output.write(line)

                            # Check if the variant is in the causal gene and if not then run the other probability
                            if not variant_in_causal_gene and variant not in variant_already_processed:

                                variant_already_processed.add(variant)

                                if random.random() <= self.options.p_snps_out_causal_gene:
                                    number_snps_outside_causal_gene = number_snps_outside_causal_gene + 1
                                    number_snps_causal_gene = number_snps_causal_gene + 1
                                    file_output.write(line)

                    ld_fh.close()

            file_output.close()
            print("Number of total causal snps: ", number_snps_causal_gene)
            print("Number of causal snps outside causal gene: ", number_snps_outside_causal_gene)

        else:
            raise Exception("Parameter Error: Do not have the correct parameters to run PEGS")

    # ----------------------------------------------------------------------------------------------------------------------
    # This method will only create causal snps from randomized causal genes but ignoring generation of causal snps outside the causal genes
    def non_infinitesimal_snp_cluster(self):
        print("This is the snps in causal gene: ", self.options.p_snps_in_causal_gene)
        if self.options.p_causal_gene != 0 and self.options.p_snps_in_causal_gene != 0:
            print("Entered Non-Infinitesimal Snp Cluster: generation of snps within causal genes only")

            import random

            window = int(self.options.window)
            proportion_causal = self.options.p_causal_gene
            number_snps_causal_gene = 0
            gene_fh = open(self.options.gene_loc)
            all_genes = []

            # To determine how many genes there are and generate proportion of causal genes
            for line in gene_fh:
                all_genes.append(line)

            random.shuffle(all_genes)
            portion = len(all_genes) * proportion_causal
            causal_genes = all_genes[0:int(portion - 1)]
            chromosome_to_gene_to_position = {}

            file_output = open(str(self.options.output_file) + ".causal_genes", "w")
            file_output.write("Gene\tChromosome\tPosition Start\tPosition End\n")

            for line in causal_genes:
                cols = line.strip().split()
                gene = cols[0]
                chromosome = cols[1]
                current_gene_start = int(cols[2])
                current_gene_end = int(cols[3])

                if chromosome in chromosome_to_gene_to_position:

                    overlap = False

                    # Checking the current gene position versus the all stored genes' positions
                    for stored_gene in chromosome_to_gene_to_position[chromosome]:

                        stored_gene_start = chromosome_to_gene_to_position[chromosome][stored_gene][0]
                        stored_gene_end = chromosome_to_gene_to_position[chromosome][stored_gene][1]

                        overlap = self.overlap_with_window(current_gene_start, current_gene_end, stored_gene_start,
                                                           stored_gene_end, window)

                        if overlap:
                            break

                    if overlap:
                        chromosome_to_gene_to_position[chromosome][gene] = [current_gene_start, current_gene_end]
                        file_output.write(
                            '{}\t{}\t{}\t{}\n'.format(gene, chromosome, current_gene_start, current_gene_end))
                    else:
                        start_with_window = current_gene_start - window
                        end_with_window = current_gene_end + window
                        chromosome_to_gene_to_position[chromosome][gene] = [start_with_window, end_with_window]
                        file_output.write(
                            '{}\t{}\t{}\t{}\n'.format(gene, chromosome, start_with_window, end_with_window))

                # If the chromosome is not currently in the array, add a new chromosome and add the gene as well with window
                else:
                    start_with_window = current_gene_start - window
                    end_with_window = current_gene_end + window
                    chromosome_to_gene_to_position[chromosome] = {}
                    chromosome_to_gene_to_position[chromosome][gene] = [start_with_window, end_with_window]
                    file_output.write('{}\t{}\t{}\t{}\n'.format(gene, chromosome, start_with_window, end_with_window))

            file_output.close()

            # Store causal snps to a file
            file_output = open(str(self.options.output_file) + ".causal_snps", "w")
            file_output.write("Chromosome\tSNP\tPosition\n")

            if True:
                for i in range(self.options.chrom_start, 23):
                    file_location = self.options.raw_file_location

                    ld_file = file_location + str(i) + ".ld"
                    ld_fh = open(ld_file)

                    # Skip the header when processing variants
                    ld_fh.readline()

                    # Set should not store duplicates
                    variant_already_processed = set()

                    # Store all variants in array from ld_file to make sure it's also in the bim file
                    # If not in bim file then ignore it
                    for line in ld_fh:
                        cols = line.strip().split()
                        chromosome_1 = cols[0]
                        position_1 = cols[1]
                        variant_1 = cols[2]
                        chromosome_2 = cols[3]
                        position_2 = cols[4]
                        variant_2 = cols[5]

                        variant_holder = [[chromosome_1, position_1, variant_1], [chromosome_2, position_2, variant_2]]

                        for variant_info in variant_holder:

                            chromosome = variant_info[0]
                            position = int(variant_info[1])
                            variant = variant_info[2]
                            line = chromosome + "\t" + variant + "\t" + str(position) + "\n"

                            if chromosome in chromosome_to_gene_to_position:
                                if variant not in variant_already_processed:
                                    for gene in chromosome_to_gene_to_position[chromosome]:
                                        if chromosome_to_gene_to_position[chromosome][gene][0] <= \
                                                position <= chromosome_to_gene_to_position[chromosome][gene][1]:

                                            variant_already_processed.add(variant)

                                            if random.random() <= self.options.p_snps_in_causal_gene:
                                                number_snps_causal_gene = number_snps_causal_gene + 1
                                                file_output.write(line)

                    ld_fh.close()

            print("Number of causal snps: ", number_snps_causal_gene)

        else:
            raise Exception("Parameter Error: Incorrect Parameters for Non-Infinitesimal Snp Cluster")

    # ----------------------------------------------------------------------------------------------------------------------
    # Non-infinitesimal Model
    # If we just need to generate the causal snps for the entire genome without worrying about the gene or positions
    # This is used for the non-infinitesimal case where causal snps are randomized throughout genome

    def non_infinitesimal(self):

        if self.options.p_causal_snps != 0:

            import random

            print("Entered Non-Infinitesimal Model: generation of causal snps without genes with probability: ",
                  self.options.p_causal_snps)
            proportion_causal = self.options.p_causal_snps

            file_output = open(str(self.options.output_file) + ".causal_genes", "w")
            file_output.write("Gene\tChromosome\tPosition Start\tPosition End\n")
            file_output.close()

            file_output = open(str(self.options.output_file) + ".causal_snps", "w")
            file_output.write("Chromosome\tSNP\tPosition\n")

            # Parse through LD file to generate causal snps
            if True:
                for i in range(self.options.chrom_start, 23):
                    print("In chromosome", i)
                    count_causal = 0
                    file_location = self.options.raw_file_location

                    ld_file = file_location + str(i) + ".ld"
                    ld_fh = open(ld_file)

                    # Ignore Header
                    ld_fh.readline()

                    # Set should not store duplicates
                    variant_already_processed = set()

                    # Store all variants in array from ld_file to make sure it's also in the bim file
                    # If not in bim file then ignore it
                    # Working under the assumption that variants should not cross over to different chromosomes
                    # Otherwise I would store all the variants in another loop

                    for line in ld_fh:
                        cols = line.strip().split()
                        chromosome_1 = cols[0]
                        position_1 = cols[1]
                        variant_1 = cols[2]
                        chromosome_2 = cols[3]
                        position_2 = cols[4]
                        variant_2 = cols[5]

                        variant_holder = [[chromosome_1, position_1, variant_1], [chromosome_2, position_2, variant_2]]

                        for variant_info in variant_holder:

                            chromosome = variant_info[0]
                            position = int(variant_info[1])
                            variant = variant_info[2]
                            line = chromosome + "\t" + variant + "\t" + str(position) + "\n"

                            if variant not in variant_already_processed:
                                variant_already_processed.add(variant)

                                if random.random() <= proportion_causal:
                                    count_causal = count_causal + 1
                                    file_output.write(line)

                    ld_fh.close()
            file_output.close()

    # ----------------------------------------------------------------------------------------------------------------------
    # infinitesimal Model
    # If we just need to generate the causal snps for the entire genome without worrying about the gene or positions
    # This is used for the infinitesimal case where every snp is a causal snps

    def infinitesimal(self):

        if self.options.p_causal_snps == 1:

            import random

            print("Entered Non-Infinitesimal Model: generation of causal snps without genes with probability: ",
                  self.options.p_causal_snps)
            proportion_causal = self.options.p_causal_snps

            file_output = open(str(self.options.output_file) + ".causal_genes", "w")
            file_output.write("Gene\tChromosome\tPosition Start\tPosition End\n")
            file_output.close()

            file_output = open(str(self.options.output_file) + ".causal_snps", "w")
            file_output.write("Chromosome\tSNP\tPosition\n")

            # Parse through LD file to generate causal snps
            if True:
                for i in range(self.options.chrom_start, 23):
                    print("In chromosome", i)
                    count_causal = 0
                    file_location = self.options.raw_file_location

                    ld_file = file_location + str(i) + ".ld"
                    ld_fh = open(ld_file)

                    # Ignore Header
                    ld_fh.readline()

                    # Set should not store duplicates
                    variant_already_processed = set()

                    # Store all variants in array from ld_file to make sure it's also in the bim file
                    # If not in bim file then ignore it
                    # Working under the assumption that variants should not cross over to different chromosomes
                    # Otherwise I would store all the variants in another loop

                    for line in ld_fh:
                        cols = line.strip().split()
                        chromosome_1 = cols[0]
                        position_1 = cols[1]
                        variant_1 = cols[2]
                        chromosome_2 = cols[3]
                        position_2 = cols[4]
                        variant_2 = cols[5]

                        variant_holder = [[chromosome_1, position_1, variant_1],
                                          [chromosome_2, position_2, variant_2]]

                        for variant_info in variant_holder:

                            chromosome = variant_info[0]
                            position = int(variant_info[1])
                            variant = variant_info[2]
                            line = chromosome + "\t" + variant + "\t" + str(position) + "\n"

                            if variant not in variant_already_processed:
                                variant_already_processed.add(variant)

                                if random.random() <= proportion_causal:
                                    count_causal = count_causal + 1
                                    file_output.write(line)

                    ld_fh.close()
            file_output.close()

    # ----------------------------------------------------------------------------------------------------------------------
    # This function is used when you want to compare the non-infinitesimal mode with snp-clustering against non-infinitesimal
    #
    # Initially all causal snps within the causal gene will be generated
    # and then the proportion will be generated given this output
    # This will calculate the proportion of causal snps / total snps = self.options.p_causal_snps
    # This value will be used to store causal snps from both methods, so it is comparable in future analysis by
    # having the same probability value. Also ignoring the p_snps outside the causal genes

    def non_infinitesimal_snp_cluster_vs_non_infinitesimal(self):

        if self.options.p_causal_gene != 0 and self.options.p_snps_in_causal_gene != 0:
            print("Entering Comparison between Non-Infinitesimal Snp Cluster versus Non-Infinitesimal")

            import random

            window = int(self.options.window)
            proportion_causal = self.options.p_causal_gene
            number_snps_in_causal_gene = 0
            gene_fh = open(self.options.gene_loc)
            all_genes = []
            # To determine how many genes there are and generate proportion of causal genes
            for line in gene_fh:
                all_genes.append(line)

            gene_fh.close()

            random.shuffle(all_genes)
            portion = len(all_genes) * proportion_causal
            causal_genes = all_genes[0:int(portion - 1)]
            chromosome_to_gene_to_position = {}

            number_snps_in_ld_file = 0

            # Going through all the genes in the chromosome
            # If the window overlaps with another gene then ignore that portion
            # At the end of this will have control_mapping which will have chromosome -> gene -> position ranges
            # This will be used later when parsing through the bim file to see which variants will be sorted and keep count
            # Of all the causal variants, this will just provide us ranges with windows

            file_output = open(str(self.options.output_file) + ".causal_genes", "w")
            file_output.write("Gene\tChromosome\tPosition Start\tPosition End\n")

            for line in causal_genes:
                cols = line.strip().split()
                gene = cols[0]
                chromosome = cols[1]
                current_gene_start = int(cols[2])
                current_gene_end = int(cols[3])

                if chromosome in chromosome_to_gene_to_position:

                    overlap = False

                    # Checking the current gene position versus the all stored genes' positions
                    for stored_gene in chromosome_to_gene_to_position[chromosome]:

                        stored_gene_start = chromosome_to_gene_to_position[chromosome][stored_gene][0]
                        stored_gene_end = chromosome_to_gene_to_position[chromosome][stored_gene][1]

                        overlap = self.overlap_with_window(current_gene_start, current_gene_end, stored_gene_start,
                                                           stored_gene_end,
                                                           window)

                        if overlap:
                            break

                    if overlap:
                        chromosome_to_gene_to_position[chromosome][gene] = [current_gene_start, current_gene_end]
                        file_output.write(
                            '{}\t{}\t{}\t{}\n'.format(gene, chromosome, current_gene_start, current_gene_end))
                    else:
                        start_with_window = current_gene_start - window
                        end_with_window = current_gene_end + window
                        chromosome_to_gene_to_position[chromosome][gene] = [start_with_window, end_with_window]
                        file_output.write(
                            '{}\t{}\t{}\t{}\n'.format(gene, chromosome, start_with_window, end_with_window))

                # If the chromosome is not currently in the array, add a new chromosome and add the gene as well with window
                else:
                    start_with_window = current_gene_start - window
                    end_with_window = current_gene_end + window
                    chromosome_to_gene_to_position[chromosome] = {}
                    chromosome_to_gene_to_position[chromosome][gene] = [start_with_window, end_with_window]
                    file_output.write('{}\t{}\t{}\t{}\n'.format(gene, chromosome, start_with_window, end_with_window))

            file_output.close()

            if True:
                for i in range(self.options.chrom_start, 23):
                    file_location = self.options.raw_file_location

                    ld_file = file_location + str(i) + ".ld"
                    ld_fh = open(ld_file)

                    # Skip the header when processing variants
                    ld_fh.readline()

                    # Set should not store duplicates
                    variant_already_processed = set()

                    # Store all variants in array from ld_file to make sure it's also in the bim file
                    # If not in bim file then ignore it
                    for line in ld_fh:
                        cols = line.strip().split()
                        chromosome_1 = cols[0]
                        position_1 = cols[1]
                        variant_1 = cols[2]
                        chromosome_2 = cols[3]
                        position_2 = cols[4]
                        variant_2 = cols[5]

                        variant_holder = [[chromosome_1, position_1, variant_1], [chromosome_2, position_2, variant_2]]

                        for variant_info in variant_holder:

                            chromosome = variant_info[0]
                            position = int(variant_info[1])
                            variant = variant_info[2]

                            if chromosome in chromosome_to_gene_to_position:
                                if variant not in variant_already_processed:
                                    number_snps_in_ld_file = number_snps_in_ld_file + 1
                                    variant_already_processed.add(variant)

                                    for gene in chromosome_to_gene_to_position[chromosome]:
                                        if chromosome_to_gene_to_position[chromosome][gene][0] <= \
                                                position <= chromosome_to_gene_to_position[chromosome][gene][1]:

                                            if random.random() <= self.options.p_snps_in_causal_gene:
                                                number_snps_in_causal_gene = number_snps_in_causal_gene + 1

                    ld_fh.close()

            # Using the proportion generated from the Non-infinitesimal SNP cluster we will generate
            # SNPS from the Non-infinitesimal model and store those snps
            prop_causal = number_snps_in_causal_gene / number_snps_in_ld_file

            # Now add the causal snps
            print("This is the proportion: ", prop_causal)

            file_output = open(str(self.options.output_file) + ".causal_snps", "w")
            file_output.write("Chromosome\tSNP\tPosition\n")

            count_causal = 0
            # Parse through LD file to generate causal snps
            if True:
                for i in range(self.options.chrom_start, 23):
                    print("In chromosome", i)

                    file_location = self.options.raw_file_location

                    ld_file = file_location + str(i) + ".ld"
                    ld_fh = open(ld_file)

                    # Ignore Header
                    ld_fh.readline()

                    # Set should not store duplicates
                    variant_already_processed = set()

                    # Store all variants in array from ld_file to make sure it's also in the bim file
                    # If not in bim file then ignore it
                    # Working under the assumption that variants should not cross over to different chromosomes
                    # Otherwise I would store all the variants in another loop

                    for line in ld_fh:
                        cols = line.strip().split()
                        chromosome_1 = cols[0]
                        position_1 = cols[1]
                        variant_1 = cols[2]
                        chromosome_2 = cols[3]
                        position_2 = cols[4]
                        variant_2 = cols[5]

                        variant_holder = [[chromosome_1, position_1, variant_1], [chromosome_2, position_2, variant_2]]

                        for variant_info in variant_holder:

                            chromosome = variant_info[0]
                            position = int(variant_info[1])
                            variant = variant_info[2]
                            line = chromosome + "\t" + variant + "\t" + str(position) + "\n"

                            if variant not in variant_already_processed:
                                variant_already_processed.add(variant)

                                if random.random() <= prop_causal:
                                    count_causal = count_causal + 1
                                    file_output.write(line)

                    ld_fh.close()
            file_output.close()

            # Print out all the variant snps, Causal Genes, Number of variant snps
            print("Number of causal snps: ", count_causal)

        else:
            raise Exception(
                "Parameter Error: Need the correct parameters to run Non-Infinitesimal Snp Cluster versus Non-Infinitesimal")
