"""
Example Run
====================================
This is an example file used to demonstrate how to run the simulations
First you need to determine which model will be used in the generateCausalGenesOrSnps file
After the model is determined and ran that will generate a file of causal genes and causal snps
The causal snps generated will be the input to running the class - Simulation
The simulation will only need to have one method - run() which will run the entire process of
turning the causal snps into simulated GWAS data

"""

# Author: Sumit Narain <narain@broadinstitute.org>
#
# License:

from generate_causal_genes_or_snps_1 import GenerateCausalGenesOrSnps
from simulation_pipeline import Simulation
import simulation_args as sim

# Instance creation
default_options = sim.arg_generate_causal_snps()
sim_default_options = sim.arg_simulations_pipeline()

my_model = GenerateCausalGenesOrSnps(default_options)
run_model = Simulation(sim_default_options)

my_model.infinitesimal()
run_model.run()


