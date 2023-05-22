import matplotlib.pyplot as plt
import numpy as np

from models import PEGS
import pegs_utils as pg
import time

# Instance creation
test_options = pg.arg_settings()
my_model = PEGS(test_options)
# Recording time
num_chains = 20
nelder_time = np.zeros(num_chains)
BFGS_time = np.zeros(num_chains)
for i in range(num_chains):
    start_time = time.time()
    my_model.gibbs(burn_in=5, gibbs_iter=10, plot=False, plot_window_regression=False, regression_method="nelder-mead")
    nelder_time[i] = time.time() - start_time
for i in range(num_chains):
    start_time = time.time()
    my_model.gibbs(burn_in=5, gibbs_iter=10, plot=False, plot_window_regression=False, regression_method="BFGS")
    BFGS_time[i] = time.time() - start_time

# plotting data
data = [nelder_time, BFGS_time]

fig = plt.figure(figsize =(10, 7))
ax = fig.add_subplot(111)
bp = ax.boxplot(data, patch_artist=True, notch='False', vert=0)
ax.set_yticklabels(['nelder_time', 'BFGS_time'])
plt.title("Execution time of optimization methods")
plt.show()
