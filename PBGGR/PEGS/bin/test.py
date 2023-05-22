from models import PEGS
import pegs_utils as pg
import time

# Instance creation
test_options = pg.arg_settings()
my_model = PEGS(test_options)

start_time = time.time()
my_model.gibbs_multi_chain(num_chains=20, burn_in=50, max_iter=1000,
                           plot_window_regression=False, infinite=False, local_test=True,
                           convergence_thr=1.0001)
print("--- %s seconds ---" % (time.time() - start_time))
my_model.plot_results()
# my_model.gibbs(burn_in=2, gibbs_iter=10, plot=True, plot_window_regression=False)
