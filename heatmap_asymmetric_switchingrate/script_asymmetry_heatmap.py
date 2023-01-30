from model import *
from msb import *
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats
import pickle
import time as tm
import sys

param = int(sys.argv[1])
param_noise = int(sys.argv[2])

with open('paramgrid_asymmetry.txt', 'rb') as f: 
  param_grid = pickle.load(f)

noise_levels = np.logspace(-6, np.log10(0.95), 100)

s = param_grid[param]['s']
tau = param_grid[param]['tau']
mu = param_grid[param]['mu']

genotype_names_3 = []
for i in range(3):
    for j in range(2):
        for k in range(5):
            genotype_names_3.append((i,j,k))
            
results_msb = []
noise = noise_levels[param_noise]
start_time = tm.time()
msb_pop = [msb_simulation(genotype_names_3, 'MSB', 
                                             s, -1, 
                                             mu, 
                                             tau, 
                                             noise, 
                                             noise2)['pop'] for noise2 in noise_levels]
print(tm.time()-start_time)

output = '/home/labs/pilpel/gabril/review_NIMR/range_gamma_0_1/asymmetry_paramgrid_0_1/'

with open(output + 'results_' + str(param) + '_' + str(param_noise)+'.txt', 'wb') as f:
    pickle.dump(msb_pop, f)