from complex_landscape_msb import *
import os
import pickle
import matplotlib.pyplot as plt
import numpy as np
import scipy.stats
import itertools
from scipy.spatial.distance import hamming
from itertools import product
import time
from collections import Counter
import sys

with open('aspergillus_landscape.txt', 'rb') as f:
    aspergillus_landscape = pickle.load(f)
with open('paramgrid_complex_adaptation.txt', 'rb') as f:
    param_grid = pickle.load(f)
    
def return_all_genotypes(landscape):
    
    landscape_tuple = {}

    for g in landscape:
        landscape_tuple[tuple([str(gi) for gi in g])] = landscape[g]
        
    genotypes = [tuple(x) for x in landscape_tuple.keys()]

    all_genotypes = []
    for i in range(2):
        for g in genotypes:
            all_genotypes.append((i, g))
            
    return(landscape_tuple, all_genotypes)

landscape_tuple_aspergillus, all_genotypes_aspergillus = return_all_genotypes(aspergillus_landscape)
switching_rates = np.logspace(-6,np.log10(0.95),100)

param = int(sys.argv[1])

mu = param_grid[param]['U']
tau = param_grid[param]['tau']
gamma1 = param_grid[param]['gamma1']
gamma2 = param_grid[param]['gamma2']

init_genotype_index = 248

all_results = []

msb_pop = simulation_msb(all_genotypes_aspergillus, landscape_tuple_aspergillus, switching_rates[gamma1], switching_rates[gamma2], mu, tau, all_genotypes_aspergillus[init_genotype_index][1])['pop']
transition_matrix = construct_transition_matrix(all_genotypes_aspergillus, landscape_tuple_aspergillus,switching_rates[gamma1], switching_rates[gamma2], param_grid[param]['U'], param_grid[param]['tau'])
fitness_vector = construct_selection_matrix(all_genotypes_aspergillus, landscape_tuple_aspergillus)

nreps = 1000
pop_size = 1000
ncat = len(all_genotypes_aspergillus)
adaptation_results = np.zeros((nreps, ncat))

for rep in range(nreps):
    adaptation_results[rep] = np.random.multinomial(pop_size, msb_pop)
adaptation_results /= np.sum(adaptation_results, axis = 1, keepdims = 1)

t = 0
time_start = time.process_time()

ngen = 5000

#evol_mpf = np.zeros((ngen, nreps))
evol_argmax = np.zeros((ngen, nreps))
#evol_mrate = np.zeros((ngen, nreps))

while t < ngen:

    adaptation_results = adaptation_results @ transition_matrix
    adaptation_results *= fitness_vector
    adaptation_results /= np.sum(adaptation_results, axis = 1, keepdims = True)

    for rep in np.arange(nreps):
        adaptation_results[rep] = np.random.multinomial(pop_size, adaptation_results[rep])
        #evol_mpf[t, rep] = np.dot(adaptation_results[rep], fitness_vector) / pop_size
        evol_argmax[t, rep] = np.argmax(adaptation_results[rep][:255]+adaptation_results[rep][255:])
        #evol_mrate[t, rep] = np.sum(adaptation_results[rep][255:])
    t += 1


elapsed_time = time.process_time() - time_start
print(elapsed_time)
evol_fittest = [Counter(evol_argmax[i])[0] for i in range(5000)]

output = '/home/labs/pilpel/gabril/review_NIMR/range_gamma_0_1/complex_adaptation_aspergillus_0_1/'
with open(output+'results_'+str(param)+'.txt', 'wb') as f:
    pickle.dump(evol_fittest, f)