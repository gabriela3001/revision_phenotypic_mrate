from complex_landscape_msb import *
from helper_functions_complex import *
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

with open('NK_landscapes_6.txt', 'rb') as f:
    NK_landscapes = pickle.load(f)
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

switching_rates = np.logspace(-6,np.log10(0.95),100)
param, NKind = int(sys.argv[1]), int(sys.argv[2])

mu = param_grid[param]['U']
tau = param_grid[param]['tau']
gamma1 = param_grid[param]['gamma1']
gamma2 = param_grid[param]['gamma2']

init_genotype_index = [11,24,57,5,17,49][NKind]
landscape_tuple_NK, all_genotypes_NK = return_all_genotypes(NK_landscapes[NKind])
all_results = []

NK_pop_sizes = [1000,1000,10000000,10000000,10000000,10000000]
NK_fittest_genotypes = [52,39,6,58,46,14]

nreps = 1000
pop_size = NK_pop_sizes[NKind]
ncat = len(all_genotypes_NK)
adaptation_results = np.zeros((nreps, ncat))


msb_pop = simulation_msb(all_genotypes_NK, landscape_tuple_NK, switching_rates[gamma1], switching_rates[gamma2], mu, tau, all_genotypes_NK[init_genotype_index][1])['pop']
transition_matrix = construct_transition_matrix(all_genotypes_NK, landscape_tuple_NK,switching_rates[gamma1], switching_rates[gamma2], param_grid[param]['U'], param_grid[param]['tau'])
fitness_vector = construct_selection_matrix(all_genotypes_NK, landscape_tuple_NK)

for rep in range(nreps):
    adaptation_results[rep] = np.random.multinomial(pop_size, msb_pop)
adaptation_results /= np.sum(adaptation_results, axis = 1, keepdims = 1)

t = 0
time_start = time.process_time()

ngen = 600


evol_mpf = np.zeros((ngen, nreps))
evol_argmax = np.zeros((ngen, nreps))
evol_mrate = np.zeros((ngen, nreps))

while t < ngen:

    adaptation_results = adaptation_results @ transition_matrix
    adaptation_results *= fitness_vector
    adaptation_results /= np.sum(adaptation_results, axis = 1, keepdims = True)

    for rep in np.arange(nreps):
        wildtype = np.argmax(adaptation_results[rep][:64]+adaptation_results[rep][64:])
        adaptation_results[rep] = population_breaking_association(all_genotypes_NK, adaptation_results[rep], wildtype)
        adaptation_results[rep] = np.random.multinomial(pop_size, adaptation_results[rep])
        evol_mpf[t, rep] = np.dot(adaptation_results[rep], fitness_vector) / pop_size
        evol_argmax[t, rep] = np.argmax(adaptation_results[rep][:64]+adaptation_results[rep][64:])
        evol_mrate[t, rep] = np.sum(adaptation_results[rep][64:])
    t += 1

    if t % 100 == 0:
        print(t)

elapsed_time = time.process_time() - time_start
evol_fittest = [Counter(evol_argmax[i])[NK_fittest_genotypes[NKind]] for i in range(600)]

output = '/home/labs/pilpel/gabril/review_NIMR/range_gamma_0_1/breaking_association_NK/'

with open(output+'results_'+str(param)+'_'+str(NKind)+'.txt', 'wb') as f:
    pickle.dump(evol_fittest, f)