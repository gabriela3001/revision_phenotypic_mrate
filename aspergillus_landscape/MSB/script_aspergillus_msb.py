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
with open('aspergillus_paramgrid_msb.txt', 'rb') as f:
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
init_genotype_index = 248

all_results = []

for gamma2ind in range(len(switching_rates)):
    print(gamma2ind)
    gamma2 = switching_rates[gamma2ind]
    msb_pop = simulation_msb(all_genotypes_aspergillus, landscape_tuple_aspergillus, switching_rates[gamma1], gamma2, mu, tau, all_genotypes_aspergillus[init_genotype_index][1])['pop']

    all_results.append(msb_pop)
    
output = '/home/labs/pilpel/gabril/review_NIMR/range_gamma_0_1/MSB_aspergillus/'
    
with open(output+'results_'+str(param)+'.txt', 'wb') as f:
    pickle.dump(all_results, f)