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

with open('NK_landscapes_6.txt', 'rb') as f:
    NK_landscapes = pickle.load(f)
with open('NK_paramgrid_msb.txt', 'rb') as f:
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
init_genotype_index = [11,24,57,5,17,49][NKind]
landscape_tuple_NK, all_genotypes_NK = return_all_genotypes(NK_landscapes[NKind])
all_results = []

for gamma2ind in range(len(switching_rates)):
    print(gamma2ind)
    gamma2 = switching_rates[gamma2ind]
    msb_pop = simulation_msb(all_genotypes_NK, landscape_tuple_NK, switching_rates[gamma1], gamma2, mu, tau, all_genotypes_NK[init_genotype_index][1])['pop']
    all_results.append(msb_pop)
    
output = '/home/labs/pilpel/gabril/review_NIMR/range_gamma_0_1/MSB_NK/'
    
with open(output+'results_'+str(param)+'_'+str(NKind)+'.txt', 'wb') as f:
    pickle.dump(all_results, f)

