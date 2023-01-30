from model import *
from msb import *
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats
import pickle
import time as tm
import sys
from sklearn.model_selection import ParameterGrid

param = int(sys.argv[1])

with open('paramgrid_adaptation_time_fourloci_oneparamset.txt', 'rb') as f:
    param_grid = pickle.load(f)
    
noise_levels = np.logspace(-6, np.log10(0.95), 100)

U = param_grid[param]['U']
tau = param_grid[param]['tau']
s = param_grid[param]['s']
pop_size = param_grid[param]['pop_size']
landscape_motifs =  param_grid[param]['landscape_motifs']
switching_rate_1 =  param_grid[param]['gamma1']
switching_rate_2 =  param_grid[param]['gamma2']

genotype_names_3 = []
for i in range(4):
    for j in range(2):
        for k in range(5):
            genotype_names_3.append((i,j,k))
            
results_msb = []
noise1 = noise_levels[switching_rate_1]
noise2 = noise_levels[switching_rate_2]
start_time = tm.time()
msb_pop_array = msb_simulation_4(genotype_names_3, 'MSB', 
                                             s, -1, 
                                             U, 
                                             tau, 
                                             noise1, 
                                             noise2)['pop']
print(tm.time()-start_time)

adapted_genotypes = np.where(np.array(np.array(genotype_names_3)[:,0])==3)[0]

nreps = 1000
ncat = len(msb_pop_array)

strain_fitness_vect = landscape_motifs

transition_matrix = mutation_rate_transition_4(U, tau, 5, noise_levels[switching_rate_1], noise_levels[switching_rate_2], beta = 5000)
fitness_vector = fitness(genotype_names_3, param_grid[param]['s'], 5, 'V', strain_fitness_vect)
msb_pop = msb_pop_array

results = np.zeros((nreps, ncat))

for rep in range(nreps):
    results[rep] = np.random.multinomial(pop_size, msb_pop)

results /= np.sum(results, axis = 1, keepdims = 1)

# we want: array of times of first appearance; fate of each mutation
t = 0
update_fates = np.array([False]*nreps)
adaptation_evol = []

start_time = tm.time()

adapted_genotypes = np.where(np.array(np.array(genotype_names_3)[:,0])==3)[0]

while t < 100000:

    t += 1

    results = results @ transition_matrix
    results *= fitness_vector

    results /= np.sum(results, axis = 1, keepdims = True)


    for rep in np.arange(nreps)[~update_fates]:

        results[rep] = np.random.multinomial(pop_size, results[rep])

        if np.sum(results[rep][adapted_genotypes])/pop_size > 0.99:
            update_fates[rep] = True

    if t%1000 == 0:
        adaptation_evol.append(sum(update_fates))
        print(t)

print(tm.time() - start_time)

print('Finished Noise ' + str(switching_rate_2) + '  ' + str(adaptation_evol))

output = '/home/labs/pilpel/gabril/review_NIMR/range_gamma_0_1/motifs_oneparamset_four_loci_0_1/'


with open(output + 'results_' + str(param) + '.txt', 'wb') as f:
    pickle.dump({'U':U, 'tau':tau, 's':s, 'pop_size':pop_size, 'landscape_motifs':landscape_motifs, 'switchingrate1':switching_rate_1,
                 'switchingrate2': switching_rate_2,'results': adaptation_evol}, f)