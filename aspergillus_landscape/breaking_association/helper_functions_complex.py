import numpy as np

def mutator_single_mut(genotypes_names, pop_vec, initial_genotype):
    p = 0
    for g in range(len(genotypes_names)):
        if (genotypes_names[g][1] != initial_genotype) and genotypes_names[g][0] == 1:
            p += pop_vec[g]
    return(p)

def baseline_single_mut(genotypes_names, pop_vec, initial_genotype):
    p = 0
    for g in range(len(genotypes_names)):
        if (genotypes_names[g][1] != initial_genotype) and genotypes_names[g][0] == 0:
            p += pop_vec[g]
    return(p)

def baseline_nomut(genotypes_names, pop_vec, initial_genotype):
    p = 0
    for g in range(len(genotypes_names)):
        if (genotypes_names[g][1] == initial_genotype) and genotypes_names[g][0] == 0:
            p += pop_vec[g]
    return(p)

def mutator_nomut(genotypes_names, pop_vec, initial_genotype):
    p = 0
    for g in range(len(genotypes_names)):
        if (genotypes_names[g][1] == initial_genotype) and genotypes_names[g][0] == 1:
            p += pop_vec[g]
    return(p)

def pS_proportion(genotypes_names, pop_vec, initial_genotype):
    return(mutator_single_mut(genotypes_names, pop_vec, initial_genotype) + baseline_single_mut(genotypes_names, pop_vec, initial_genotype))

def pM_proportion(genotypes_names, pop_vec, initial_genotype):
    return(mutator_single_mut(genotypes_names, pop_vec, initial_genotype) + mutator_nomut(genotypes_names, pop_vec, initial_genotype))

def pR_proportion(genotypes_names, pop_vec, initial_genotype):
    return(baseline_single_mut(genotypes_names, pop_vec, initial_genotype)/baseline_nomut(genotypes_names, pop_vec, initial_genotype))

def population_limit_pM(genotypes_names, old_pop, init_genotype_index, new_pM):
    
    true_pM = pM_proportion(genotypes_names, old_pop, genotypes_names[init_genotype_index][1])
    true_pS = pS_proportion(genotypes_names, old_pop, genotypes_names[init_genotype_index][1])
    true_M1 = mutator_single_mut(genotypes_names, old_pop, genotypes_names[init_genotype_index][1])
    true_m1 = baseline_single_mut(genotypes_names, old_pop, genotypes_names[init_genotype_index][1])
    true_association = true_M1 - true_pM*true_pS
    
    new_M1 = np.max([true_association+new_pM*true_pS,0])
    new_m1 = np.max([true_pS-new_M1, 0])
    new_M0 = np.max([new_pM-new_M1,0])
    new_m0 = [1-new_M1-new_m1-new_M0 if 1-new_M1-new_m1-new_M0>0 else 0][0]
    
    #print(new_M1, new_m1, new_M0, new_m0)

    new_population = np.zeros(len(genotypes_names))

    genotypes_M1 = [g for g in range(len(genotypes_names)) if genotypes_names[g][0] == 1 and genotypes_names[g][1] != genotypes_names[init_genotype_index][1]]
    genotypes_m1 = [g for g in range(len(genotypes_names)) if genotypes_names[g][0] == 0 and genotypes_names[g][1] != genotypes_names[init_genotype_index][1]]
    genotypes_M0 = [g for g in range(len(genotypes_names)) if genotypes_names[g][0] == 1 and genotypes_names[g][1] == genotypes_names[init_genotype_index][1]]
    genotypes_m0 = [g for g in range(len(genotypes_names)) if genotypes_names[g][0] == 0 and genotypes_names[g][1] == genotypes_names[init_genotype_index][1]]

    new_population[genotypes_M1] = new_M1*(old_pop[genotypes_M1]/np.sum(old_pop[genotypes_M1]))
    new_population[genotypes_m1] = new_m1*(old_pop[genotypes_m1]/np.sum(old_pop[genotypes_m1]))
    new_population[genotypes_M0] = new_M0*(old_pop[genotypes_M0]/np.sum(old_pop[genotypes_M0]))
    new_population[genotypes_m0] = new_m0*(old_pop[genotypes_m0]/np.sum(old_pop[genotypes_m0]))

    return(new_population/np.sum(new_population))

def population_breaking_association(genotypes_names, old_pop, init_genotype_index):
    
    true_pM = pM_proportion(genotypes_names, old_pop, genotypes_names[init_genotype_index][1])
    true_pS = pS_proportion(genotypes_names, old_pop, genotypes_names[init_genotype_index][1])
    true_M1 = mutator_single_mut(genotypes_names, old_pop, genotypes_names[init_genotype_index][1])
    true_m1 = baseline_single_mut(genotypes_names, old_pop, genotypes_names[init_genotype_index][1])
    true_association = true_M1 - true_pM*true_pS
    
    new_M1 = np.max([true_pM*true_pS,0])
    new_m1 = np.max([(1-true_pM*true_pS),0])
    new_M0 = np.max([true_pM*(1-true_pS),0])
    new_m0 = np.max([(1-true_pM)*(1-true_pS),0])
    
    #print(new_M1, new_m1, new_M0, new_m0)

    new_population = np.zeros(len(genotypes_names))

    genotypes_M1 = [g for g in range(len(genotypes_names)) if genotypes_names[g][0] == 1 and genotypes_names[g][1] != genotypes_names[init_genotype_index][1]]
    genotypes_m1 = [g for g in range(len(genotypes_names)) if genotypes_names[g][0] == 0 and genotypes_names[g][1] != genotypes_names[init_genotype_index][1]]
    genotypes_M0 = [g for g in range(len(genotypes_names)) if genotypes_names[g][0] == 1 and genotypes_names[g][1] == genotypes_names[init_genotype_index][1]]
    genotypes_m0 = [g for g in range(len(genotypes_names)) if genotypes_names[g][0] == 0 and genotypes_names[g][1] == genotypes_names[init_genotype_index][1]]

    new_population[genotypes_M1] = new_M1*(old_pop[genotypes_M1]/np.sum(old_pop[genotypes_M1]))
    new_population[genotypes_m1] = new_m1*(old_pop[genotypes_m1]/np.sum(old_pop[genotypes_m1]))
    new_population[genotypes_M0] = new_M0*(old_pop[genotypes_M0]/np.sum(old_pop[genotypes_M0]))
    new_population[genotypes_m0] = new_m0*(old_pop[genotypes_m0]/np.sum(old_pop[genotypes_m0]))

    return(new_population/np.sum(new_population))
