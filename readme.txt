ReadMe File for the code used to generate results described in the paper: "Phenotype switching of the mutation rate facilitates adaptive evolution".

Each folder contains a wrapper, a main script, and other files needed to run the main script, such as parameter grids, fitness values for fitness landscapes, and helper functions. 

The folder "heatmap_asymmetric_switchingrate" contains code used to generate results for the simple model. The main script computes the mutation-selection balance frequencies of each possible pheno-genotype of the simple model. 
To run this script for many parameter sets, use the wrapper which passes the index of each parameter set in the parameter grid to the script and submits it to the computing cluster. 
Results of this script were used to plot Figure 2, Figure 3, Figure S2-S5, Figure S9, Figure S13-S14, and Figure S16.

The folder "gradual_change" contains code used to generate results for the model extension of gradual change in mutation rate. Use the wrapper to run the main script for all parameter sets. Results of this script were used to plot Figures S10-S11.

The folder "different_U" contains code used to compare results of the simple model for extreme values of beta. Use the wrapper to run the main script for all parameter sets. beta_1 is the script for beta = 1, and beta_1M is the script for beta = 1,000,000. Results of this script were used to generate Figure S12.

The folder "landscape_motifs" contains code used to study the rates of adaptation on different landscape motifs. Scripts for motifs of length 3 are in the folder "three_loci" and for motifs of length 4 are in the folder "four_loci". Results of this script were used to plot Figure 3, and Figure S17 and S18.

The folders "NK_landscape" and "aspergillus_landscape" contains code that was used to generate all results for complex adaptation. Each folder contains subfolders for MSB results, complex adaptation, as well as breaking association, eliminating hitchhiking, and limiting pM. This code was used to generated results plotted in Figure 5, and Figures S21-S29.