# Conserved-current loop extrusion (CCLE)
This repository contains all the necessary codes for running the CCLE simulation, carrying out experiment-simulation Hi-C comparison, and calculating various metrics for goodness-of-fit.

CCLE_main:
This folder comprehensively contains the codes for the CCLE simulation program. 
The only input file is the file of ChIP-seq data of the protein of interest. The most important output of the program is the simulated Hi-C map.
Please refer to the comments in the "Main.m" file for more details about parameter adjustment or explanation.
This folder also contains subfolders with the yeast ChIP-seq data used in the paper.

HiC_compare_exp_sim:
This folder contains the codes that generate Hi-C map comparison between the experimental and the simulated Hi-C. 
The main code to run is named "HiC_compare_exp_sim.m".
Please read the comments in the main code file for more instructions.

HiC_compare_exp_exp:
This folder contains the codes that generate Hi-C map comparison between two experimental Hi-C maps. 
The main code to run is named "HiC_compare_exp_exp.m".
Please read the comments in the main code file for more instructions.

ChIP-seq_correlation:
This folder contains the program that generates the ChIP-seq correlation plots. 
Please read to the comments in the file "chipseq_correlation" for more instructions.

Hi-C expeiments:
This folder contains a few examples of the experimental Hi-C maps used in the paper. 
These example experimental Hi-C can be used as tests or templates.
