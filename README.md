# HNN_paper

Codes used for numerical studies and real data analysis in the paper "Hierarchical nuclear norm penalization for multi-view data integration" by Sangyoon Yi, Raymond K.W. Wong and Irina Gaynanova.

## 1. File description

### 1.1 functions

**DBFB.rcpp** - Rcpp implementation of dual block-coordinate forward-backward algorithm for hierarchical nuclear norm (HNN) penalization 
	
**all_source.R** - R code for our main method and several auxiliary functions (constructing grid, refitting and modified BCV) 

**gen_data.R** - data generation used in simulations

**myslide.R**, **myunifacs.R** - wrapper for running JIVE, BIDIFAC and BIDIFAC+ in R 

### 1.2 simulations

Inside this folder, there are "main" and "appendix" folders having functions to run simulation stuides in the main paper and appendix. In each folder, the R files ending with HNN and others shows codes for HNN and other methods except AJIVE. "dat" folder has simulated data sets we used for numerical summaires for reproducibility.

### 1.3 GTEx

**GTEx_all.RData, GTEx_data.mat** - Preprocessed Genotype-Tissue Expression (GTEx) data by [Li and Jung (2017)] (https://github.com/reagan0323/SIFA)

**fit_AJIVE.m, fit_BIDIFACs.R, fit_HNN.R, fit_JIVE.R, fit_SLIDE.R** - applying each method to the GTEx data

### 1.4 all_sim_AJIVE.m
Matlab code to reproduce simulation results of AJIVE. The code depends on using AJIVE as [implemented in MATLAB here](https://github.com/MeileiJiang/AJIVE_Project).


## 2. Example

