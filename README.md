# HNN_paper

Codes used for numerical studies and real data analysis in the paper "Hierarchical nuclear norm penalization for multi-view data integration" by Sangyoon Yi, Raymond K.W. Wong and Irina Gaynanova.

## 1. File description

### 1.1 "functions" folder

**DBFB.rcpp** - Rcpp implementation of dual block-coordinate forward-backward algorithm for hierarchical nuclear norm penalization 
	
**all_source.R** - R code for our main method and several auxiliary functions 

**gen_data.R** - data generation used in simulations

**myslide.R**, **myunifacs.R** - wrapper for running JIVE, BIDIFAC and BIDIFAC+ in R 

### 1.2 "simulations" folder

Inside "simulations" folder, there are "main" and "appendix" folders having functions to run simulations for HNN and other competing methods except AJIVE. In each folder, "dat" folder has simulated data sets for reproducibility.

## 2. Example

