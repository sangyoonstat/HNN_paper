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

**GTEx_all.RData, GTEx_data.mat** - Preprocessed Genotype-Tissue Expression (GTEx) data in [Li and Jung (2017)](https://github.com/reagan0323/SIFA)

**fit_AJIVE.m, fit_BIDIFACs.R, fit_HNN.R, fit_JIVE.R, fit_SLIDE.R** - applying each method to the GTEx data

### 1.4 all_sim_AJIVE.m
MATLAB code to reproduce simulation results of AJIVE using the one [implemented in here](https://github.com/MeileiJiang/AJIVE_Project).


## 2. Example

### 2.0 Get ready to run HNN

```{r}

##################################
#### Setup parallel computing #### 
##################################

require(doMC) ; num_cores <- detectCores() ; registerDoMC(cores = num_cores)

#############################
##### Load the codes ########
#############################

source("./functions/gen_data.R")
source("./functions/all_source.R")

### gam : step parameter for dual forward-backward algorithm
### nseq : number of sequence for each tuning parameter
### max_iter: maximum number of iteration for dual forward-backward algorithm
### eps: numerical tolerance for dual forward-backward algorithm
### nrow_fold, ncol_fold: the number of row and column fold for BCV
### tol: numerical tolerance for singular value

gam <- 0.45 ; nseq <- 20 ; max_iter <- 5000 ; eps <- 10^-5 
nrow_fold <- ncol_fold <- 2 ; tol <- 10^-3 

```


### 2.1 Generate data
```{r}

########################
###  Data generation ###
########################

### n: sample size
### p1, p2 : number of variables in 1st and 2nd view
### r0 : rank of joint structure
### r1, r2: rank of individual structure
### snr: signal to noise ratio
### smin, smax: lower and upper bounds used to generate singular values

n <- 150 ; p1 <- 50 ; p2 <- 50 ; pvec <- c(p1, p2) 
r0 <- r1 <- r2 <- 2 ; snr <- 1 ; smin <- 1.0 ; smax <- 1.5 

set.seed(2022)

sjoint <- runif(r0, smin, smax) ; s1ind <- runif(r1, smin, smax) ; s2ind <- runif(r2, smin, smax) # generate singular values
dat <- dgen_orth_d2(n, pvec, sjoint, s1ind, s2ind, rep(snr, length(pvec)), orthogonalV = T) # generate data (orthogonal case)

X_list <- cfold_list <- list() 
X_list[[1]] <- dat$X[,c(1:p1)] # X_1
X_list[[2]] <- dat$X[,c((p1+1):(p1+p2))] # X_2
true_M <- dat$true_M # true signal

# generate BCV folds 
rfold <- sample(rep(seq_len(nrow_fold), length.out = n))

for(d in 1:length(pvec)){
  cfold_list[[d]] <- sample(rep(seq_len(ncol_fold), length.out = pvec[d]))
}

```

### 2.2 Standardize data and construct gird

```{r}

# standardize the simulated data
std <- stdX(X_list, center = F) ; Z_list <- std$Z_list 
sim_fb_vec <- std$fb_vec # squared frobenius norms after column-centering
sim_mean_mat <- cbind(std$mean_mat[[1]], std$mean_mat[[2]]) # matrix of column means


# set up the tuning grid
grid_obj <- set_grid_d2(Z_list, tol, nseq, qprob)
df_tune <- grid_obj$df_tune ; nlamb <- nrow(df_tune)
ratio_indiv <- grid_obj$ratio_indiv

```


### 2.3 Calculate HNN estimator over tuning gird

```{r}
################################
##### Solution path of HNN  ####
################################

sol_path <- foreach(l = 1:nlamb) %dopar% {
  
  fit <- tryCatch(algo1_d2(Z_list, df_tune[l,1], df_tune[l,2], ratio_indiv, 
                           sim_fb_vec, sim_mean_mat, gam, max_iter, eps, tol), error = function(e) NA)
  
  hnn_res <- summary_fn(fit$Mhat-sim_mean_mat, true_M, pvec, tol)
  
  hnn_res
  
}


# Save the concatenated ranks and scaled Frobenius norm over solution path

hnn_rank_vec <- hnn_sfb_vec <- rep(NA, nlamb)

for(l in 1:nlamb){
  hnn_rank_vec[l] <- sol_path[[l]]$ranks[1] ; hnn_sfb_vec[l] <- sum(sol_path[[l]]$sfb)
}

```


### 2.4 Run modified BCV

```{r}


##############################
###### BCV to tune HNN #######
##############################

bcv_err_array <- array(NA, dim = c(nlamb, nrow_fold*ncol_fold, length(pvec)))
bcv_dat <- make_bcv_split(X_list, rfold, cfold_list)

for(b in 1:(nrow_fold*ncol_fold)){
  
  X_jk_list <- bcv_dat$X_jk[[b]] ; X_njk_list <- bcv_dat$X_njk[[b]] 
  X_jnk_list <- bcv_dat$X_jnk[[b]] ; X_njnk_list <- bcv_dat$X_njnk[[b]]
  
  bobj <- foreach(l = 1:nlamb) %dopar% {
    tryCatch(bcv_single_tune_d2(X_njnk_list, X_jnk_list, X_jk_list, X_njk_list, df_tune[l,1], df_tune[l,2], 
                                ratio_indiv, gam, max_iter, eps, tol), error = function(e) NA)}
  
  for(l in 1:nlamb){
    bcv_err_array[l,b,] <- bobj[[l]]
  }
  
}


```

### 2.5 Choosing tuning parameters

```{r}
# choose tuning parameter
sel <- choose_opt(bcv_err_array, hnn_rank_vec) 
df_tune[sel$se,] # one-SE rule
df_tune[which.min(hnn_sfb_vec),] # best in terms of Frobenius norm error

sol_path[[sel$se]] # resulting ranks and Frobenius norm error by one-SE rule
sol_path[[which.min(hnn_sfb_vec)]] # resulting ranks and Frobenius norm error by best 

```
