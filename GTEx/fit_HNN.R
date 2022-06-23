
rm(list = ls())

library(doMC) ; num_cores <- detectCores() ; registerDoMC(cores = num_cores) ; num_bcv_rep <- 10

source("./functions/all_source.R")
load("./GTEx/GTEx_all.RData")

gam <- 0.45 ; nseq <- 10 ; qprob <- 1.0 ; max_iter <- 5000 ; eps <- 10^-5 ; tol <- 10^-3 

# standardize the simulated data
std <- stdX(X_list, center = F) ; Z_list <- std$Z_list 
fb_vec <- std$fb_vec # squared frobenius norms after column-centering
mean_mat <- cbind(std$mean_mat[[1]], std$mean_mat[[2]], std$mean_mat[[3]]) # matrix of column means
X_centered <- cbind(X_list[[1]], X_list[[2]], X_list[[3]]) - mean_mat
pvec <- c(ncol(X_list[[1]]), ncol(X_list[[2]]), ncol(X_list[[3]]))

# set up the tuning grid
grid_obj <- set_grid_d3(Z_list, tol, nseq, qprob)
df_tune <- grid_obj$df_tune ; nlamb <- nrow(df_tune)
ratio_pair <- grid_obj$ratio_pair ; ratio_indiv <- grid_obj$ratio_indiv


################################
##### Solution path of HNN  ####
################################


sol_path <- foreach(l = 1:nlamb) %dopar% {
  
  fit <- tryCatch(algo1_d3(Z_list, df_tune[l,1], df_tune[l,2], df_tune[l,3], 
                           ratio_pair, ratio_indiv, fb_vec, mean_mat, gam, max_iter, eps, tol), error = function(e) NA)
  
  if(sum(is.na(fit))==0){
    
    Mhat_decentered <- fit$Mhat-mean_mat
    percents<- c(round(100*sum(Mhat_decentered**2)/sum(X_centered**2), 2),
                     round(100*sum(Mhat_decentered[,c(1:pvec[1])]**2)/sum(X_centered[,c(1:pvec[1])]**2), 2),
                     round(100*sum(Mhat_decentered[,c((pvec[1]+1):(pvec[1]+pvec[2]))]**2)/sum(X_centered[,c((pvec[1]+1):(pvec[1]+pvec[2]))]**2), 2),
                     round(100*sum(Mhat_decentered[,-c(1:(pvec[1]+pvec[2]))]**2)/sum(X_centered[,-c(1:(pvec[1]+pvec[2]))]**2), 2))
    
    sgvs <- sgv_all(fit$Mhat-mean_mat, pvec, tol) 
    ranks <- sapply(c(1:length(sgvs)), FUN = function(i){length(sgvs[[i]])})
    
    hnn_res <- list("percents" =  percents, "ranks" = ranks)
    
  } else{
    
    hnn_res <- list("percents" = NA, "ranks" = NA)
    
  }
  
  hnn_res
  
}


bcv_err_array <- array(NA, dim = c(nlamb, nrow_fold*ncol_fold, length(pvec), num_bcv_rep))
cfold_list <- list()

############################################
######### Repeat BCV independently #########
############################################


for(i in 1:num_bcv_rep){
  
  set.seed(seednum + i)
  
  # generate BCV folds 
  rfold <- sample(rep(seq_len(nrow_fold), length.out = n))
  
  for(d in 1:length(pvec)){
    cfold_list[[d]] <- sample(rep(seq_len(ncol_fold), length.out = pvec[d]))
  }
  
  bcv_dat <- make_bcv_split(Z_list, rfold, cfold_list)
  
  for(b in 1:(nrow_fold*ncol_fold)){
    
    X_jk_list <- bcv_dat$X_jk[[b]] ; X_njk_list <- bcv_dat$X_njk[[b]] 
    X_jnk_list <- bcv_dat$X_jnk[[b]] ; X_njnk_list <- bcv_dat$X_njnk[[b]]
    
    bobj <- foreach(l = 1:nlamb) %dopar% {
      tryCatch(bcv_single_tune_d3(X_njnk_list, X_jnk_list, X_jk_list, X_njk_list,
                                  df_tune[l,1], df_tune[l,2], df_tune[l,3],
                                  ratio_pair, ratio_indiv, gam, max_iter, eps, tol), error = function(e) NA)}
    
    for(l in 1:nlamb){
      bcv_err_array[l,b,,i] <- bobj[[l]]
    }
    
  }
  
  print(paste(paste("BCV split", sep = " ", i), sep = " ", "is done at", Sys.time()))
  
}



chosen_1se_vec <- chosen_min_vec <- c()

for(i in 1:num_bcv_rep){
  sel <- choose_opt(bcv_err_array[,,,i], hnn_rank_vec)
  chosen_1se_vec[i] <- sel$se 
  chosen_min_vec[i] <- sel$min 
}


table(chosen_1se_vec) # Candidate : 705th or 706th in df_tune

sol_path[[705]]$ranks[1] # rank (M) = 54
sol_path[[706]]$ranks[1] # rank (M) = 49  

# choose 706th since it gives more lower concatenated rank than 705th

# Calculate the final HNN estimator
l <- 706
fit <- tryCatch(algo1_d3(Z_list, df_tune[l,1], df_tune[l,2], df_tune[l,3], 
                         ratio_pair, ratio_indiv, fb_vec, mean_mat, gam, max_iter, eps, tol), error = function(e) NA)

hnn_Mhat <- fit$Mhat

save(hnn_Mhat, file = "./GTEx/fit_HNN.RData")

