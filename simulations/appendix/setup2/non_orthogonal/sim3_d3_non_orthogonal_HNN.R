
rm(list = ls())

seednum <- 2021 ; num_sim <- 130 ; setting <- 3

library(doMC) ; num_cores <- detectCores() ; registerDoMC(cores = num_cores)

source("./functions/gen_data.R")
source("./functions/all_source.R")

n <- 100 ; p1 <- p2 <- p3 <- 100 ; pvec <- c(p1, p2, p3)

if(setting==2){
  r0 <- 2 ; r12 <- 2 ; r13 <- 2 ; r23 <- 2 ; r1 <- 2 ; r2 <- 2 ; r3 <- 2 ; snr <- 2.0
} else if(setting==3){
  r0 <- 1 ; r12 <- 1 ; r13 <- 3 ; r23 <- 5 ; r1 <- 2 ; r2 <- 2 ; r3 <- 2 ; snr <- 2.0
} 


gam <- 0.45 ; nseq <- 10 ; qprob <- 1.0 ; max_iter <- 5000 ; eps <- 10^-5 
nrow_fold <- ncol_fold <- 2 ; tol <- 10^-3 
smin <- 1.0 ; smax <- 1.5 

X_list <- cfold_list <- list() 
cfold_vec <- chosen_1se_vec <- chosen_min_vec <- chosen_best_vec <- c()

p.ind <- list()
p.ind[[1]] <- c(1:pvec[1])
p.ind[[2]] <- c((pvec[1]+1):(pvec[1]+pvec[2]))
p.ind[[3]] <- c((pvec[1]+pvec[2]+1):(pvec[1]+pvec[2]+pvec[3]))


hnn_res_array <- array(NA, dim = c(num_sim, 1 + 1 + length(pvec), 3))

dimnames(hnn_res_array) <- list(paste("replication", sep = " ", c(1:num_sim)), 
                                c("sfb", "rankM", "rankM1", "rankM2", "rankM3"), 
                                c("hnn_1se", "hnn_min", "hnn_best"))


##########################################
##### Repeat independent simulations #####
##########################################


for(i in 1:num_sim){
  
  
  set.seed(seednum + i)
  
  #########################
  #### Data generation ####
  #########################
  
  sjoint <- runif(r0, smin, smax) 
  s12 <- runif(r12, smin, smax) ; s13 <- runif(r13, smin, smax) ; s23 <- runif(r23, smin, smax)
  s1ind <- runif(r1, smin, smax) ; s2ind <- runif(r2, smin, smax) ; s3ind <- runif(r3, smin, smax)
  
  dat <- dgen_d3_ver2(n, pvec, sjoint, s12, s13, s23, s1ind, s2ind, s3ind, rep(snr, length(pvec)), scaleU = T, orthogonalV = T)
  
  
  true_M <- dat$true_M 
  
  
  # generate BCV folds 
  rfold <- sample(rep(seq_len(nrow_fold), length.out = n))
  
  for(d in 1:length(pvec)){
    cfold_list[[d]] <- sample(rep(seq_len(ncol_fold), length.out = pvec[d]))
    X_list[[d]] <- dat$X[,p.ind[[d]]]
    if(d==1){
      cfold_vec <- cfold_list[[d]] 
    } else{
      cfold_vec <- c(cfold_vec, cfold_list[[d]])
    }
  }
  
  # standardize the simulated data
  std <- stdX(X_list, center = F) ; Z_list <- std$Z_list 
  sim_fb_vec <- std$fb_vec # squared frobenius norms after column-centering
  sim_mean_mat <- cbind(std$mean_mat[[1]], std$mean_mat[[2]], std$mean_mat[[3]]) # matrix of column means
  
  # set up the tuning grid
  grid_obj <- set_grid_d3(Z_list, tol, nseq, qprob)
  df_tune <- grid_obj$df_tune ; nlamb <- nrow(df_tune)
  ratio_pair <- grid_obj$ratio_pair ; ratio_indiv <- grid_obj$ratio_indiv
  
  
  ################################
  ##### Solution path of HNN  ####
  ################################
  
  
  sol_path <- foreach(l = 1:nlamb) %dopar% {
    
    fit <- tryCatch(algo1_d3(Z_list, df_tune[l,1], df_tune[l,2], df_tune[l,3], 
                             ratio_pair, ratio_indiv, sim_fb_vec, sim_mean_mat, gam, max_iter, eps, tol), error = function(e) NA)
    
    if(sum(is.na(fit))==0){
      hnn_res <- summary_fn(fit$Mhat-sim_mean_mat, true_M, pvec, tol)
    } else{
      hnn_res <- list("sfb" = NA, "ranks" = NA)
    }
    
    hnn_res
    
  }
  
  hnn_rank_vec <- hnn_sfb_vec <- rep(NA, nlamb)
  
  for(l in 1:nlamb){
    hnn_rank_vec[l] <- sol_path[[l]]$ranks[1] ; hnn_sfb_vec[l] <- sum(sol_path[[l]]$sfb)
  }
  
  
  ##############################
  ###### BCV to tune HNN #######
  ##############################
  
  
  bcv_err_array <- array(NA, dim = c(nlamb, nrow_fold*ncol_fold, length(pvec)))
  bcv_dat <- make_bcv_split(X_list, rfold, cfold_list)
  
  
  for(b in 1:(nrow_fold*ncol_fold)){
    
    X_jk_list <- bcv_dat$X_jk[[b]] ; X_njk_list <- bcv_dat$X_njk[[b]] 
    X_jnk_list <- bcv_dat$X_jnk[[b]] ; X_njnk_list <- bcv_dat$X_njnk[[b]]
    
    bobj <- foreach(l = 1:nlamb) %dopar% {
      tryCatch(bcv_single_tune_d3(X_njnk_list, X_jnk_list, X_jk_list, X_njk_list,
                                  df_tune[l,1], df_tune[l,2], df_tune[l,3],
                                  ratio_pair, ratio_indiv, gam, max_iter, eps, tol), error = function(e) NA)}
    
    for(l in 1:nlamb){
      bcv_err_array[l,b,] <- bobj[[l]]
    }
    
  }
  
  
  ############################################
  ### Choose the optimal tuning parameters ###
  ############################################
  
  
  sel <- choose_opt(bcv_err_array, hnn_rank_vec)
  chosen_1se <- chosen_1se_vec[i] <- sel$se 
  chosen_min <- chosen_min_vec[i] <- sel$min 
  chosen_best <- chosen_best_vec[i] <- which.min(hnn_sfb_vec)
  
  hnn_res_array[i,1,"hnn_1se"] <- hnn_sfb_vec[chosen_1se] 
  hnn_res_array[i,1,"hnn_min"] <- hnn_sfb_vec[chosen_min] 
  hnn_res_array[i,1,"hnn_best"] <- hnn_sfb_vec[chosen_best] 
  
  hnn_res_array[i,c(2:(2+length(pvec))),"hnn_1se"] <- sol_path[[chosen_1se]]$ranks 
  hnn_res_array[i,c(2:(2+length(pvec))),"hnn_min"] <- sol_path[[chosen_min]]$ranks 
  hnn_res_array[i,c(2:(2+length(pvec))),"hnn_best"] <- sol_path[[chosen_best]]$ranks 
  
  print(paste(paste("Replication", sep = " ", i), sep = " ", "is done at", Sys.time()))
  
  
}


save(hnn_res_array, chosen_1se_vec, chosen_best_vec, chosen_min_vec, file = paste("./sim", sep = "", setting, "_d3_non_orthogonal_HNN.RData"))

