
rm(list = ls())

seednum <- 2021 ; num_sim <- 120 

require(doMC) ; num_cores <- detectCores() ; registerDoMC(cores = num_cores)

source("./functions/gen_data.R")
source("./functions/all_source.R")

n <- 150 ; p1 <- 50 ; p2 <- 50 ; pvec <- c(p1, p2)
r0 <- 2 ; theta <- c(30, 40, 50, 60) ; snr <- 1

gam <- 0.45 ; nseq <- 20 ; qprob <- 1.0 ; max_iter <- 5000 ; eps <- 10^-5 
nrow_fold <- ncol_fold <- 2 ; tol <- 10^-3 
smin <- 1.0 ; smax <- 1.5 

X_list <- cfold_list <- list() 
cfold_vec <- chosen_1se_vec <- chosen_min_vec <- chosen_best_vec <- c()

hnn_res_array <- array(NA, dim = c(num_sim, 1 + 1 + length(pvec), 3))

dimnames(hnn_res_array) <- list(paste("replication", sep = " ", c(1:num_sim)), 
                                c("sfb", "rankM", "rankM1", "rankM2"), 
                                c("hnn_1se", "hnn_min", "hnn_best"))


##########################################
##### Repeat independent simulations #####
##########################################


for(i in 1:num_sim){
  
  
  set.seed(seednum + i)
  
  ########################
  ## 1. Data generation ##
  ########################
  
  # 1.1 Generate each scores 
  pairs <- gen_two_scores(n, theta) ; allU <- gen_others_d2(r0, pairs)
  
  # 1.2 Generate data
  sjoint <- runif(r0, smin, smax) ; s1ind <- runif(length(theta), smin, smax) ; s2ind <- runif(length(theta), smin, smax)
  dat <- dgen_d2_givenU(allU, pvec, sjoint, s1ind, s2ind, rep(snr, length(pvec)), orthogonalV = T)  
  
  true_M <- dat$true_M 
  X_list[[1]] <- dat$X[,c(1:p1)] 
  X_list[[2]] <- dat$X[,c((p1+1):(p1+p2))] 

  # generate BCV folds 
  rfold <- sample(rep(seq_len(nrow_fold), length.out = n))
  
  for(d in 1:length(pvec)){
    cfold_list[[d]] <- sample(rep(seq_len(ncol_fold), length.out = pvec[d]))
    if(d==1){
      cfold_vec <- cfold_list[[d]] 
    } else{
      cfold_vec <- c(cfold_vec, cfold_list[[d]])
    }
  }
  
  # standardize the simulated data
  std <- stdX(X_list, center = F) ; Z_list <- std$Z_list 
  sim_fb_vec <- std$fb_vec # squared frobenius norms after column-centering
  sim_mean_mat <- cbind(std$mean_mat[[1]], std$mean_mat[[2]]) # matrix of column means
  
  # set up the tuning grid
  grid_obj <- set_grid_d2(Z_list, tol, nseq, qprob)
  df_tune <- grid_obj$df_tune ; nlamb <- nrow(df_tune)
  ratio_indiv <- grid_obj$ratio_indiv
  
  
  ################################
  ##### Solution path of HNN  ####
  ################################
  
  
  sol_path <- foreach(l = 1:nlamb) %dopar% {
    
    fit <- tryCatch(algo1_d2(Z_list, df_tune[l,1], df_tune[l,2], ratio_indiv, 
                             sim_fb_vec, sim_mean_mat, gam, max_iter, eps, tol), error = function(e) NA)
    
    if(sum(is.na(fit))==0){
      hnn_res <- summary_fn(fit$Mhat-sim_mean_mat, true_M, pvec, tol)
    } else{
      hnn_res <- list("sfb" = NA, "ranks" = NA)
    }
    
    hnn_res
    
  }
  
  # Save the concatenated ranks and scaled Frobenius norm over solution path
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
      tryCatch(bcv_single_tune_d2(X_njnk_list, X_jnk_list, X_jk_list, X_njk_list, df_tune[l,1], df_tune[l,2], 
                                  ratio_indiv, gam, max_iter, eps, tol), error = function(e) NA)}
    
    for(l in 1:nlamb){
      bcv_err_array[l,b,] <- bobj[[l]]
    }
    
  }
  
  
  ############################################
  ### Choose the optimal tuning parameters ###
  ############################################
  
  sel <- choose_opt(bcv_err_array, hnn_rank_vec) 
  chosen_1se <- chosen_1se_vec[i] <- sel$se # one-SE rule
  chosen_min <- chosen_min_vec[i] <- sel$min # minimum rule
  chosen_best <- chosen_best_vec[i] <- which.min(hnn_sfb_vec) # best in terms of Frobenius norm error
  
  hnn_res_array[i,1,"hnn_1se"] <- hnn_sfb_vec[chosen_1se] 
  hnn_res_array[i,1,"hnn_min"] <- hnn_sfb_vec[chosen_min] 
  hnn_res_array[i,1,"hnn_best"] <- hnn_sfb_vec[chosen_best] 
  
  hnn_res_array[i,c(2:(2+length(pvec))),"hnn_1se"] <- sol_path[[chosen_1se]]$ranks 
  hnn_res_array[i,c(2:(2+length(pvec))),"hnn_min"] <- sol_path[[chosen_min]]$ranks 
  hnn_res_array[i,c(2:(2+length(pvec))),"hnn_best"] <- sol_path[[chosen_best]]$ranks 
  
  print(paste(paste("Replication", sep = " ", i), sep = " ", "is done at", Sys.time()))
  
}


save(hnn_res_array, chosen_1se_vec, chosen_best_vec, chosen_min_vec, file = "./sim_d2_non_orthogonal_HNN.RData")



