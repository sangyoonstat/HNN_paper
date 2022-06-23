
rm(list = ls())

seednum <- 2021 ; num_sim <- 120 ; jive_orth <- F

library(r.jive) ; library(doMC) ; num_cores <- detectCores() ; registerDoMC(cores = num_cores)

source("./functions/gen_data.R")
source("./functions/all_source.R")
source("./functions/my_slide.R")
source("./functions/my_unifacs.R")


n <- 150 ; p1 <- 50 ; p2 <- 50 ; pvec <- c(p1, p2)
r0 <- 2 ; theta <- c(30, 40, 50, 60) ; snr <- 1

true_jive <- true_ajive <- c()
true_jive[1] <- r0 ; true_jive[2] <- true_jive[3] <- length(theta) 

p.ind <- list()
p.ind[[1]] <- c(1:pvec[1])
p.ind[[2]] <- c((pvec[1]+1):(pvec[1]+pvec[2]))

p.ind.list <- list()
p.ind.list[[1]] <- unlist(p.ind[1:2])
p.ind.list[[2]] <- unlist(p.ind[c(1)])
p.ind.list[[3]] <- unlist(p.ind[c(2)])

nrow_fold <- ncol_fold <- 2 ; tol <- 10^-3 ; smin <- 1.0 ; smax <- 1.5 

X_list <- cfold_list <- list() 
cfold_vec <- chosen_1se_vec <- chosen_min_vec <- chosen_best_vec <- c()

####################################################################
##### Function to implement all competing methods except AJIVE #####
####################################################################

# jive_orth: indicates whether the argument ``orthIndiv" in jive would be set to be true or not, see r.jive
# true_jive: 1st element = true rank of joint space // others = true ranks of individuals

fit_others <- function(i, seednum, true_jive, jive_orth){
  
  set.seed(seednum + i)
  
  #########################
  #### Data generation ####
  #########################
  
  # 1.1 Generate each scores 
  pairs <- gen_two_scores(n, theta) ; allU <- gen_others_d2(r0, pairs)
  
  # 1.2 Generate data
  sjoint <- runif(r0, smin, smax) ; s1ind <- runif(length(theta), smin, smax) ; s2ind <- runif(length(theta), smin, smax)
  dat <- dgen_d2_givenU(allU, pvec, sjoint, s1ind, s2ind, rep(snr, length(pvec)), orthogonalV = T)  
  
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
  sim_mean_mat <- cbind(std$mean_mat[[1]], std$mean_mat[[2]]) # matrix of column means
  
  
  ###################################
  ###### SLIDE with 2 x 2 BCV #######
  ###################################
  
  out_slide = slide(cbind(Z_list[[1]], Z_list[[2]]), pvec = pvec, 
                    fold_id_n = rfold, fold_id_p = cfold_vec, center = F)
  Mtilde <- out_slide$model$U%*%t(out_slide$model$V) 
  slide_Mhat <- cbind(sim_fb_vec[1]*Mtilde[,p.ind[[1]]], sim_fb_vec[2]*Mtilde[,p.ind[[2]]])
  slide_res <- summary_fn(slide_Mhat, true_M, pvec, tol) ; rm(slide_Mhat)
  
  ####################
  ###### JIVE ########
  ####################
  
  Z_list2 <- list() 
  
  for(d in 1:length(pvec)){
    Z_list2[[d]] <- t(Z_list[[d]])
  }
  
  jive_out <- jive(Z_list2, method = "perm", scale = F, orthIndiv = jive_orth, est = T, showProgress = F)
  
  jive_Mhat <- cbind(sim_fb_vec[1]*t(jive_out$joint[[1]] + jive_out$individual[[1]]), 
                     sim_fb_vec[2]*t(jive_out$joint[[2]] + jive_out$individual[[2]]))
  
  jive_res <- summary_fn(jive_Mhat, true_M, pvec, tol) ; rm(jive_Mhat)
  
  ###################################
  ####### JIVE with true ranks ######
  ###################################
  
  jive_true_out <- jive(Z_list2, rankJ = true_jive[1], rankA = true_jive[-1], 
                        method = "given", scale = F, orthIndiv = jive_orth, est = T, showProgress = F)
  
  jive_true_Mhat <- cbind(sim_fb_vec[1]*t(jive_true_out$joint[[1]] + jive_true_out$individual[[1]]), 
                          sim_fb_vec[2]*t(jive_true_out$joint[[2]] + jive_true_out$individual[[2]]))
  
  jive_true_res <- summary_fn(jive_true_Mhat, true_M, pvec, tol) ; rm(jive_true_Mhat)
  
  #######################
  ###### UNIFAC #########
  #######################
  
  unifac_fit <- my_UNIFAC(X_list, center = F, seednum = seednum + i)
  unifac_Mhat <- cbind(unifac_fit$S[[1]], unifac_fit$S[[2]])
  unifac_res <- summary_fn(unifac_Mhat, true_M, pvec, tol) ; rm(unifac_Mhat) 
  
  
  ########################
  ###### UNIFAC+ #########
  ########################
  
  X00 <- unifac_fit$X00 ; sigma_vec <- unifac_fit$sigma_vec
  
  ###### res.g : UNIFAC+ given structure
  
  res.g <- unifac.plus.given(X0 = t(X00), p.ind = p.ind, p.ind.list = p.ind.list)
  
  for(k in 1:length(res.g$S)){
    if(k==1){
      unifac_plus_given_mat <- res.g$S[[1]]
    } else{
      unifac_plus_given_mat <- unifac_plus_given_mat + res.g$S[[k]]
    }
  }
  
  unifac_plus_given_Mhat <- cbind(sigma_vec[1]*t(unifac_plus_given_mat)[,p.ind[[1]]],
                                  sigma_vec[2]*t(unifac_plus_given_mat)[,p.ind[[2]]])
  
  unifac_plus_res <- summary_fn(unifac_plus_given_Mhat, true_M, pvec, tol)
  
  ####################################
  ### refitting UNIFAC and UNIFAC+ ###
  ####################################
  
  Z_list2 <- list()
  Z_list2[[1]] <- unifac_fit$X00[,p.ind[[1]]] ; Z_list2[[2]] <- unifac_fit$X00[,p.ind[[2]]]
  
  unifac_Mtilde <- cbind(unifac_fit$S[[1]]/sigma_vec[1], unifac_fit$S[[2]]/sigma_vec[2])
  unifac_refit <- refit_indiv(Z_list2, unifac_Mtilde, pvec, tol) 
  unifac_refitted_Mhat <- cbind(unifac_fit$sigma_vec[1]*unifac_refit[[1]], unifac_fit$sigma_vec[2]*unifac_refit[[2]])
  unifac_refit_res <- summary_fn(unifac_refitted_Mhat, true_M, pvec, tol)
  
  unifac_plus_refit <- refit_indiv(Z_list2, t(unifac_plus_given_mat), pvec, tol) 
  unifac_plus_refitted_Mhat <- cbind(sigma_vec[1]*unifac_plus_refit[[1]], 
                                     sigma_vec[2]*unifac_plus_refit[[2]])
  unifac_plus_refit_res <- summary_fn(unifac_plus_refitted_Mhat, true_M, pvec, tol)
  
  print(paste(paste("Replication", sep = " ", i), sep = " ", "is done at", Sys.time()))
  
  
  return(list("slide" = slide_res, "jive" = jive_res, "jive_true" = jive_true_res, 
              "unifac" = unifac_res, "unifac_refit" = unifac_refit_res, 
              "unifac_plus" = unifac_plus_res, "unifac_plus_refit" = unifac_plus_refit_res))
  
}


others_res <- foreach(i = 1:num_sim) %dopar% {
  tryCatch(fit_others(i, seednum, true_jive, jive_orth), error = function(e) NA)}



save(others_res, file = "./sim_d2_non_orthogonal_others.RData")



