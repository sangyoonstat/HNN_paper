
# the package "denoiseR" is to set up the tuning grid in the function "set_grid" below
require(denoiseR) 

# load some functions written in Rcpp
require(Rcpp)
require(RcppArmadillo)
sourceCpp("./functions/DBFB.cpp")

############################################
##### 1. function to standardize data ######
############################################

## Input
# X_list : data from each view in each element of X_list
# If center = F, do scaling only

## Output
# Z_list : list where each element is standardized data for each view
# mean_mat: list where each element is resulting matrix for each view after column centering
# fb_vec : vector of scaled Frobenius norms for scaling each view
# If center = F, each element of fb_vec is the squared Frobenius norm of the original data
# If center = T, each element of fb_vec is the squared Frobenius norm of the column-centered data

stdX <- function(X_list, center = F){
  
  D <- length(X_list) ; Z_list <- mean_mat <- list() ; fb_vec <- c() 
  
  if(center){
    for(d in 1:D){
      col_means <- colMeans(X_list[[d]])
      mean_mat[[d]] <- matrix(col_means, nrow = nrow(X_list[[d]]), ncol = ncol(X_list[[d]]), byrow = T)
      centered_X <- X_list[[d]] - mean_mat[[d]] 
      fb_vec[d] <- sqrt(sum(centered_X**2)) ; Z_list[[d]] <- centered_X/fb_vec[d] 
    }
  } else{
    for(d in 1:D){
      mean_mat[[d]] <- matrix(0, nrow = nrow(X_list[[d]]), ncol = ncol(X_list[[d]])) 
      fb_vec[d] <- sqrt(sum(X_list[[d]]**2)) ; Z_list[[d]] <- X_list[[d]]/fb_vec[d] 
    }
  }
  
  return(list(Z_list = Z_list, mean_mat = mean_mat, fb_vec = fb_vec))
  
}

############################################################
#### 2. Some useful functions to summarize the results #####
############################################################

#######################################################################################
########## 2.1 Calculate cosine of principal angles between two column spaces #########
#######################################################################################

cal_angles <- function(mat1, mat2, tol){
  lobj <- tSVD_left(cbind(mat1, mat2), tol) ; Uall <- lobj$u ; rtilde <- lobj$r
  if(rtilde>0){
    est_left1 <- tSVD_left(mat1, tol) ; est_left2 <- tSVD_left(mat2, tol)
    if(est_left1$r>0 & est_left2$r>0){
      angles <- sgv_single(crossprod(tSVD_left(tcrossprod(Uall) %*% est_left1$u, tol)$u, 
                                     tSVD_left(tcrossprod(Uall) %*% est_left2$u, tol)$u), tol)
      #angles <- sgv_single(crossprod(tcrossprod(Uall) %*% est_left1$u, tcrossprod(Uall) %*% est_left2$u), tol)
    } else{
      angles <- 0 
    }
  } else{
    angles <- 0
  }
  return(angles)
}


#####################################################################################
########## 2.2 calculate the ranks of each combinations for D = 2 and 3 #############
#####################################################################################

## 1. Input
# mat : input matrix
# pvec : vector of p_d (the number of variables for view d)
# tol : numerical cutoff for small singular values

## 2. Output
# If length(pvec) = 2 (two view case), show the singular values of the concatenated (denoted as "sgv")
# and those of view 1 (denoted as "sgv1") and view 2 (denoted as "sgv2")
# If length(pvec) = 3 (two view case), show the singular values of the concatenated (denoted as "sgv")
# and those of individual view (denoted as "sgv1", "sgv2", "sgv3")
# and those of pairwise view (denoted as "sgv12", "sgv13", "sgv23")

sgv_all <- function(mat, pvec, tol){
  
  if(length(pvec)==2){
    p1 <- pvec[1] ; sgv <- sgv_single(mat, tol) 
    sgv1 <- sgv_single(mat[,c(1:p1)], tol) ; sgv2 <- sgv_single(mat[,-c(1:p1)], tol) 
    return(list("sgv" = sgv, "sgv1" = sgv1, "sgv2" = sgv2))
  } else if(length(pvec)==3){
    p1 <- pvec[1] ; p2 <- pvec[2]
    sgv <- sgv_single(mat, tol) ; sgv1 <- sgv_single(mat[,c(1:p1)], tol)
    sgv2 <- sgv_single(mat[,c((p1+1):(p1+p2))], tol) ; sgv3 <- sgv_single(mat[,-c(1:(p1+p2))], tol)
    sgv12 <- sgv_single(mat[,c(1:(p1+p2))], tol) ; sgv13 <- sgv_single(mat[,-c((p1+1):(p1+p2))], tol)
    sgv23 <- sgv_single(mat[,-c(1:p1)], tol)
    return(list("sgv" = sgv, "sgv1" = sgv1, "sgv2" = sgv2, "sgv3" = sgv3))    
  }
  
}

####################################################################################
########## 2.3 Calculate the scaled Frobenius norm errors for D = 2 and 3 ##########
####################################################################################

## 1. Input
# Mhat : estimator for trueM
# trueM : true signal
# pvec : vector of p_d (the number of variables for view d)

## 2. Output
# If length(pvec) = 2 (two view case), 
# show the scaled Frobenius norm errors for view 1 (denoted as "fb1") and view 2 (denoted as "fb2")
# If length(pvec) = 3 (three view case), 
# show the scaled Frobenius norm errors for view 1 (denoted as "fb1") and view 2 (denoted as "fb2") and view 3 (denoted as "fb3")

cal_sfb <- function(Mhat, trueM, pvec){
  
  if(length(pvec)==2){
    p1 <- pvec[1] ; M1hat <- Mhat[,c(1:p1)] ; M2hat <- Mhat[,-c(1:p1)] 
    trueM1 <- trueM[,c(1:p1)] ; trueM2 <- trueM[,-c(1:p1)] 
    fb1 <- sum((M1hat-trueM1)**2)/sum(trueM1**2) ; fb2 <- sum((M2hat-trueM2)**2)/sum(trueM2**2) 
    return(c(fb1, fb2))      
  } else if(length(pvec)==3){
    p1 <- pvec[1] ; p2 <- pvec[2]
    M1hat <- Mhat[,c(1:p1)] ; M2hat <- Mhat[,c((p1+1):(p1+p2))] ; M3hat <- Mhat[,-c(1:(p1+p2))]
    trueM1 <- trueM[,c(1:p1)] ; trueM2 <- trueM[,c((p1+1):(p1+p2))] ; trueM3 <- trueM[,-c(1:(p1+p2))]
    fb1 <- sum((M1hat-trueM1)**2)/sum(trueM1**2) 
    fb2 <- sum((M2hat-trueM2)**2)/sum(trueM2**2) 
    fb3 <- sum((M3hat-trueM3)**2)/sum(trueM3**2)
    return(c(fb1, fb2, fb3))    
  }
  
}

########################################################
##### 2.4 function to calculate generalized inverse ####
########################################################

## 1. Input
# mat : input matrix
# tol : numerical cutoff for small singular values

## 2. Output
# resulting generalized inverse

cal_ginv <- function(mat, tol){
  sv <- full_tSVD(mat, tol) ; ind <- sv$d>tol ; r <- sum(ind)
  if(r>0){
    U <- sv$u[,1:r] ; V <- sv$v[,1:r]
    B <- diag(1/sv$d[ind], nrow = r, ncol = r)
    ginv_mat <- V%*%B%*%t(U)
  } else{
    ginv_mat <- matrix(0, nrow = ncol(mat), ncol = nrow(mat))
  }
  return(ginv_mat)
}

########################################################################
#####  3. Functions to find the tuning parameter by the SURE  ##########
########################################################################

###################################################################
### 3.1 Calculate the SURE formula in Cand{`e`}s et al. (2013) ####
###################################################################

cal_SURE <- function(n, p, sgv, lambda, noise_sd){
  
  RSS <- sum(sapply(c(1:length(sgv)), FUN = function(i){min(sgv[i]^2, lambda^2)}))
  ind <- which(sgv>lambda)
  v <- sapply(ind, FUN = function(i){sum(sgv[i]*(sgv[i]-lambda)/(sgv[i]^2-sgv[-i]^2))})  
  
  if(length(v)==0|length(unique(sgv))!=length(sgv)){
    div_est <- 0
  } else{
    div_est <- abs(n-p)*sum(1-lambda/sgv[ind]) + length(ind) + 2*sum(v)
  }
  sure_est <- -n*p*noise_sd^2+RSS+2*noise_sd^2*div_est
  
  return(sure_est)
}

#########################################################################################
#### 3.2 Find the optimal one which minimizes the SURE criteria based on grid search ####
#########################################################################################

find_SURE <- function(X, nlamb, tol){
  n <- nrow(X) ; p <- ncol(X) ; sgv <- sgv_single(X, tol)
  lam_seq <- seq(0, max(sgv), len = nlamb)
  shat <- estim_sigma(X, method = "MAD", center = FALSE)
  SURE_vec <- sapply(c(1:nlamb), FUN = function(l){cal_SURE(n, p, sgv, lam_seq[l], shat)})
  lambda <- lam_seq[which.min(SURE_vec)]
  return(lambda)
}

###################################################################################
########## 4. functions to set up the tuning grid using the ratios ################
###################################################################################

################################################
### 4.1 Make the triangle-type region in 2D  ###
################################################

tri_2D <- function(seq_x, seq_y, keep_close = F){
  
  len_x <- length(seq_x) ; max_x <- max(seq_x) ; max_y <- max(seq_y)
  
  # Calculate the line passing the 3 points, (max_x, 0) and (0, max_y)
  border_line <- function(x){
    -(max_y/max_x)*x + max_y
  }
  
  for(nx in 1:len_x){
    
    if(nx==1){
      df <- data.frame("x" = seq_x[nx], "y" = seq_y)
    } else if(nx==len_x){
      df <- rbind.data.frame(df, data.frame("x" = seq_x[nx], "y" = 0))      
    } else{
      bd <- border_line(seq_x[nx]) ; new_y <- seq_y[seq_y<=bd]  
      
      if(keep_close){
        new_y <- c(new_y, bd)        
      } else if((bd-max(new_y))>0.05){
        new_y <- c(new_y, bd) 
      }
      
      if(nx%%2==0){
        new_y <- sort(new_y, decreasing = T)  
      } 
      
      df <- rbind.data.frame(df, data.frame("x" = seq_x[nx], "y" = new_y))
      
    }
    
  }
  
  return(df)
  
}


################################################
### 4.2 Make the triangle-type region in 3D  ###
################################################

tri_3D <- function(seq_x, seq_y, seq_z, keep_close = F){
  
  max_x <- max(seq_x) ; max_y <- max(seq_y) ; max_z <- max(seq_z)
  
  # Calculate the plane passing the 3 points, (max_x, 0, 0), (0, max_y, 0) and (0, 0, max_z)
  border_plane <- function(x, y){
    (-max_z/max_x)*x + (-max_z/max_y)*y + max_z
  }
  
  df_xy <- tri_2D(seq_x, seq_y, keep_close)
  
  # j <- 41
  
  for(j in 1:nrow(df_xy)){
    
    bd <- border_plane(df_xy[j,1], df_xy[j,2]) ; new_z <- seq_z[seq_z<=abs(bd)]
    
    if(sum(seq_z<=abs(bd))>0){
      
      if(keep_close){
        new_z <- c(new_z, bd)        
      } else if((bd-max(new_z))>0.05){
        new_z <- c(new_z, bd) 
      }
      
      if(j==1){
        df <- data.frame("x" = df_xy[j,1], "y" = df_xy[j,2], "z" = new_z)
      } else{
        df <- rbind.data.frame(df, data.frame("x" = df_xy[j,1], "y" = df_xy[j,2], "z" = new_z))
      }
      
    }
    
  }
  
  return(df)
  
}

########################################################################################
########## 4.3 functions to produce the tuning grid when D = 2 or D = 3 ################
########################################################################################

## 1. Input
# Z_list : list of standardized data for each view
# tol : numerical cutoff for small singular values
# nseq : size of sequence for each tuning parameter
# qprob : percentage multiplied by max_tau, max_kappa and max_lambda (see below)

## 2. Output
# df_tune : data frame of the 2 tuning parameters when D = 2 (3 tuning parameters when D = 3) in the ratio formulation
# sure_indiv : vector of cutoffs of singular value soft-thresholding chosen by SURE for each individual
# sure_pair : vector of cutoffs of singular value soft-thresholding chosen by SURE for each pair
# ratio_indiv : vector of ratios of individual cutoff by SURE to the sum of all individual cutoffs by SURE
# ratio_pair : vector of ratios of pairwise cutoff by SURE to the sum of all pairwise cutoffs by SURE
# max_tau, max_kappa, max_lambda  : upper bound for each tuning parameter 
# When D = 2, both "sure_pair" and "ratio_pair" are printed as vector of NAs

####################
#### for D = 2 #####
####################

set_grid_d2 <- function(Z_list, tol, nseq, qprob = 1.0){
  
  D <- length(Z_list) ; sure_indiv <- max_indiv <- ratio_indiv <- rep(NA, D)
  
  # calculate the cutoffs of singular value soft-thresholding chosen by SURE for each individual (sure_inidv)
  # calculate calculate the maximum singular values for each individual (max_indiv)
  
  for(d in 1:D){
    sure_indiv[d] <- find_SURE(Z_list[[d]], nlamb = 1000, tol)
    max_indiv[d] <- max(sgv_single(Z_list[[d]], tol)) 
  }
  
  # set up the grids for each tuning parameter 
  # lamb_seq : corresponding to the concatenated nuclear norm
  # tau_seq : common factor to the individual nuclear norms

  # calculate upper limit of lamb_seq
  max_sgv <- max(sgv_single(cbind(Z_list[[1]], Z_list[[2]]), tol))
  
  # calculate upper limit of tau_seq
  sum_indiv <- sum(sure_indiv) 
  max_tau <- max(max_indiv[1]*sum_indiv/sure_indiv[1], max_indiv[2]*sum_indiv/sure_indiv[2])
  
  # set up the grid 
  lamb_seq <- c(0, exp(seq(-5, log(qprob*max_sgv), len = nseq)))
  tau_seq <- c(0, exp(seq(-5, log(qprob*max_tau), len = nseq)))     
  
  # consider all the possible combinations of the above sequences
  df_tune <- tri_2D(tau_seq, lamb_seq, keep_close = F)
  colnames(df_tune) <- c("tau", "lambda")

  # ratio of the individual SURE for view 1 to the sum of the individual SUREs
  ratio_indiv[1] <- sure_indiv[1]/sum_indiv ; ratio_indiv[2] <- sure_indiv[2]/sum_indiv 
  
  return(list("df_tune" = df_tune, "sure_indiv" = sure_indiv, "ratio_indiv" = ratio_indiv, 
              "max_tau" = max_tau, "max_lambda" = max_sgv))
  
}

####################
#### for D = 3 #####
####################

set_grid_d3 <- function(Z_list, tol, nseq, qprob = 1.0){
  
  D <- length(Z_list) ; sure_indiv <- sure_pair <- max_indiv <- max_pair <- ratio_indiv <- ratio_pair <- rep(NA, D)
  
  # calculate the cutoffs of singular value soft-thresholding chosen by SURE for each individual (sure_inidv)
  # calculate calculate the maximum singular values for each individual (max_indiv)
  
  for(d in 1:D){
    sure_indiv[d] <- find_SURE(Z_list[[d]], nlamb = 1000, tol)
    max_indiv[d] <- max(sgv_single(Z_list[[d]], tol)) 
  }
  
  # calculate the cutoffs of singular value soft-thresholding chosen by SURE for each pairwise
  # 1st element : for view 1 & 2 // 2nd element : for view 1 & 3 // 3rd element : for view 2 & 3
  
  sure_pair[1] <- find_SURE(cbind(Z_list[[1]], Z_list[[2]]), nlamb = 1000, tol)
  sure_pair[2] <- find_SURE(cbind(Z_list[[1]], Z_list[[3]]), nlamb = 1000, tol)
  sure_pair[3] <- find_SURE(cbind(Z_list[[2]], Z_list[[3]]), nlamb = 1000, tol)
  
  # calculate calculate the maximum singular values for each pairwise
  # 1st element : for view 1 & 2 // 2nd element : for view 1 & 3 // 3rd element : for view 2 & 3
  max_pair[1] <- max(sgv_single(cbind(Z_list[[1]], Z_list[[2]]), tol)) 
  max_pair[2] <- max(sgv_single(cbind(Z_list[[1]], Z_list[[3]]), tol)) 
  max_pair[3] <- max(sgv_single(cbind(Z_list[[2]], Z_list[[3]]), tol)) 
  
  # set up the grids for each tuning parameter 
  # lamb_seq : corresponding to the concatenated nuclear norm
  # tau_seq : common factor to the individual nuclear norms
  # kappa_seq : common factor to the pairwise nuclear norms
  
  # calculate upper limit of lamb_seq
  max_sgv <- max(sgv_single(cbind(Z_list[[1]], Z_list[[2]], Z_list[[3]]), tol))
  
  # calculate upper limit of tau_seq
  sum_indiv <- sum(sure_indiv) 
  max_tau <- max(max_indiv[1]*sum_indiv/sure_indiv[1], max_indiv[2]*sum_indiv/sure_indiv[2], max_indiv[3]*sum_indiv/sure_indiv[3])
  
  # calculate upper limit of kappa_seq
  sum_pair <- sum(sure_pair)
  max_kappa <- max(max_pair[1]*sum_pair/sure_pair[1],max_pair[2]*sum_pair/sure_pair[2], max_pair[3]*sum_pair/sure_pair[3])

  # set up the grid 
  lamb_seq <- c(0, exp(seq(-5, log(qprob*max_sgv), len = nseq)))
  kappa_seq <- c(0, exp(seq(-5, log(qprob*max_kappa), len = nseq)))
  tau_seq <- c(0, exp(seq(-5, log(qprob*max_tau), len = nseq)))     
  
  # consider all the possible combinations of the above sequences
  df_tune <- tri_3D(tau_seq, kappa_seq, lamb_seq, keep_close = F)
  colnames(df_tune) <- c("tau", "kappa", "lambda")
  
  # ratio of the individual SURE for view 1 to the sum of the individual SUREs
  ratio_indiv[1] <- sure_indiv[1]/sum_indiv ; ratio_indiv[2] <- sure_indiv[2]/sum_indiv ; ratio_indiv[3] <- 1-ratio_indiv[1]-ratio_indiv[2]
  
  # ratio of the individual SURE for view 1 to the sum of the individual SUREs
  ratio_pair[1] <- sure_pair[1]/sum_pair ; ratio_pair[2] <- sure_pair[2]/sum_pair; ratio_pair[3] <- 1-ratio_pair[1]-ratio_pair[2]
  
  return(list("df_tune" = df_tune, "sure_indiv" = sure_indiv, "sure_pair" = sure_pair, 
              "ratio_indiv" = ratio_indiv, "ratio_pair" = ratio_pair,
              "max_tau" = max_tau, "max_kappa" = max_kappa, "max_lambda" = max_sgv))
  
}


##############################################
#### 5. function for individual refitting ####
##############################################

refit_indiv <- function(Z_list, Mtilde, pvec, tol){
  
  indexd <- proj_list <- list()
  n <- nrow(Mtilde) ; D <- length(pvec) ; cpvec <- c(0, cumsum(pvec))
  lobj <- tSVD_left(Mtilde, tol) ; rtilde <- lobj$r 
  
  for(d in 1:D){
    indexd[[d]] <- (cpvec[d]+1):(cpvec[d+1])
  }
  
  # conduct refitting and truncation as shown in the short note
  if(rtilde>=1){
    Uall <- lobj$u
    for(d in 1:D){
      viewd <- tSVD_left(Mtilde[, indexd[[d]]], tol) ; rd <- viewd$r ; Ud <- viewd$u
      if(rd>=1){
        # Get the new basis out
        Udnew <- tSVD_left(tcrossprod(Uall) %*% Ud, tol)$u
        # Refit individually 
        proj_list[[d]] <- tcrossprod(Udnew) %*% Z_list[[d]]
      } else{
        proj_list[[d]] <- matrix(0, nrow = n, ncol = length(indexd[[d]]))
      }
    }
  } else{
    for(d in 1:D){
      proj_list[[d]] <- matrix(0, nrow = n, ncol = length(indexd[[d]]))     
    }
  }
  
  return(proj_list)
  
}


##############################################################################################
##########  6. function to implement the DBFB algorithm and individual refitting #############
##############################################################################################

####################
#### for D = 2 #####
####################

algo1_d2 <- function(Z_list, tau, lambda, ratio_indiv, fb_vec, mean_mat, gam, max_iter, eps, tol){
  
  pvec <- c(ncol(Z_list[[1]]), ncol(Z_list[[2]]))
  Mhat_array <- array(NA, dim = c(nrow(Z_list[[1]]), sum(pvec), 2))
  
  # 1. Run the dual block forward-backward algorithm  
  fit <- DBFB_d2(cbind(Z_list[[1]], Z_list[[2]]), pvec, tau, lambda, ratio_indiv, gam, max_iter, eps) 
  # 2. Refitting  
  indiv_refit <- refit_indiv(Z_list, fit$Mtilde, pvec, tol) 
  # 3. Get back the final estimates after back-scaling and back-centering
  Mhat <- cbind(fb_vec[1]*indiv_refit[[1]], fb_vec[2]*indiv_refit[[2]]) + mean_mat  
  
  return(list(Mhat = Mhat, Mtilde = fit$Mtilde, num_iter = fit$num_iter))
}

####################
#### for D = 3 #####
####################

algo1_d3 <- function(Z_list, tau, kappa, lambda, ratio_pair, ratio_indiv, fb_vec, mean_mat, gam, max_iter, eps, tol){
  
  pvec <- c(ncol(Z_list[[1]]), ncol(Z_list[[2]]), ncol(Z_list[[3]]))
  Mhat_array <- array(NA, dim = c(nrow(Z_list[[1]]), sum(pvec), 3))
  # 1. Run the dual block forward-backward algorithm  
  fit <- DBFB_d3(cbind(Z_list[[1]], Z_list[[2]], Z_list[[3]]), pvec, tau, kappa, lambda, ratio_pair, ratio_indiv, gam, max_iter, eps) 
  # 2. Refitting  
  indiv_refit <- refit_indiv(Z_list, fit$Mtilde, pvec, tol) 
  # 3. Get back the final estimates after back-scaling and back-centering
  Mhat <- cbind(fb_vec[1]*indiv_refit[[1]], fb_vec[2]*indiv_refit[[2]], fb_vec[3]*indiv_refit[[3]]) + mean_mat     
  
  return(list(Mhat = Mhat, Mtilde = fit$Mtilde, num_iter = fit$num_iter))
}


##################################################
######### 7. Several functions for BCV  ##########
##################################################

########################################
### 7.1 functions to get BCV splits ####
########################################

## 1. Input
# Z_list : list where each element is standardized data for each view
# rfold : vector of row fold indices
# cfold_list : list having vector of column fold indices for each view

make_bcv_split <- function(X_list, rfold, cfold_list){
  
  ind_df <- expand.grid("j" = sort(unique(rfold)), "k" = sort(unique(cfold_list[[1]])))
  
  res_list <- list("X_jk" = list(), "X_njk" = list(), "X_jnk" = list(), "X_njnk" = list())
  
  for(b in 1:nrow(ind_df)){
    
    X_jk_list <- X_njk_list <- X_jnk_list <- X_njnk_list <- list() 
    
    j <- ind_df[b, 1] ; k <- ind_df[b, 2]
    
    for(d in 1:length(cfold_list)){
      cfold <- cfold_list[[d]] ; Xd <- X_list[[d]] 
      X_njnk_list[[d]] <- Xd[rfold!=j, cfold!=k] ; X_jk_list[[d]] <- Xd[rfold==j, cfold==k] 
      X_jnk_list[[d]] <- Xd[rfold==j, cfold!=k] ; X_njk_list[[d]] <- Xd[rfold!=j, cfold==k] 
    }
    
    res_list$X_jk[[b]] <- X_jk_list ; res_list$X_njk[[b]] <- X_njk_list 
    res_list$X_jnk[[b]] <- X_jnk_list ; res_list$X_njnk[[b]] <- X_njnk_list    
    
  }
  
  return(res_list) 
}


##################################################
#### 7.2 calculate BCV error for single-fold  ####
##################################################

####################
#### for D = 2 #####
####################

bcv_single_tune_d2 <- function(X_njnk_list, X_jnk_list, X_jk_list, X_njk_list, tau, lambda, ratio_indiv, gam, max_iter, eps, tol){
  
  train_pvec <- c()
  train_pvec[1] <- ncol(X_njnk_list[[1]]) ; train_pvec[2] <- ncol(X_njnk_list[[2]]) 
  cpvec <- c(0, cumsum(train_pvec))
  
  # standardize the training data X_njnk
  std_obj <- stdX(X_njnk_list, T) 
  mean_njnk_mat <- cbind(std_obj$mean_mat[[1]], std_obj$mean_mat[[2]]) 
  Z_njnk_list <- std_obj$Z_list ; fb_njnk_vec <- std_obj$fb_vec
  
  # get the estimator
  fit <- algo1_d2(Z_njnk_list, tau, lambda, ratio_indiv, fb_njnk_vec, mean_njnk_mat, gam, max_iter, eps, tol)  
  
  indiv_err <- c()
  
  for(d in 1:2){
    ginv_Mhat <- cal_ginv(fit$Mhat[,(cpvec[d]+1):(cpvec[d+1])], tol)
    denom <- sum((X_jk_list[[d]] - matrix(colMeans(X_jk_list[[d]]), nrow = nrow(X_jk_list[[d]]), ncol = ncol(X_jk_list[[d]]), byrow = T))**2)
    indiv_err[d] <- sum((X_jk_list[[d]]-X_jnk_list[[d]]%*%ginv_Mhat%*%X_njk_list[[d]])**2)/denom
  }
    
  return(indiv_err)
  
}



bcv_single_tune_d2_non <- function(X_njnk_list, X_jnk_list, X_jk_list, X_njk_list, tau, lambda, ratio_indiv, gam, max_iter, eps, tol){
  
  train_pvec <- c()
  train_pvec[1] <- ncol(X_njnk_list[[1]]) ; train_pvec[2] <- ncol(X_njnk_list[[2]]) 
  cpvec <- c(0, cumsum(train_pvec))
  
  # standardize the training data X_njnk
  std_obj <- stdX(X_njnk_list, T) 
  mean_njnk_mat <- cbind(std_obj$mean_mat[[1]], std_obj$mean_mat[[2]]) 
  Z_njnk_list <- std_obj$Z_list ; fb_njnk_vec <- std_obj$fb_vec
  
  # get the estimator
  fit <- algo1_d2(Z_njnk_list, tau, lambda, ratio_indiv, fb_njnk_vec, mean_njnk_mat, gam, max_iter, eps, tol)  
  
  indiv_err <- c()
  
  for(d in 1:2){
    ginv_Mhat <- cal_ginv(fit$Mhat[,(cpvec[d]+1):(cpvec[d+1])], tol)
    denom <- sum((X_jk_list[[d]])**2)
    indiv_err[d] <- sum((X_jk_list[[d]]-X_jnk_list[[d]]%*%ginv_Mhat%*%X_njk_list[[d]])**2)/denom
  }
  
  return(indiv_err)
  
}

####################
#### for D = 3 #####
####################


check_est_in_bcv_d3 <- function(X_njnk_list, M_njnk_list, tau, kappa, lambda, ratio_pair, ratio_indiv, gam, max_iter, eps, tol){
  
  train_pvec <- c()
  train_pvec[1] <- ncol(X_njnk_list[[1]]) ; train_pvec[2] <- ncol(X_njnk_list[[2]]) ; train_pvec[3] <- ncol(X_njnk_list[[3]]) 
  
  # standardize the training data X_njnk
  std_obj <- stdX(X_njnk_list, T) 
  mean_njnk_mat <- cbind(std_obj$mean_mat[[1]], std_obj$mean_mat[[2]], std_obj$mean_mat[[3]])
  Z_njnk_list <- std_obj$Z_list ; fb_njnk_vec <- std_obj$fb_vec
  
  # get the estimator
  fit <- algo1_d3(Z_njnk_list, tau, kappa, lambda, ratio_pair, ratio_indiv, fb_njnk_vec, mean_njnk_mat, gam, max_iter, eps, tol)

  return(cal_sfb(fit$Mhat, cbind(M_njnk_list[[1]],M_njnk_list[[2]],M_njnk_list[[3]]), train_pvec))
  
}


bcv_single_tune_d3 <- function(X_njnk_list, X_jnk_list, X_jk_list, X_njk_list, tau, kappa, lambda,  
                               ratio_pair, ratio_indiv, gam, max_iter, eps, tol){

  train_pvec <- c()
  train_pvec[1] <- ncol(X_njnk_list[[1]]) ; train_pvec[2] <- ncol(X_njnk_list[[2]]) ; train_pvec[3] <- ncol(X_njnk_list[[3]]) 
  cpvec <- c(0, cumsum(train_pvec))
  
  # standardize the training data X_njnk
  std_obj <- stdX(X_njnk_list, T) 
  mean_njnk_mat <- cbind(std_obj$mean_mat[[1]], std_obj$mean_mat[[2]], std_obj$mean_mat[[3]])
  Z_njnk_list <- std_obj$Z_list ; fb_njnk_vec <- std_obj$fb_vec
  
  # get the estimator
  fit <- algo1_d3(Z_njnk_list, tau, kappa, lambda, ratio_pair, ratio_indiv, fb_njnk_vec, mean_njnk_mat, gam, max_iter, eps, tol)
  
  # calculate BCV error
  indiv_err <- c()
  
  for(d in 1:3){
    ginv_Mhat <- cal_ginv(fit$Mhat[,(cpvec[d]+1):(cpvec[d+1])], tol)
    denom <- sum((X_jk_list[[d]] - matrix(colMeans(X_jk_list[[d]]), nrow = nrow(X_jk_list[[d]]), ncol = ncol(X_jk_list[[d]]), byrow = T))**2)
    indiv_err[d] <- sum((X_jk_list[[d]]-X_jnk_list[[d]]%*%ginv_Mhat%*%X_njk_list[[d]])**2)/denom
  }
  
  return(indiv_err)
  
}



##################################################
#### 7.3 choose the optimal tuning parameter  ####
##################################################

# bcv_err_array : nlamb by nrow_fold*ncol_fold by D array 

choose_opt <- function(bcv_err_array, rank_vec){
  
  nlamb <- length(rank_vec) ; avg_err <- se_vec <- c()
  
  for(l in 1:nlamb){
    fold_err <- apply(bcv_err_array[l,,], 1, mean, na.rm = T) 
    avg_err[l] <- mean(fold_err) ; se_vec[l] <- sqrt(var(fold_err)/length(fold_err))  
  }
  
  # we apply one-SE rule
  # if there are multiple tuning parameters with the same overall rank, choose the one with the smallest BCV error
  
  #set.seed(2022)
  #avg_err <- runif(25) ; rank_vec <- sample(c(1:10), 25, replace = T)
  #se_vec <- runif(25, min = 0, max = 0.2)
  
  opt_min <- order(avg_err)[1] # the one which minimizes the BCV error
  se_ind <- which(avg_err<=avg_err[opt_min] + se_vec[opt_min])
  small_rank <- se_ind[which(rank_vec[se_ind]==min(rank_vec[se_ind]))]
  
  if(sum(small_rank)>1){
    opt_se <- small_rank[which.min(avg_err[small_rank])]
  } else{
    opt_se <- small_rank
  }
  
  return(list("se" = opt_se, "min" = opt_min, "avg_err" = avg_err, "se_vec" = se_vec))
  
}

#######################################################################################################
######### function to provide scaled squared Frobenius norm, ranks and scaled chordal distances #######
#######################################################################################################

summary_fn <- function(est_mat, true_mat, pvec, tol){
  
  sfb <- sum(cal_sfb(est_mat, true_mat, pvec))
  sgvs <- sgv_all(est_mat, pvec, tol)
  ranks <- sapply(c(1:length(sgvs)), FUN = function(i){length(sgvs[[i]])})
  #dists <- chordal_dist(est_mat, true_mat, pvec, tol)  
  
  return(list("sfb" = sfb, "ranks" = ranks))
  
}
