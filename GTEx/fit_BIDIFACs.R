

rm(list = ls())

load("./GTEx/GTEx_all.RData")
source("./functions/my_unifacs.R")
source("./functions/all_source_mac.R")

tol <- 10^-3

center = F ; eps = 10^-5 ; max_iter = 1000 ; seednum = NULL

#############################
######### 1.UNIFAC ##########
#############################

#pdat <- prep_dat(X_list, center)
#sigma_vec <- pdat$sigma_vec ; mean_mat <- pdat$mean_mat

sigma_vec <- rep(1, 3)
X00 <- cbind(X_list[[1]], X_list[[2]], X_list[[3]]) ; n <- pdat$n ; pvec <- pdat$pvec ; D <- length(pvec)

# Create index for each view
start.ind.p <- c(1, cumsum(pvec)[1:(D-1)] + 1) ; end.ind.p <- cumsum(pvec)

# Set the values of the lambdas
lambda.R <- sqrt(n) + sqrt(sum(pvec)) ; lambda.I <- sqrt(n) + sqrt(pvec) 

# Initialization 

R00 <- replicate(sum(pvec), rnorm(n)) ; I00 <- replicate(sum(pvec), rnorm(n))

# Run the ISSVT algorithm
R00.nuc <- NA ; I00.nuc <-c() ; bool <- TRUE ; count <- 1 ; crit0 <- 0

while(bool){
  
  crit0.old <- crit0
  
  # Update R
  sft1 <- softSVD(X00-I00, lambda.R) ; R00 <- sft1$out; R00.nuc <- sft1$nuc
  
  # Update I
  for (j in 1:D){
    ind2 <- start.ind.p[j]:end.ind.p[j]
    sft2 <- softSVD(X00[,ind2]-R00[,ind2], lambda.I[j])
    I00[,ind2] <- sft2$out ; I00.nuc[j] <- sft2$nuc
  }
  
  crit0 = frob_sq(X00-R00-I00) + 2*sum(lambda.R*R00.nuc) + 2*sum(lambda.I*I00.nuc)
  
  if (abs(crit0.old-crit0)<eps){ bool=FALSE }
  else if (count==max_iter){ bool=FALSE}
  else{ count = count+1 }
  
}


S00_list <- R00_list <- I00_list <- list()

for (j in 1:D){
  ind2 <- start.ind.p[j]:end.ind.p[j]
  R00_list[[j]] <- R00[,ind2]*sigma_vec[j] 
  I00_list[[j]] <- I00[,ind2]*sigma_vec[j] 
  S00_list[[j]] <- R00_list[[j]] + I00_list[[j]] 
  #+ mean_mat[,ind2]
}

unifac_Mhat <- cbind(S00_list[[1]], S00_list[[2]], S00_list[[3]])

###################################
######### 2.UNIFAC_refit ##########
###################################

unifac_refit_Mhat <- refit_indiv(X_list, unifac_Mhat, pvec, tol)


###############################
######### 3.UNIFAC+  ##########
###############################


p.ind <- list()
p.ind[[1]] <- c(1:pvec[1])
p.ind[[2]] <- c((pvec[1]+1):(pvec[1]+pvec[2]))
p.ind[[3]] <- c((pvec[1]+pvec[2]+1):(pvec[1]+pvec[2]+pvec[3]))

p.ind.list <- list()
p.ind.list[[1]] <- unlist(p.ind[1:3])
p.ind.list[[2]] <- unlist(p.ind[1:2])
p.ind.list[[3]] <- unlist(p.ind[c(1,3)])
p.ind.list[[4]] <- unlist(p.ind[c(2,3)])
p.ind.list[[5]] <- unlist(p.ind[c(1)])
p.ind.list[[6]] <- unlist(p.ind[c(2)])
p.ind.list[[7]] <- unlist(p.ind[c(3)])


res.g <- unifac.plus.given(X0 = t(X00), p.ind = p.ind, p.ind.list = p.ind.list)

for(k in 1:length(res.g$S)){
  if(k==1){
    unifac_plus_mat <- res.g$S[[1]]
  } else{
    unifac_plus_mat <- unifac_plus_mat + res.g$S[[k]]
  }
}

unifac_plus_Mhat <- cbind(sigma_vec[1]*t(unifac_plus_mat)[,p.ind[[1]]],
                          sigma_vec[2]*t(unifac_plus_mat)[,p.ind[[2]]],
                          sigma_vec[3]*t(unifac_plus_mat)[,p.ind[[3]]])


####################################
######### 4.UNIFAC+_refit ##########
####################################

unifac_plus_Mtilde <- cbind(t(unifac_plus_mat)[,p.ind[[1]]], t(unifac_plus_mat)[,p.ind[[2]]], t(unifac_plus_mat)[,p.ind[[3]]])

unifac_plus_Mhat <- refit_indiv(X_list, unifac_plus_Mtilde, pvec, tol)

save(unifac_Mhat, unifac_refit_Mhat, unifac_plus_Mhat, unifac_plus_Mhat, file = "~/Desktop/GTEx/fit_BIDIFACs.RData")
