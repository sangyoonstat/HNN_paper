
rm(list = ls())

library(r.jive)
load("./GTEx/GTEx_all.RData")
source("./functions/all_source_mac.R")

# standardize the simulated data
std <- stdX(X_list, center = F) ; Z_list <- std$Z_list ; tol <- 10^-3
fb_vec <- std$fb_vec # squared frobenius norms after column-centering
mean_mat <- cbind(std$mean_mat[[1]], std$mean_mat[[2]], std$mean_mat[[3]]) # matrix of column means
X_centered <- cbind(X_list[[1]], X_list[[2]], X_list[[3]]) - mean_mat
n <- nrow(X_list[[1]])
pvec <- c(ncol(X_list[[1]]), ncol(X_list[[2]]), ncol(X_list[[3]]))


####################
###### JIVE ########
####################

Z_list2 <- list() 

for(d in 1:length(pvec)){
  Z_list2[[d]] <- t(Z_list[[d]])
}

jive_out <- jive(Z_list2, method = "perm", scale = F, est = T, showProgress = F)

jive_M1hat <- t(jive_out$joint[[1]] + jive_out$individual[[1]])
jive_M2hat <- t(jive_out$joint[[2]] + jive_out$individual[[2]])
jive_M3hat <- t(jive_out$joint[[3]] + jive_out$individual[[3]])

jive_Mhat <- cbind(fb_vec[1]*jive_M1hat, fb_vec[2]*jive_M2hat, fb_vec[3]*jive_M3hat)


save(jive_Mhat, file = "./GTEx/fit_JIVE.RData")


#################################################
###### 1. Apply JIVE to 1st and 2nd views #######
#################################################


Z_list2 <- list() 
Z_list2[[1]] <- t(Z_list[[1]]) ; Z_list2[[2]] <- t(Z_list[[2]])
jive_out <- jive(Z_list2, method = "perm", scale = F, est = T, showProgress = F)

length(sgv_single(jive_out$joint[[1]], tol)) # joint rank between view 1 and view 2
length(sgv_single(jive_out$individual[[1]], tol)) # individual rank of view 1
length(sgv_single(jive_out$individual[[2]], tol)) # individual rank of view 2

#################################################
###### 1. Apply JIVE to 1st and 3rd views #######
#################################################

Z_list2 <- list() 
Z_list2[[1]] <- t(Z_list[[1]]) ; Z_list2[[2]] <- t(Z_list[[3]])
jive_out <- jive(Z_list2, method = "perm", scale = F, est = T, showProgress = F)
length(sgv_single(jive_out$joint[[1]], tol)) # joint rank between view 1 and view 3
length(sgv_single(jive_out$individual[[1]], tol)) # individual rank of view 1
length(sgv_single(jive_out$individual[[2]], tol)) # individual rank of view 3


#################################################
###### 3. Apply JIVE to 2nd and 3rd views #######
#################################################

Z_list2 <- list() 
Z_list2[[1]] <- t(Z_list[[2]]) ; Z_list2[[2]] <- t(Z_list[[3]])
jive_out <- jive(Z_list2, method = "perm", scale = F, est = T, showProgress = F)
length(sgv_single(jive_out$joint[[1]], tol)) # joint rank between view 2 and view 3
length(sgv_single(jive_out$individual[[1]], tol)) # individual rank of view 2
length(sgv_single(jive_out$individual[[2]], tol)) # individual rank of view 3



