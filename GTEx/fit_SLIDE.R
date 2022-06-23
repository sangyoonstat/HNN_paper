
rm(list = ls())

seednum <- 2021 ; library(SLIDE) 

load("./GTEx/GTEx_all.RData")
source("./functions/all_source_mac.R")

# standardize the simulated data
std <- stdX(X_list, center = F) ; Z_list <- std$Z_list ; tol <- 10^-3
fb_vec <- std$fb_vec # squared frobenius norms after column-centering
#mean_mat <- cbind(std$mean_mat[[1]], std$mean_mat[[2]], std$mean_mat[[3]]) # matrix of column means
#X_centered <- cbind(X_list[[1]], X_list[[2]], X_list[[3]]) - mean_mat
n <- nrow(X_list[[1]])
pvec <- c(ncol(X_list[[1]]), ncol(X_list[[2]]), ncol(X_list[[3]]))


tol <- 10^-3 ; nrow_fold <- ncol_fold <- 2 


out_slide = tryCatch(slide(cbind(Z_list[[1]], Z_list[[2]], Z_list[[3]]), pvec = pvec,
                           n_fold = 2, p_fold = 2, center = F), error = function(e) NA)


Mtilde <- out_slide$model$U%*%t(out_slide$model$V) 

slide_Mhat <- cbind(fb_vec[1]*Mtilde[,c(1:pvec[1])], 
                    fb_vec[2]*Mtilde[,c((pvec[1]+1):(pvec[1]+pvec[2]))],
                    fb_vec[3]*Mtilde[,-c(1:(pvec[1]+pvec[2]))])

save(slide_Mhat, file = "./GTEx/fit_SLIDE.RData")
