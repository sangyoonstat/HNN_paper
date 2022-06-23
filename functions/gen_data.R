


gen_vec <- function(v, theta){
  m <- length(v) ; w <- runif(m) ; w <- w-mean(w)
  w <- ((diag(1,m) - tcrossprod(v))%*%w)/sqrt(sum(((diag(1,n) - tcrossprod(v))%*%w)**2))
  return(cos(theta)*v + sqrt(1-cos(theta)^2)*w)
}


gen_two_scores <- function(n, theta){
  
  U_list <- list()
  
  r <- length(theta)
  ## 1. Fix the score corresponding to the first view
  U1 <- scale(matrix(runif(n), nrow = n, ncol = r), scale = F) ; U_list[[1]] <- qr.Q(qr(U1))  
  
  ## 2. Generate the score for the second view
  for(k in 1:r){
    w <- gen_vec(U_list[[1]][,k], (pi/180)*theta[k])
    if(k==1){
      U2 <- w  
    } else{
      U2 <- cbind(U2, w)
    }
  }
  
  U_list[[2]] <- qr.Q(qr(U2))
  
  return(U_list)
  
}


gen_others_d2 <- function(r0, U_list){
  n <- nrow(U_list[[1]]) ; U1 <- U_list[[1]] ; U2 <- U_list[[2]] ; Umat <- cbind(U1, U2)
  # generate the joint scores
  U0 <- scale(matrix(runif(n*r0), nrow = n, ncol = r0), scale = F)
  U0 <- qr.Q(qr((diag(1, n) - Umat%*%solve(crossprod(Umat))%*%t(Umat))%*%U0))
  return(list("U0" = U0, "U1" = U1, "U2" = U2))
}



gen_pair_scores <- function(n, theta1, theta2){
  
  U_list <- list()
  
  ## 1. Fix the score corresponding to the partially-shared scores of the 1st and 2nd views
  r <- length(theta1)
  U12 <- scale(matrix(runif(n), nrow = n, ncol = length(theta1)), scale = F) 
  U_list[[1]] <- qr.Q(qr(U12)) 
  
  ## 2. Generate the partially-shared scores of the 1st and 3rd views
  for(k in 1:r){
    w <- gen_vec(U_list[[1]][,k], (pi/180)*theta1[k])
    if(k==1){
      U13 <- w  
    } else{
      U13 <- cbind(U13, w)
    }
  }
  
  U_list[[2]] <- qr.Q(qr(U13))
  
  ## 3. Generate the partially-shared scores of the 2nd and 3rd views
  for(k in 1:r){
    w <- gen_vec(U_list[[2]][,k], (pi/180)*theta2[k])
    if(k==1){
      U23 <- w  
    } else{
      U23 <- cbind(U23, w)
    }
  }
  
  U_list[[3]] <- qr.Q(qr(U23))
  
  return(U_list)
  
}


gen_others_d3 <- function(r0, r1, r2, r3, U_list){
  
  n <- nrow(U_list[[1]]) ; U12 <- U_list[[1]] ; U13 <- U_list[[2]] ; U23 <- U_list[[3]]
  
  # generate the joint scores
  Umat <- cbind(U12, U13, U23)
  U0 <- scale(matrix(runif(n*r0), nrow = n, ncol = r0), scale = F)
  U0 <- qr.Q(qr((diag(1, n) - Umat%*%solve(crossprod(Umat))%*%t(Umat))%*%U0))
  
  # generate the individual scores for the 1st view
  Umat <- cbind(U0, U12, U13)
  U1 <- scale(matrix(runif(n*r1), nrow = n, ncol = r1), scale = F)
  U1 <- qr.Q(qr((diag(1, n) - Umat%*%solve(crossprod(Umat))%*%t(Umat))%*%U1))
  
  # generate the individual scores for the 2nd view
  Umat <- cbind(U0, U12, U23)
  U2 <- scale(matrix(runif(n*r2), nrow = n, ncol = r2), scale = F)
  U2 <- qr.Q(qr((diag(1, n) - Umat%*%solve(crossprod(Umat))%*%t(Umat))%*%U2))
  
  # generate the individual scores for the 3rd view
  Umat <- cbind(U0, U13, U23)
  U3 <- scale(matrix(runif(n*r3), nrow = n, ncol = r3), scale = F)
  U3 <- qr.Q(qr((diag(1, n) - Umat%*%solve(crossprod(Umat))%*%t(Umat))%*%U3))  
  
  return(list("U0" = U0, "U12" = U12, "U13" = U13, "U23" = U23, "U1" = U1, "U2" = U2, "U3" = U3))
}


################################################################################################################################################################


dgen_d2_givenU <- function(U_list, pvec, sjoint, s1ind, s2ind, snr_vec, orthogonalV = T){
  
  for(i in 1:length(U_list)){
    if(i==1){
      U <- U_list[[i]]
    } else{
      U <- cbind(U, U_list[[i]])
    }
  }
  
  n <- nrow(U) ; p_total <- sum(pvec) ; D <- length(pvec) 
  r0 <- length(sjoint) ; rvec <- c(length(s1ind), length(s2ind))
  r_total <- r0 + sum(rvec)
  
  # generate noise part of X
  X <- matrix(0, n, p_total)
  
  # generate dataset-specific loadings
  Vtmp <- list()
  for (d in 1:D){
    Vtmp[[d]] <- matrix(runif(pvec[d]*(r0+rvec[d])), pvec[d], (r0+rvec[d]))
    # orthogonalize the loadings
    Vtmp[[d]] <- qr.Q(qr(Vtmp[[d]]))
  }
  
  # Put everything together into one V
  V <- matrix(0, nrow = p_total, ncol = r_total)
  
  for (d in 1:D){
    
    # Row index
    if (d == 1){
      index <- c(1:pvec[1])
    } else{
      index <- c((pvec[1]+1):p_total)
    }
    
    # Column index
    if (d == 1){
      col_indexd = c(1:(r0+rvec[1]))
    } else{
      col_indexd = c(1:r0,(r0+rvec[1]+1):(r0+sum(rvec)))
    }
    
    V[index,col_indexd] = Vtmp[[d]]
    
  }
  
  # Orthogonallize V
  if (orthogonalV == T){
    V <- qr.Q(qr(V))
  }
  
  # Form true M
  V <- V%*%diag(c(sjoint, s1ind, s2ind))
  X <- tcrossprod(U, V) ; true_M <- X
  
  # Datasets have different scales or different dimensions
  sigmavec = rep(0, D) ; err_list <- list()
  
  for (d in 1:D){
    # Row index for dth dataset
    if (d == 1){
      index <- c(1:pvec[1])
    } else{
      index <- c((pvec[1]+1):p_total)
    }
    
    sigmavec[d] <- sqrt(sum(X[,index]^2) / (n * pvec[d] * snr_vec[d]))
    err_list[[d]] <- matrix(rnorm(n * pvec[d], sd = sigmavec[d]), n, pvec[d])
    X[,index] <- X[,index] + err_list[[d]]
  }
  
  # Combine joint loadings and individual loadings into one matrix V
  return(list(X = X, pvec = pvec, sigmavec = sigmavec, true_M = true_M, err_list = err_list)) 
}


dgen_d3_givenU <- function(U_list, pvec, sjoint, s12, s13, s23, s1ind, s2ind, s3ind, snr_vec, orthogonalV = T){
  
  for(i in 1:length(U_list)){
    if(i==1){
      U <- U_list[[i]]
    } else{
      U <- cbind(U, U_list[[i]])
    }
  }
  
  n <- nrow(U) ; p_total <- sum(pvec) ; pcum <- cumsum(pvec) ; D <- length(pvec) 
  r0 <- length(sjoint) ; rmat <- matrix(NA, nrow = D, ncol = D)
  
  diag(rmat) <- c(length(s1ind), length(s2ind), length(s3ind))
  rmat[1,2] <- rmat[2,1] <- length(s12)
  rmat[1,3] <- rmat[3,1] <- length(s13)
  rmat[2,3] <- rmat[3,2] <- length(s23)
  
  r_total <- r0 + sum(rmat[upper.tri(rmat, diag=T)])
  
  # generate noise part of X
  X <- matrix(0, n, p_total)
  
  # generate dataset-specific loadings
  Vtmp <- list()
  for (d in 1:D){
    Vtmp[[d]] <- matrix(runif(pvec[d]*(r0+sum(rmat[d,]))), pvec[d],(r0+sum(rmat[d,])))
    # orthogonalize the loadings
    Vtmp[[d]] <- qr.Q(qr(Vtmp[[d]]))
  }
  
  # Put everything together into one V
  V <- matrix(0, nrow = p_total, ncol = r_total)
  
  for (d in 1:D){
    
    # Row index
    if (d == 1){
      index <- c(1:pvec[1])
    } else{
      index <- c((pcum[d-1]+1):pcum[d])
    }
    
    # Column index, order is r0, r12, r13, r23, r1, r2, r3
    if (d==1){
      col_indexd = c(1:(r0+rmat[1,2]+rmat[1,3]),(r0+rmat[1,2]+rmat[1,3]+rmat[2,3]+1):(r0+rmat[1,2]+rmat[1,3]+rmat[2,3]+rmat[1,1]))
    } else if (d==2){
      col_indexd = c(1:(r0+rmat[1,2]),(r0+rmat[1,2]+rmat[1,3]+1):(r0+rmat[1,2]+rmat[1,3]+rmat[2,3]), (r0+rmat[1,2]+rmat[1,3]+rmat[2,3]+rmat[1,1]+1):(r0+rmat[1,2]+rmat[1,3]+rmat[2,3]+rmat[1,1]+rmat[2,2]))
    } else{
      col_indexd = c(1:r0,(r0+rmat[1,2]+1):(r0+rmat[1,2]+rmat[1,3]+rmat[2,3]),(r0+rmat[1,2]+rmat[1,3]+rmat[2,3]+rmat[1,1]+rmat[2,2]+1):(r0+rmat[1,2]+rmat[1,3]+rmat[2,3]+rmat[1,1]+rmat[2,2]+rmat[3,3]))
    }
    V[index,col_indexd]=Vtmp[[d]]
  }
  
  # Orthogonallize V
  if (orthogonalV == T){
    V <- qr.Q(qr(V))
  }
  # Add the scale between joint, semi-joint and individual by fixing singular values
  # Column index, order is r0, r12, r13, r23, r1, r2, r3
  V <- V%*%diag(c(sjoint,s12, s13, s23, s1ind, s2ind, s3ind))
  
  # Form X
  X <- tcrossprod(U, V) ; true_M <- X
  
  # Datasets have different scales or different dimensions
  sigmavec = rep(0, D) ; err_list <- list()
  
  for (d in 1:D){
    # Row index for dth dataset
    if (d == 1){
      index <- c(1:pvec[1])
    } else{
      index <- c((pcum[d-1]+1):pcum[d])
    }
    sigmavec[d] <- sqrt(sum(X[,index]^2) / (n * pvec[d] * snr_vec[d]))
    err_list[[d]] <- matrix(rnorm(n * pvec[d], sd = sigmavec[d]), n, pvec[d])
    X[,index] <- X[,index] + err_list[[d]]
  }
  
  # Combine joint loadings and individual loadings into one matrix V
  return(list(X = X, pvec = pvec, sigmavec = sigmavec, true_M = true_M, err_list = err_list)) 
}



################################################################################################################################################################



dgen_orth_d2 <- function(n, pvec, sjoint, s1ind, s2ind, snr_vec, orthogonalV = T){
  
  p_total <- sum(pvec) ; D <- length(pvec) 
  r0 <- length(sjoint) ; rvec <- c(length(s1ind), length(s2ind))
  r_total <- r0 + sum(rvec)
  
  # generate noise part of X
  X <- matrix(0, n, p_total)
  
  # Column-center U so that the overall mean is zero
  U = matrix(runif(n*r_total), n, r_total)
  U = scale(U, scale = F)
  
  # Orthogonalize the scores
  U <- qr.Q(qr(U))
  
  # generate dataset-specific loadings
  Vtmp <- list()
  for (d in 1:D){
    Vtmp[[d]] <- matrix(runif(pvec[d]*(r0+rvec[d])), pvec[d], (r0+rvec[d]))
    # orthogonalize the loadings
    Vtmp[[d]] <- qr.Q(qr(Vtmp[[d]]))
  }
  
  # Put everything together into one V
  V <- matrix(0, nrow = p_total, ncol = r_total)
  
  for (d in 1:D){
    
    # Row index
    if (d == 1){
      index <- c(1:pvec[1])
    } else{
      index <- c((pvec[1]+1):p_total)
    }
    
    # Column index
    if (d == 1){
      col_indexd = c(1:(r0+rvec[1]))
    } else{
      col_indexd = c(1:r0,(r0+rvec[1]+1):(r0+sum(rvec)))
    }
    
    V[index,col_indexd] = Vtmp[[d]]
    
  }
  
  # Orthogonallize V
  if (orthogonalV == T){
    V <- qr.Q(qr(V))
  }
  
  # Form true M
  V <- V%*%diag(c(sjoint, s1ind, s2ind))
  X <- tcrossprod(U, V) ; true_M <- X
  
  # Datasets have different scales or different dimensions
  sigmavec = rep(0, D) ; err_list <- list()
  
  for (d in 1:D){
    # Row index for dth dataset
    if (d == 1){
      index <- c(1:pvec[1])
    } else{
      index <- c((pvec[1]+1):p_total)
    }
    
    sigmavec[d] <- sqrt(sum(X[,index]^2) / (n * pvec[d] * snr_vec[d]))
    err_list[[d]] <- matrix(rnorm(n * pvec[d], sd = sigmavec[d]), n, pvec[d])
    X[,index] <- X[,index] + err_list[[d]]
  }
  
  # Combine joint loadings and individual loadings into one matrix V
  return(list(X = X, pvec = pvec, sigmavec = sigmavec, true_M = true_M, U = U, err_list = err_list)) 
  
}



dgen_d2 <- function(n, pvec, sjoint, s1ind, s2ind, snr_vec, orthogonalV = T){
  
  p_total <- sum(pvec) ; D <- length(pvec) 
  r0 <- length(sjoint) ; r1 <- length(s1ind) ; r2 <- length(s2ind) ; rvec <- c(r1, r2)
  r_total <- r0 + sum(rvec)
  
  # generate noise part of X
  X <- matrix(0, n, p_total)
  
  # Column-center U so that the overall mean is zero
  U0 <- scale(matrix(runif(n*r0, 0, 1), n, r0), scale = F)
  U1 <- scale(matrix(runif(n*r1, 0, 1), n, r1), scale = F)
  U2 <- scale(matrix(runif(n*r2, 0, 1), n, r2), scale = F)
  
  conc_U0_U1 <- qr.Q(qr(cbind(U0, U1))) ; conc_U0_U2 <- qr.Q(qr(cbind(U0, U2))) 
  
  U0 <- as.matrix(conc_U0_U1[,c(1:r0)])
  U1 <- as.matrix(conc_U0_U1[,-c(1:r0)])
  U2 <- as.matrix(conc_U0_U2[,-c(1:r0)])
  
  U <- cbind(U0, U1, U2)
  
  # generate dataset-specific loadings
  Vtmp <- list()
  for (d in 1:D){
    Vtmp[[d]] <- matrix(runif(pvec[d]*(r0+rvec[d])), pvec[d], (r0+rvec[d]))
    # orthogonalize the loadings
    Vtmp[[d]] <- qr.Q(qr(Vtmp[[d]]))
  }
  
  # Put everything together into one V
  V <- matrix(0, nrow = p_total, ncol = r_total)
  
  for (d in 1:D){
    
    # Row index
    if (d == 1){
      index <- c(1:pvec[1])
    } else{
      index <- c((pvec[1]+1):p_total)
    }
    
    # Column index
    if (d == 1){
      col_indexd = c(1:(r0+rvec[1]))
    } else{
      col_indexd = c(1:r0,(r0+rvec[1]+1):(r0+sum(rvec)))
    }
    
    V[index,col_indexd] = Vtmp[[d]]
    
  }
  
  # Orthogonallize V
  if (orthogonalV == T){
    V <- qr.Q(qr(V))
  }
  
  # Form true M
  V <- V%*%diag(c(sjoint, s1ind, s2ind))
  X <- tcrossprod(U, V) ; true_M <- X
  
  # Datasets have different scales or different dimensions
  sigmavec = rep(0, D) ; err_list <- list()
  
  for (d in 1:D){
    # Row index for dth dataset
    if (d == 1){
      index <- c(1:pvec[1])
    } else{
      index <- c((pvec[1]+1):p_total)
    }
    
    sigmavec[d] <- sqrt(sum(X[,index]^2) / (n * pvec[d] * snr_vec[d]))
    err_list[[d]] <- matrix(rnorm(n * pvec[d], sd = sigmavec[d]), n, pvec[d])
    X[,index] <- X[,index] + err_list[[d]]
  }
  
  # Combine joint loadings and individual loadings into one matrix V
  return(list(X = X, pvec = pvec, sigmavec = sigmavec, true_M = true_M, U = U, err_list = err_list)) 
  
}



################################################################################################################################################################



dgen_orth_d3 <- function(n, pvec, sjoint, s12, s13, s23, s1ind, s2ind, s3ind, snr_vec, scaleU = T, orthogonalV = T){
  
  p_total <- sum(pvec) ; pcum <- cumsum(pvec) ; cvec <- c(1,1,1) ; D <- 3 
  r0 <- length(sjoint) ; rmat <- matrix(NA, nrow = 3, ncol = 3)
  
  diag(rmat) <- c(length(s1ind), length(s2ind), length(s3ind))
  rmat[1,2] <- rmat[2,1] <- length(s12)
  rmat[1,3] <- rmat[3,1] <- length(s13)
  rmat[2,3] <- rmat[3,2] <- length(s23)
  
  r_total <- r0 + sum(rmat[upper.tri(rmat, diag=T)])
  
  # generate noise part of X
  X <- matrix(0, n, p_total)
  
  # generate all the scores
  U0 <- matrix(runif(n*r0, 0, 1), n, r0) 
  U12 <- matrix(runif(n*rmat[1,2], 0, 1), n, rmat[1,2]) 
  U13 <- matrix(runif(n*rmat[1,3], 0, 1), n, rmat[1,3]) 
  U23 <- matrix(runif(n*rmat[2,3], 0, 1), n, rmat[2,3]) 
  U1 <- matrix(runif(n*rmat[1,1], 0, 1), n, rmat[1,1]) 
  U2 <- matrix(runif(n*rmat[2,2], 0, 1), n, rmat[2,2]) 
  U3 <- matrix(runif(n*rmat[3,3], 0, 1), n, rmat[3,3]) 
  
  # Column-center U so that the overall mean is zero
  if(scaleU){
    U0 <- scale(U0, scale = F) ; U12 <- scale(U12, scale = F) ; U13 <- scale(U13, scale = F) ; U23 <- scale(U23, scale = F) 
    U1 <- scale(U1, scale = F) ; U2 <- scale(U2, scale = F) ; U3 <- scale(U3, scale = F)
  }  
  
  U <- qr.Q(qr(cbind(U0, U12, U13, U23, U1, U2, U3)))
  
  # generate dataset-specific loadings
  Vtmp <- list()
  for (d in 1:D){
    Vtmp[[d]] <- matrix(runif(pvec[d]*(r0+sum(rmat[d,]))), pvec[d],(r0+sum(rmat[d,])))
    # orthogonalize the loadings
    Vtmp[[d]] <- qr.Q(qr(Vtmp[[d]]))
  }
  
  # Put everything together into one V
  V <- matrix(0, nrow = p_total, ncol = r_total)
  
  for (d in 1:D){
    
    # Row index
    if (d == 1){
      index <- c(1:pvec[1])
    } else{
      index <- c((pcum[d-1]+1):pcum[d])
    }
    
    # Column index, order is r0, r12, r13, r23, r1, r2, r3
    if (d==1){
      col_indexd = c(1:(r0+rmat[1,2]+rmat[1,3]),(r0+rmat[1,2]+rmat[1,3]+rmat[2,3]+1):(r0+rmat[1,2]+rmat[1,3]+rmat[2,3]+rmat[1,1]))
    } else if (d==2){
      col_indexd = c(1:(r0+rmat[1,2]),(r0+rmat[1,2]+rmat[1,3]+1):(r0+rmat[1,2]+rmat[1,3]+rmat[2,3]), (r0+rmat[1,2]+rmat[1,3]+rmat[2,3]+rmat[1,1]+1):(r0+rmat[1,2]+rmat[1,3]+rmat[2,3]+rmat[1,1]+rmat[2,2]))
    } else{
      col_indexd = c(1:r0,(r0+rmat[1,2]+1):(r0+rmat[1,2]+rmat[1,3]+rmat[2,3]),(r0+rmat[1,2]+rmat[1,3]+rmat[2,3]+rmat[1,1]+rmat[2,2]+1):(r0+rmat[1,2]+rmat[1,3]+rmat[2,3]+rmat[1,1]+rmat[2,2]+rmat[3,3]))
    }
    V[index,col_indexd]=Vtmp[[d]]
  }
  
  # Orthogonallize V
  if (orthogonalV == T){
    V <- qr.Q(qr(V))
  }
  # Add the scale between joint, semi-joint and individual by fixing singular values
  # Column index, order is r0, r12, r13, r23, r1, r2, r3
  V <- V%*%diag(c(sjoint,s12, s13, s23, s1ind, s2ind, s3ind))
  
  # Form X
  X <- tcrossprod(U, V) ; true_M <- X
  
  # Datasets have different scales or different dimensions
  sigmavec = rep(0, D) ; err_list <- list()
  
  for (d in 1:D){
    # Row index for dth dataset
    if (d == 1){
      index <- c(1:pvec[1])
    } else{
      index <- c((pcum[d-1]+1):pcum[d])
    }
    sigmavec[d] <- sqrt(sum(X[,index]^2) / (n * pvec[d] * snr_vec[d]))
    err_list[[d]] <- matrix(rnorm(n * pvec[d], sd = sigmavec[d]), n, pvec[d])
    X[,index] <- X[,index] + err_list[[d]]
  }
  
  # Combine joint loadings and individual loadings into one matrix V
  return(list(X = X, pvec = pvec, sigmavec = sigmavec, U = U, V = V, true_M = true_M, 
              err_list = err_list)) 
}




dgen_d3 <- function(n, pvec, sjoint, s12, s13, s23, s1ind, s2ind, s3ind, snr_vec, orth_indiv_shared = T, scaleU = T, orthogonalV = T){
  
  p_total <- sum(pvec) ; pcum <- cumsum(pvec) ; cvec <- c(1,1,1) ; D <- 3 
  r0 <- length(sjoint) ; rmat <- matrix(NA, nrow = 3, ncol = 3)
  
  diag(rmat) <- c(length(s1ind), length(s2ind), length(s3ind))
  rmat[1,2] <- rmat[2,1] <- length(s12)
  rmat[1,3] <- rmat[3,1] <- length(s13)
  rmat[2,3] <- rmat[3,2] <- length(s23)
  
  r_total <- r0 + sum(rmat[upper.tri(rmat, diag=T)])
  
  # generate noise part of X
  X <- matrix(0, n, p_total)
  
  # generate all the scores
  U0 <- matrix(runif(n*r0, 0, 1), n, r0) 
  U12 <- matrix(runif(n*rmat[1,2], 0, 1), n, rmat[1,2]) 
  U13 <- matrix(runif(n*rmat[1,3], 0, 1), n, rmat[1,3]) 
  U23 <- matrix(runif(n*rmat[2,3], 0, 1), n, rmat[2,3]) 
  U1 <- matrix(runif(n*rmat[1,1], 0, 1), n, rmat[1,1]) 
  U2 <- matrix(runif(n*rmat[2,2], 0, 1), n, rmat[2,2]) 
  U3 <- matrix(runif(n*rmat[3,3], 0, 1), n, rmat[3,3]) 
  
  # Column-center U so that the overall mean is zero
  if(scaleU){
    U0 <- scale(U0, scale = F) ; U12 <- scale(U12, scale = F) ; U13 <- scale(U13, scale = F) ; U23 <- scale(U23, scale = F) 
    U1 <- scale(U1, scale = F) ; U2 <- scale(U2, scale = F) ; U3 <- scale(U3, scale = F)
  }  
  
  # orthogonality between 
  # 1. globally shared scores and individual scores
  # 2. globally shared scores and pairwisely shared scores
  conc_U0_U1 <- qr.Q(qr(cbind(U0, U1))) ; U1 <- as.matrix(conc_U0_U1[,-c(1:r0)])
  conc_U0_U2 <- qr.Q(qr(cbind(U0, U2))) ; U2 <- as.matrix(conc_U0_U2[,-c(1:r0)])
  conc_U0_U3 <- qr.Q(qr(cbind(U0, U3))) ; U3 <- as.matrix(conc_U0_U3[,-c(1:r0)])
  U0 <- as.matrix(conc_U0_U1[,c(1:r0)])
  
  if(orth_indiv_shared){
    conc_U0_U12 <- qr.Q(qr(cbind(U0, U1, U2, U12))) ; U12 <- as.matrix(conc_U0_U12[,-c(1:(r0+rmat[1,1]+rmat[2,2]))])
    conc_U0_U13 <- qr.Q(qr(cbind(U0, U1, U3, U13))) ; U13 <- as.matrix(conc_U0_U13[,-c(1:(r0+rmat[1,1]+rmat[3,3]))])
    conc_U0_U23 <- qr.Q(qr(cbind(U0, U2, U3, U23))) ; U23 <- as.matrix(conc_U0_U23[,-c(1:(r0+rmat[2,2]+rmat[3,3]))]) 
  } else{
    conc.U0.U12 <- qr.Q(qr(cbind(U0, U12))) ; U12 <- as.matrix(conc.U0.U12[,-c(1:r0)])
    conc.U0.U13 <- qr.Q(qr(cbind(U0, U13))) ; U13 <- as.matrix(conc.U0.U13[,-c(1:r0)])
    conc.U0.U23 <- qr.Q(qr(cbind(U0, U23))) ; U23 <- as.matrix(conc.U0.U23[,-c(1:r0)])
  }
  
  U <- cbind(U0, U12, U13, U23, U1, U2, U3)
  
  # generate dataset-specific loadings
  Vtmp <- list()
  for (d in 1:D){
    Vtmp[[d]] <- matrix(runif(pvec[d]*(r0+sum(rmat[d,]))), pvec[d],(r0+sum(rmat[d,])))
    # orthogonalize the loadings
    Vtmp[[d]] <- qr.Q(qr(Vtmp[[d]]))
  }
  
  # Put everything together into one V
  V <- matrix(0, nrow = p_total, ncol = r_total)
  
  for (d in 1:D){
    
    # Row index
    if (d == 1){
      index <- c(1:pvec[1])
    } else{
      index <- c((pcum[d-1]+1):pcum[d])
    }
    
    # Column index, order is r0, r12, r13, r23, r1, r2, r3
    if (d==1){
      col_indexd = c(1:(r0+rmat[1,2]+rmat[1,3]),(r0+rmat[1,2]+rmat[1,3]+rmat[2,3]+1):(r0+rmat[1,2]+rmat[1,3]+rmat[2,3]+rmat[1,1]))
    } else if (d==2){
      col_indexd = c(1:(r0+rmat[1,2]),(r0+rmat[1,2]+rmat[1,3]+1):(r0+rmat[1,2]+rmat[1,3]+rmat[2,3]), (r0+rmat[1,2]+rmat[1,3]+rmat[2,3]+rmat[1,1]+1):(r0+rmat[1,2]+rmat[1,3]+rmat[2,3]+rmat[1,1]+rmat[2,2]))
    } else{
      col_indexd = c(1:r0,(r0+rmat[1,2]+1):(r0+rmat[1,2]+rmat[1,3]+rmat[2,3]),(r0+rmat[1,2]+rmat[1,3]+rmat[2,3]+rmat[1,1]+rmat[2,2]+1):(r0+rmat[1,2]+rmat[1,3]+rmat[2,3]+rmat[1,1]+rmat[2,2]+rmat[3,3]))
    }
    V[index,col_indexd]=Vtmp[[d]]
  }
  
  # Orthogonallize V
  if (orthogonalV == T){
    V <- qr.Q(qr(V))
  }
  # Add the scale between joint, semi-joint and individual by fixing singular values
  # Column index, order is r0, r12, r13, r23, r1, r2, r3
  V <- V%*%diag(c(sjoint,s12, s13, s23, s1ind, s2ind, s3ind))
  
  # Form X
  X <- tcrossprod(U, V) ; true_M <- X
  
  # Datasets have different scales or different dimensions
  sigmavec = rep(0, D) ; err_list <- list()
  
  for (d in 1:D){
    # Row index for dth dataset
    if (d == 1){
      index <- c(1:pvec[1])
    } else{
      index <- c((pcum[d-1]+1):pcum[d])
    }
    sigmavec[d] <- sqrt(sum(X[,index]^2) / (n * pvec[d] * snr_vec[d]))
    err_list[[d]] <- matrix(rnorm(n * pvec[d], sd = sigmavec[d]), n, pvec[d])
    X[,index] <- X[,index] + err_list[[d]]
  }
  
  # Combine joint loadings and individual loadings into one matrix V
  return(list(X = X, pvec = pvec, sigmavec = sigmavec, U = U, V = V, snr = snr, true_M = true_M, 
              err_list = err_list)) 
}


dgen_d3_ver2 <- function(n, pvec, sjoint, s12, s13, s23, s1ind, s2ind, s3ind, snr_vec, scaleU = T, orthogonalV = T){
  
  p_total <- sum(pvec) ; pcum <- cumsum(pvec) ; cvec <- c(1,1,1) ; D <- 3 
  r0 <- length(sjoint) ; rmat <- matrix(NA, nrow = 3, ncol = 3)
  
  diag(rmat) <- c(length(s1ind), length(s2ind), length(s3ind))
  rmat[1,2] <- rmat[2,1] <- length(s12)
  rmat[1,3] <- rmat[3,1] <- length(s13)
  rmat[2,3] <- rmat[3,2] <- length(s23)
  
  r_total <- r0 + sum(rmat[upper.tri(rmat, diag=T)])
  
  # generate noise part of X
  X <- matrix(0, n, p_total)
  
  # generate all the scores
  U0 <- matrix(runif(n*r0, 0, 1), n, r0) 
  U12 <- matrix(runif(n*rmat[1,2], 0, 1), n, rmat[1,2]) 
  U13 <- matrix(runif(n*rmat[1,3], 0, 1), n, rmat[1,3]) 
  U23 <- matrix(runif(n*rmat[2,3], 0, 1), n, rmat[2,3]) 
  U1 <- matrix(runif(n*rmat[1,1], 0, 1), n, rmat[1,1]) 
  U2 <- matrix(runif(n*rmat[2,2], 0, 1), n, rmat[2,2]) 
  U3 <- matrix(runif(n*rmat[3,3], 0, 1), n, rmat[3,3]) 
  
  # Column-center U so that the overall mean is zero
  if(scaleU){
    U0 <- scale(U0, scale = F) ; U12 <- scale(U12, scale = F) ; U13 <- scale(U13, scale = F) ; U23 <- scale(U23, scale = F) 
    U1 <- scale(U1, scale = F) ; U2 <- scale(U2, scale = F) ; U3 <- scale(U3, scale = F)
  }  
  
  # orthogonality between 
  # 1. globally shared scores and individual scores
  # 2. globally shared scores and pairwisely shared scores
  U0 <- as.matrix(qr.Q(qr(U0)))
  conc_U0_U12 <- qr.Q(qr(cbind(U0, U12))) ; U12 <- as.matrix(conc_U0_U12[,-c(1:r0)])
  conc_U0_U13 <- qr.Q(qr(cbind(U0, U13))) ; U13 <- as.matrix(conc_U0_U13[,-c(1:r0)])
  conc_U0_U23 <- qr.Q(qr(cbind(U0, U23))) ; U23 <- as.matrix(conc_U0_U23[,-c(1:r0)])
  
  mat1 <- cbind(U0, U12, U13) ; U1 <- qr.Q(qr((diag(1,n) - mat1%*%solve(crossprod(mat1))%*%t(mat1))%*%U1))
  mat2 <- cbind(U0, U12, U23) ; U2 <- qr.Q(qr((diag(1,n) - mat2%*%solve(crossprod(mat2))%*%t(mat2))%*%U2))
  mat3 <- cbind(U0, U13, U23) ; U3 <- qr.Q(qr((diag(1,n) - mat3%*%solve(crossprod(mat3))%*%t(mat3))%*%U3))
  
  U <- cbind(U0, U12, U13, U23, U1, U2, U3)
  
  # generate dataset-specific loadings
  Vtmp <- list()
  for (d in 1:D){
    Vtmp[[d]] <- matrix(runif(pvec[d]*(r0+sum(rmat[d,]))), pvec[d],(r0+sum(rmat[d,])))
    # orthogonalize the loadings
    Vtmp[[d]] <- qr.Q(qr(Vtmp[[d]]))
  }
  
  # Put everything together into one V
  V <- matrix(0, nrow = p_total, ncol = r_total)
  
  for (d in 1:D){
    
    # Row index
    if (d == 1){
      index <- c(1:pvec[1])
    } else{
      index <- c((pcum[d-1]+1):pcum[d])
    }
    
    # Column index, order is r0, r12, r13, r23, r1, r2, r3
    if (d==1){
      col_indexd = c(1:(r0+rmat[1,2]+rmat[1,3]),(r0+rmat[1,2]+rmat[1,3]+rmat[2,3]+1):(r0+rmat[1,2]+rmat[1,3]+rmat[2,3]+rmat[1,1]))
    } else if (d==2){
      col_indexd = c(1:(r0+rmat[1,2]),(r0+rmat[1,2]+rmat[1,3]+1):(r0+rmat[1,2]+rmat[1,3]+rmat[2,3]), (r0+rmat[1,2]+rmat[1,3]+rmat[2,3]+rmat[1,1]+1):(r0+rmat[1,2]+rmat[1,3]+rmat[2,3]+rmat[1,1]+rmat[2,2]))
    } else{
      col_indexd = c(1:r0,(r0+rmat[1,2]+1):(r0+rmat[1,2]+rmat[1,3]+rmat[2,3]),(r0+rmat[1,2]+rmat[1,3]+rmat[2,3]+rmat[1,1]+rmat[2,2]+1):(r0+rmat[1,2]+rmat[1,3]+rmat[2,3]+rmat[1,1]+rmat[2,2]+rmat[3,3]))
    }
    V[index,col_indexd]=Vtmp[[d]]
  }
  
  # Orthogonallize V
  if (orthogonalV == T){
    V <- qr.Q(qr(V))
  }
  # Add the scale between joint, semi-joint and individual by fixing singular values
  # Column index, order is r0, r12, r13, r23, r1, r2, r3
  V <- V%*%diag(c(sjoint,s12, s13, s23, s1ind, s2ind, s3ind))
  
  # Form X
  X <- tcrossprod(U, V) ; true_M <- X
  
  # Datasets have different scales or different dimensions
  sigmavec = rep(0, D) ; err_list <- list()
  
  for (d in 1:D){
    # Row index for dth dataset
    if (d == 1){
      index <- c(1:pvec[1])
    } else{
      index <- c((pcum[d-1]+1):pcum[d])
    }
    sigmavec[d] <- sqrt(sum(X[,index]^2) / (n * pvec[d] * snr_vec[d]))
    err_list[[d]] <- matrix(rnorm(n * pvec[d], sd = sigmavec[d]), n, pvec[d])
    X[,index] <- X[,index] + err_list[[d]]
  }
  
  # Combine joint loadings and individual loadings into one matrix V
  return(list(X = X, pvec = pvec, sigmavec = sigmavec, U = U, V = V, true_M = true_M, 
              err_list = err_list)) 
}



#############################################################################################################################


sim_data_2 <- function(n, pvec, sig_var, r, r1, r2){
  
  # joint
  U1 <- matrix(rnorm(r*pvec[1]),ncol=r,nrow=pvec[1])
  U2 <- matrix(rnorm(r*pvec[2]),ncol=r,nrow=pvec[2])
  V <- matrix(rnorm(r*n),ncol=r,nrow=n)
  
  # individual (View 1)
  U1.i <- matrix(rnorm(r1*pvec[1]),ncol=r1,nrow=pvec[1])
  V1.i <- matrix(rnorm(r1*n),ncol=r1,nrow=n)
  
  # individual (View 2)
  U2.i <- matrix(rnorm(r2*pvec[2]),ncol=r2,nrow=pvec[2])
  V2.i <- matrix(rnorm(r2*n),ncol=r2,nrow=n)
  
  J1 <- sig_var*U1%*%t(V) ; J2 <- sig_var*U2%*%t(V)  
  I1 <- sig_var*U1.i%*%t(V1.i); I2 <- sig_var*U2.i%*%t(V2.i) 
  
  S1 <- (J1+I1) ; S2 <- (J2+I2) 
  
  X1 <- S1 + matrix(rnorm(pvec[1]*n),nrow=pvec[1],ncol=n)
  X2 <- S2 + matrix(rnorm(pvec[2]*n),nrow=pvec[2],ncol=n)
  
  X0 <- rbind(X1,X2) ; S0 <- rbind(S1, S2)
  
  return(list(X0 = X0, S0 = S0))
  
}


sim_data_3 <- function(n, pvec, sig_var, r, r1, r2, r3, r12, r13, r23){
  
  # joint
  U1 <- matrix(rnorm(r*pvec[1]),ncol=r,nrow=pvec[1])
  U2 <- matrix(rnorm(r*pvec[2]),ncol=r,nrow=pvec[2])
  U3 <- matrix(rnorm(r*pvec[3]),ncol=r,nrow=pvec[3])
  V <- matrix(rnorm(r*n),ncol=r,nrow=n)
  
  # partial (View 1 and View 2)
  U1.i12 <- matrix(rnorm(r12*pvec[1]),ncol=r12,nrow=pvec[1])
  U2.i12 <- matrix(rnorm(r12*pvec[2]),ncol=r12,nrow=pvec[2])
  V.i12 <- matrix(rnorm(r12*n),ncol=r12,nrow=n)
  
  # partial (View 1 and View 3)
  U1.i13 <- matrix(rnorm(r13*pvec[1]),ncol=r13,nrow=pvec[1])
  U3.i13 <- matrix(rnorm(r13*pvec[3]),ncol=r13,nrow=pvec[3])
  V.i13 <- matrix(rnorm(r13*n),ncol=r13,nrow=n)
  
  # partial (View 2 and View 3)
  U2.i23 <- matrix(rnorm(r23*pvec[2]),ncol=r23,nrow=pvec[2])
  U3.i23 <- matrix(rnorm(r23*pvec[3]),ncol=r23,nrow=pvec[3])
  V.i23 <- matrix(rnorm(r23*n),ncol=r23,nrow=n)
  
  # individual (View 1)
  U1.i <- matrix(rnorm(r1*pvec[1]),ncol=r1,nrow=pvec[1])
  V1.i <- matrix(rnorm(r1*n),ncol=r1,nrow=n)
  
  # individual (View 2)
  U2.i <- matrix(rnorm(r2*pvec[2]),ncol=r2,nrow=pvec[2])
  V2.i <- matrix(rnorm(r2*n),ncol=r2,nrow=n)
  
  # individual (View 3)
  U3.i <- matrix(rnorm(r3*pvec[3]),ncol=r3,nrow=pvec[3])
  V3.i <- matrix(rnorm(r3*n),ncol=r3,nrow=n)
  
  J1 <- sig_var*U1%*%t(V) ; J2 <- sig_var*U2%*%t(V) ; J3 <- sig_var*U3%*%t(V)
  I1 <- sig_var*U1.i%*%t(V1.i); I2 <- sig_var*U2.i%*%t(V2.i) ; I3 <- sig_var*U3.i%*%t(V3.i)
  I1.12 <- sig_var*U1.i12%*%t(V.i12) ; I2.12 <- sig_var*U2.i12%*%t(V.i12)
  I1.13 <- sig_var*U1.i13%*%t(V.i13) ; I3.13 <- sig_var*U3.i13%*%t(V.i13)
  I2.23 <- sig_var*U2.i23%*%t(V.i23) ; I3.23 <- sig_var*U3.i23%*%t(V.i23)
  
  S1 <- (J1+I1+I1.12+I1.13) ; S2 <- (J2+I2+I2.12+I2.23) ; S3 <- (J3+I3+I3.13+I3.23)
  
  X1 <- S1 + matrix(rnorm(pvec[1]*n),nrow=pvec[1],ncol=n)
  X2 <- S2 + matrix(rnorm(pvec[2]*n),nrow=pvec[2],ncol=n)
  X3 <- S3 + matrix(rnorm(pvec[3]*n),nrow=pvec[3],ncol=n)
  
  X0 <- rbind(X1,X2,X3) ; S0 <- rbind(S1, S2, S3)
  
  return(list(X0 = X0, S0 = S0))
  
}