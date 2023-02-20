
###############################################################################################################
######## The below code is the modified version of the code from https://github.com/lockEF/bidifac ############
###############################################################################################################


require(rARPACK)
require(denoiseR)

frob_sq=function(X){ sum(X^2, na.rm=T) }

softSVD=function(X, lambda){
  svdX=svd(X)
  nuc=pmax(svdX$d-lambda,0)
  out=tcrossprod(svdX$u, tcrossprod( svdX$v,diag(nuc) ))
  return(list(out=out, nuc=sum(nuc)))
}

prep_dat <- function(X_list, center = T){
  
  D <- length(X_list) ; sigma_vec <- pvec <- c() 

  for(d in 1:D){
    
    pvec[d] <- ncol(X_list[[d]])
    
    if(center){
      col_means <- colMeans(X_list[[d]])
      mean_mat <- matrix(col_means, nrow = nrow(X_list[[d]]), ncol = ncol(X_list[[d]]), byrow = T)
      centered_X <- X_list[[d]] - mean_mat 
    } else{
      mean_mat <- matrix(0, nrow = nrow(X_list[[d]]), ncol = ncol(X_list[[d]]))
      centered_X <- X_list[[d]] - mean_mat 
    }

    sigma_vec[d] <- estim_sigma(centered_X, method = "MAD", center = F) ; Zmat <- centered_X/sigma_vec[d] 
    
    if(d==1){
      out_mean <- mean_mat ; out_Z <- Zmat 
    } else{
      out_mean <- cbind(out_mean, mean_mat) ; out_Z <- cbind(out_Z, Zmat)
    }
    
  }
  
  return(list(Z = out_Z, mean_mat = out_mean, sigma_vec = sigma_vec, n = nrow(X_list[[1]]), pvec = pvec))
  
}


######################################################################################################################################################

my_UNIFAC <- function(X_list, center = F, eps = 10^-3, max_iter = 1000, seednum = NULL){
  
  if (!is.null(seednum)){set.seed(seednum)}
  
  # Center and scale the data by MAD
  pdat <- prep_dat(X_list, center)
  sigma_vec <- pdat$sigma_vec ; mean_mat <- pdat$mean_mat
  X00 <- pdat$Z ; n <- pdat$n ; pvec <- pdat$pvec ; D <- length(pvec)
  
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
    S00_list[[j]] <- R00_list[[j]] + I00_list[[j]] + mean_mat[,ind2]
  }
  
  
  return(list(S = S00_list, R = R00_list, I = I00_list, mean_mat = mean_mat, X00 = X00, sigma_vec = sigma_vec))
  
}

######################################################################################################################################################

unifac.plus <- function(X0,p.ind,max.comb=20,num.comp=20,max.iter=500,conv.thresh=0.001,temp.iter=100){
  n.source <- length(p.ind)
  X0.resid <- X0
  S <- list()
  pen <- c()
  for(i in 1:max.comb){
    S[[i]]=0
    pen[i]=0}
  obj.vec <- c(sum(X0^2))
  #obj.vec <- c()
  temp.fac <- svd(X0)$d[1]/(sum(sqrt(dim(X0))))-1
  obj.prev <- 10^10
  for(jj in 1:max.iter){
    #print(jj)
    if(jj < temp.iter){
      lambda <- 1+(temp.iter-1-jj)/temp.iter*temp.fac
    }
    p.ind.list <- list()
    for(k in 1:max.comb){
      temp.source <- c(1:n.source)
      X0.temp <- X0.resid+S[[k]]
      X0.c <- tcrossprod(X0.temp)
      X0.rl <- list()
      for(i in 1:n.source) X0.rl[[i]] <- crossprod(X0.temp[p.ind[[i]],])
      cur.p.ind <- c()
      sse.pen.tot <- c()
      SSE.pen.temp <- c()
      cur.p.ind <- c()
      cur.p.ind.list <- list()
      cur_cross=0
      for(i in 1:n.source){
        SSE.pen.temp <- c()
        for(j in temp.source){
          temp.p.ind <- c(cur.p.ind,p.ind[[j]]) 
          pl <- length(temp.p.ind)
          if(pl<=n){
            a <- sqrt(eigs_sym(X0.c[temp.p.ind,temp.p.ind], num.comp, which = "LM",retvec=FALSE)$values)
          }  
          if(pl>n){
            a <- sqrt(eigs_sym(cur_cross+X0.rl[[j]], num.comp, which = "LM",retvec=FALSE)$values)
          }  
          SSE.pen.temp[j] <- sum(pmax(a-lambda*(sqrt(pl)+sqrt(n)),0)^2)
        }    
        indmin <- which.max(SSE.pen.temp)
        sse.pen.tot[i] <- SSE.pen.temp[indmin]
        cur.p.ind <- c(cur.p.ind,p.ind[[indmin]])
        cur.p.ind.list[[i]] <- cur.p.ind
        cur_cross <- cur_cross+X0.rl[[indmin]]
        temp.source <- temp.source[temp.source!=indmin]
      }
      Order = order(sse.pen.tot, decreasing=TRUE)
      for(i in 1:length(sse.pen.tot)){
        cur.p.ind <- cur.p.ind.list[[Order[i]]]
        CurPen <- sse.pen.tot[Order[i]]
        if(!(list(sort(cur.p.ind))%in%p.ind.list)) break
        #print('GGG!')
      }
      if(CurPen<=0){
        pen = pen[-k]
        S[[k]] <- NULL
        pen[max.comb]=0
        S[[max.comb]]=0
        X0.resid=X0.temp
        break
      }
      p.ind.list[[k]] <- sort(cur.p.ind)
      a <- sqrt(eigs_sym(X0.c[cur.p.ind,cur.p.ind], num.comp, which = "LM",retvec=FALSE)$values)
      nc <- sum(a>(lambda*(sqrt(length(cur.p.ind))+sqrt(n))))
      SVD <- svd(X0.temp[cur.p.ind,], nu=nc,nv=nc)
      s.vals <- pmax(SVD$d[1:nc]-lambda*(sqrt(length(cur.p.ind))+sqrt(n)),0)
      Diag <- diag(s.vals,nrow=nc,ncol=nc)
      Est <- SVD$u%*%Diag%*%t(SVD$v)
      S[[k]] <- array(rep(0,prod(dim(X0))),dim=dim(X0))
      S[[k]][cur.p.ind,] <- Est
      X0.resid[cur.p.ind,] <- X0.temp[cur.p.ind,]-Est
      X0.resid = X0-Reduce('+',S)
      pen[k] <- lambda*(sqrt(n)+sqrt(length(cur.p.ind)))*(sum(s.vals))
      obj.vec <- c(obj.vec,sum(X0.resid^2)+2*sum(pen))
      #     if(abs(obj.vec[jj+1]-obj.vec[jj])<0.001) break
    }
    obj.cur <- sum(X0.resid^2)+2*sum(pen)
    if(abs(obj.cur-obj.prev)<conv.thresh) break
    obj.prev <- sum(X0.resid^2)+2*sum(pen)
  }
  Sums <- matrix(nrow=k-1,ncol=n.source)
  for(kk in 1:(k-1)){for(j in 1:n.source){
    Sums[kk,j] = sum(norm(S[[kk]][p.ind[[j]],])^2)
  }}
  return(list(S=S,p.ind.list=p.ind.list,Sums=Sums,obj.vec=obj.vec))
}

######################################################################################################################################################


unifac.plus.given <- function(X0,p.ind,p.ind.list,num.comp=20,max.iter=500,conv.thresh=0.001){
  n.source <- length(p.ind)
  X0.resid <- X0
  S <- list()
  pen <- c()
  max.comb <- length(p.ind.list)
  for(i in 1:max.comb){
    S[[i]]=0
    pen[i]=0}
  obj.vec <- c(sum(X0^2))
  #obj.vec <- c()
  temp.fac <- svd(X0)$d[1]/(sum(sqrt(dim(X0))))-1
  obj.prev <- 10^10
  for(jj in 1:max.iter){
    #print(jj)
    if(jj < 100){
      lambda <- 1+(99-jj)/100*temp.fac
    }
    for(k in 1:max.comb){
      X0.temp <- X0.resid+S[[k]]
      X0.r <- crossprod(X0.temp[p.ind.list[[k]],])
      #   print(paste('*',k))
      a <- sqrt(eigs_sym(X0.r, num.comp, which = "LM",retvec=FALSE)$values)
      nc <- sum(a>(lambda*(sqrt(length(p.ind.list[[k]]))+sqrt(n))))
      S[[k]] <- array(rep(0,prod(dim(X0))),dim=dim(X0))
      pen[k] = 0
      #     Est =
      if(nc>0){
        SVD <- svd(X0.temp[p.ind.list[[k]],], nu=nc,nv=nc)
        s.vals <- pmax(SVD$d[1:nc]-lambda*(sqrt(length(p.ind.list[[k]]))+sqrt(n)),0)
        Diag <- diag(s.vals,nrow=nc,ncol=nc)
        Est <- SVD$u%*%Diag%*%t(SVD$v)
        S[[k]] <- array(rep(0,prod(dim(X0))),dim=dim(X0))
        S[[k]][p.ind.list[[k]],] <- Est
        pen[k] <- lambda*(sqrt(n)+sqrt(length(p.ind.list[[k]])))*(sum(s.vals))
      }
      X0.resid = X0-Reduce('+',S)
      obj.vec <- c(obj.vec,sum(X0.resid^2)+2*sum(pen))
      #     if(abs(obj.vec[jj+1]-obj.vec[jj])<0.001) break
    }
    obj.cur <- sum(X0.resid^2)+2*sum(pen)
    if(abs(obj.cur-obj.prev)<conv.thresh) break
    obj.prev <- sum(X0.resid^2)+2*sum(pen)
  }
  
  Sums <- matrix(nrow=max.comb,ncol=n.source)
  for(kk in 1:max.comb){for(j in 1:n.source){
    Sums[kk,j] = sum(norm(S[[kk]][p.ind[[j]],])^2)
  }}
  return(list(S=S,p.ind.list=p.ind.list,Sums=Sums,obj.vec=obj.vec))
}
