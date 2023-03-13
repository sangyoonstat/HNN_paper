
require(SLIDE) # https://github.com/irinagain/SLIDE

################################################################################################################
######## The below code is the modified version of the code from https://github.com/irinagain/SLIDE ############
################################################################################################################

slide_BCV <- function(X, pvec, structure_list, fold_id_n, fold_id_p, k_max = 1000,
                      eps = 1e-06, center = T) {
  d <- length(pvec)
  n <- nrow(X)
  p <- ncol(X)
  pcum <- c(0, cumsum(pvec))
  n_fold <- length(unique(fold_id_n)) ; p_fold <- length(unique(fold_id_p))
  # Get length of structure list
  n_structure = length(structure_list)
  
  # Get the folds For samples
  #fold_id_n <- sample(rep(seq_len(n_fold), length.out = n))
  # For measurements
  #fold_id_p <- c()
  #for (i in 1:d) {
  #  fold_id_p <- c(fold_id_p, sample(rep(seq_len(p_fold), length.out = pvec[i])))
  #}
  
  # Since bcv submatrices have lower dimensions than original ones, need to ensure that the tried ranks (number of columns in structure_list) do not exceed the dimensions of subfolds
  
  # Save prediction error
  error <- matrix(0, n_fold * p_fold, n_structure)
  iter <- 0
  for (k in 1:n_fold) {
    for (j in 1:p_fold) {
      iter <- iter + 1
      # Update the pvec
      pvec_foldj = pvec
      for (i in 1:d) {
        index <- c((pcum[i] + 1):pcum[i + 1])
        pvec_foldj[i] = sum(fold_id_p[index] != j)
      }
      pcum_fold = c(0, cumsum(pvec_foldj))
      # Standardize each dataset
      out_s <- standardizeX(X = X[fold_id_n != k, fold_id_p != j], pvec = pvec_foldj,
                            center = center)
      Xfold <- out_s$X
      
      # Consider all structure from the list, fit the model, evaluate BCV error
      for (l in 1:n_structure) {
        # Find U and V based on the given structure
        out <- slide_givenS(X = Xfold, pvec = pvec_foldj, S = structure_list[[l]],
                            k_max = k_max, eps = eps)
        Vadj <- out$V
        for (i in 1:d) {
          index <- c((pcum_fold[i] + 1):pcum_fold[i + 1])
          Vadj[index, ] <- out$V[index, ] * sqrt(out_s$norms[i])
        }
        # Perform SVD on the output
        out_svd <- svd(tcrossprod(out$U, Vadj))
        if (max(out_svd$d) < eps) {
          error[iter, l] <- d
        } else {
          # because of est_givenranks, only have nonzero columns already
          if (center == F) {
            # No centering was done
            tildeV = crossprod(X[fold_id_n != k, fold_id_p == j], out$U)
            tildeU = X[fold_id_n == k, fold_id_p != j] %*% Vadj %*% solve(crossprod(Vadj))
            # Calculate the error
            pcum_notfold <- pcum - pcum_fold
            for (i in 1:d) {
              index <- c((pcum_notfold[i] + 1):pcum_notfold[i + 1])
              error[iter, l] = error[iter, l] + sum((X[fold_id_n == k, which(fold_id_p == j)[index]] - tcrossprod(tildeU, tildeV[index, ]))^2)/sum(X[fold_id_n ==k, which(fold_id_p == j)[index]]^2)
            }
          } else {
            n_newf = sum(fold_id_n != k)  # number of samples in the submatrix for model fitting
            # Centering was done, stored in out_s$Xmean
            tildeV = cbind(colMeans(X[fold_id_n != k, fold_id_p == j]) * sqrt(n_newf),
                           crossprod(X[fold_id_n != k, fold_id_p == j], out$U))
            Vnew = cbind(sqrt(n_newf) * out_s$Xmean, Vadj)
            # Adjust for potential singularity here
            Vsvd = svd(Vnew)
            Vproj = Vsvd$u[, Vsvd$d > eps] %*% diag(1/Vsvd$d[Vsvd$d > eps]) %*%
              t(Vsvd$v[, Vsvd$d > eps])
            tildeU = X[fold_id_n == k, fold_id_p != j] %*% Vproj
            # Calculate the error
            pcum_notfold <- pcum - pcum_fold
            for (i in 1:d) {
              index <- c((pcum_notfold[i] + 1):pcum_notfold[i + 1])
              error[iter, l] = error[iter, l] + sum((X[fold_id_n == k, which(fold_id_p == j)[index]] - tcrossprod(tildeU, tildeV[index, ]))^2)/sum(scale(X[fold_id_n == k, which(fold_id_p == j)[index]], scale = F)^2)
            }
          }
        }
      }  # end for lamda_seq
    }  # end for p folds
  }  # end for n folds
  # Calculate average prediction error for each tuning parameter
  error_sum <- colSums(error)
  # Find largest lambda corresponding to minimal average error, here assume that
  # lambda_seq is sorted from largest to smallest
  id_min <- min(c(1:n_structure)[error_sum == min(error_sum)])
  structure_min <- structure_list[[id_min]]
  return(list(structure_list = structure_list, structure_min = structure_min, error_sum = error_sum,
              error = error, fold_id_n = fold_id_n, fold_id_p = fold_id_p))
}


slide <- function(X, pvec, n_lambda = 50, lambda_min = 0.01, fold_id_n, fold_id_p, center = T, k_max = 5000, eps = 1e-06, ratio_max = NULL) {
  
  n_fold <- length(unique(fold_id_n)) ; p_fold <- length(unique(fold_id_p))
  
  # Center and scale the original dataset
  out_s <- standardizeX(X, pvec, center = center)
  
  # Make restriction on the maximal allowable rank as a ratio of n, p in each dataset
  if (is.null(ratio_max)){
    ratio_max = min((n_fold - 1)/n_fold, (p_fold - 1)/p_fold)
  }
  
  # Form the list of candidate structures
  out_struct <- create_structure_list(X = out_s$X, pvec = pvec, lambda_max = max(out_s$svec), standardized = T, lambda_min = lambda_min, n_lambda = n_lambda, ratio_max = ratio_max)
  
  # Select the structure from the list using bcv procedure
  out_bcv <- slide_BCV(X = out_s$X, pvec = pvec, structure_list = out_struct$Slist, fold_id_n = fold_id_n,
                       fold_id_p = fold_id_p, k_max = k_max, eps = eps, center = center)
  
  # Fit slide model on a selected structure S
  out_slide <- slide_givenS(X = out_s$X, pvec = pvec, S = out_bcv$structure_min, k_max = k_max,
                            eps = eps, standardized = T)
  
  return(list(out_s = out_s, structure_list = out_struct$Slist, S = out_bcv$structure_min, model = out_slide,
              fold_id_n = out_bcv$fold_id_n, fold_id_p = out_bcv$fold_id_p, bcv_err = out_bcv$error_sum))
}
