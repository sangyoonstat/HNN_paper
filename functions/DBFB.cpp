#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;

/******************************************************************************************************************
 ****** 1. Do truncated SVD where singular values smaller than specified tolerance level are set to be zero *******
 ******************************************************************************************************************/

/************************************************************************
 ****** 1.1 Print singular value larger than given tolerance only *******
 ************************************************************************/

// [[Rcpp::export]]
arma::vec sgv_single(const arma::mat& X, double tol){
  
  arma::vec sgv;
  
  // SVD of X (Question: svd_econ or svd?)
  arma::svd(sgv, X);
  
  // Set singular values smaller than cutoff to be 0
  arma::uvec ind = arma::find(sgv>tol);
  
  // Return the resulting matrix after soft-thresholding singular values   
  return sgv(ind);
  
}

/*****************************************************************************
 ****** 1.2 Print left singular matrix and its rank from truncated SVD *******
 *****************************************************************************/

// [[Rcpp::export]]
Rcpp::List tSVD_left(const arma::mat& X, double tol){
  
  arma::mat U, V;
  arma::vec sgv;
  
  // SVD of X (Question: svd_econ or svd?)
  arma::svd(U, sgv, V, X, "std");
  arma::uvec ind = arma::find(sgv>tol);
  
  double rtilde = ind.n_elem;
  arma::mat Utilde = U.cols(ind);
  
  return Rcpp::List::create(Rcpp::Named("u") = Utilde, Rcpp::Named("r") = rtilde);
  
}

/***********************************************************************************************
 ****** 1.3 Print singular values as well as left and right singular matrices separately *******
 ***********************************************************************************************/

// [[Rcpp::export]]
Rcpp::List full_tSVD(const arma::mat& X, double tol){
  
  arma::mat U, V;
  arma::vec sgv;
  
  // SVD of X (Question: svd_econ or svd?)
  arma::svd(U, sgv, V, X, "std");
  
  // Set singular values smaller than cutoff to be 0
  arma::uvec ind = arma::find(sgv>tol);
  
  // Return the resulting matrix after soft-thresholding singular values   
  return Rcpp::List::create(Rcpp::Named("u") = U, Rcpp::Named("v") = V, Rcpp::Named("d") = sgv);
  
}

/*************************************************************************************
 ****** 2. function to calculate matrix after singular value soft-thresholding *******
 *************************************************************************************/

// [[Rcpp::export]]
arma::mat rcpp_sft(const arma::mat& X, double cutoff){
  
  arma::mat U, V;
  arma::vec sgv;
  
  // SVD of X 
  arma::svd(U, sgv, V, X, "std");
  
  // Set singular values smaller than cutoff to be 0
  // arma::mat D(r, r, arma::fill::zeros);
  arma::uvec ind = arma::find(sgv>cutoff);
  
  // Return the resulting matrix after soft-thresholding singular values   
  return U.cols(ind)*diagmat(sgv(ind)-cutoff)*V.cols(ind).t();
  
}

/************************************************************************************************
 ****** 3. Function to update dual variables in the dual block forward-backward algorithm *******
 ************************************************************************************************/

/*********************************************************
 ****** 3.1 function to update the dual variables  *******
 *********************************************************/  

// [[Rcpp::export]]
arma::mat update_dual(const arma::mat& X, double cutoff){
  return X - rcpp_sft(X, cutoff);
}  

/*********************************************************************************************
 ****** 3.2 function to implement the dual block forward-backward algorithm for D = 2  *******
 *********************************************************************************************/  

// [[Rcpp::export]]
Rcpp::List DBFB_d2(const arma::mat& Z, arma::vec pvec, double tau, double lambda, arma::vec ratio_indiv, 
                   double gamma, int max_iter, double eps){
  
  // Configure dimensions 
  int n = Z.n_rows;
  int p1 = pvec[0];
  int p2 = pvec[1];
 
 arma::mat upd_D1, upd_D2, upd_D12, gam_M;
 arma::mat upd_M(n, p1 + p2, arma::fill::zeros);
 
 // Create vectors for indexing (we can probably require these as input arguments?)
 arma::uvec ind1 = arma::regspace<arma::uvec>(0, p1-1);
 arma::uvec ind2 = arma::regspace<arma::uvec>(p1, p1+p2-1);
 
 // Initialize dual variables with zero matrices
 arma::mat curr_D1(n, p1, arma::fill::zeros);
 arma::mat curr_D2(n, p2, arma::fill::zeros);
 arma::mat curr_D12(n, p1+p2, arma::fill::zeros);

 // Calculate the tuning parameters by the ratio formulation  
 double lambda1 = tau*ratio_indiv[0];
 double lambda2 = tau*ratio_indiv[1];

 // Some auxiliary objects used in obj
 int conv_indc = 0;
 int num_iter = -1;
 Rcpp::NumericVector diff_val(max_iter, NA_REAL);
 
 // Initial M is simply the concatenated data
 arma::mat curr_M = Z;
 
 // Run the dual-block forward-backward algorithm until convergence  
 
 while(conv_indc<1){
   
   num_iter = num_iter + 1;
   
   gam_M = gamma*curr_M;
   
   // Step 1 : Update the dual variables
   
   // Individuals
   upd_D1 = update_dual(curr_D1 + gam_M.cols(ind1), lambda1);
   upd_D2 = update_dual(curr_D2 + gam_M.cols(ind2), lambda2);

   // Concatenated
   upd_D12 = update_dual(curr_D12 + gam_M, lambda);

   // Step 2 : Update the primal variables
   upd_M.cols(ind1) = curr_M.cols(ind1) - (upd_D1 - curr_D1) - (upd_D12.cols(ind1) - curr_D12.cols(ind1));
   upd_M.cols(ind2) = curr_M.cols(ind2) - (upd_D2 - curr_D2) - (upd_D12.cols(ind2) - curr_D12.cols(ind2));

   // decide convergence of algorithm based on scaled Frobenius norm of upd_M-curr_M
   // (or calculate the convergence criterion which Irina used?)
   diff_val[num_iter] = norm(upd_M-curr_M, "fro");
   
   if(diff_val[num_iter]<=eps){
     conv_indc = 1;
     break;
   } else{
     if(num_iter==max_iter){
       // Stop the algorithm once we arrive at the maximum number of iteration        
       break;
     } else{
       // Get ready for the next iteration
       curr_M = upd_M; 
       curr_D1 = upd_D1; 
       curr_D2 = upd_D2;
       curr_D12 = upd_D12; 
     }
   }
   
 }
 
 return Rcpp::List::create(Rcpp::Named("Mtilde") = upd_M, Rcpp::Named("diff_val") = diff_val, 
                           Rcpp::Named("num_iter") = num_iter, Rcpp::Named("conv_indc") = conv_indc);
   
}



/*********************************************************************************************
 ****** 3.3 function to implement the dual block forward-backward algorithm for D = 3  *******
 *********************************************************************************************/  

// [[Rcpp::export]]
Rcpp::List DBFB_d3(const arma::mat& Z, arma::vec pvec, double tau, double kappa, double lambda,   
                   arma::vec ratio_pair, arma::vec ratio_indiv, double gamma, int max_iter, double eps){
  
  // Configure dimensions 
  int n = Z.n_rows;
  int p1 = pvec[0];
  int p2 = pvec[1];
  int p3 = pvec[2];
  
  arma::mat upd_D1, upd_D2, upd_D3, upd_D12, upd_D13, upd_D23, upd_D123, gam_M;
  arma::mat upd_M(n, p1 + p2 + p3, arma::fill::zeros);
  
  // Create vectors for indexing (we can probably require these as input arguments?)
  arma::uvec ind1 = arma::regspace<arma::uvec>(0, p1-1);
  arma::uvec ind2 = arma::regspace<arma::uvec>(p1, p1+p2-1);
  arma::uvec ind3 = arma::regspace<arma::uvec>(p1+p2, p1+p2+p3-1);
  arma::uvec ind12 = join_cols(ind1, ind2);
  arma::uvec ind13 = join_cols(ind1, ind3);
  arma::uvec ind23 = join_cols(ind2, ind3);
  arma::uvec ind13_second = arma::regspace<arma::uvec>(p1, p1+p3-1);
  arma::uvec ind23_first = arma::regspace<arma::uvec>(0, p2-1);
  arma::uvec ind23_second = arma::regspace<arma::uvec>(p2, p2+p3-1);
  
  // Initialize dual variables with zero matrices
  arma::mat curr_D1(n, p1, arma::fill::zeros);
  arma::mat curr_D2(n, p2, arma::fill::zeros);
  arma::mat curr_D3(n, p3, arma::fill::zeros);
  arma::mat curr_D12(n, p1+p2, arma::fill::zeros);
  arma::mat curr_D13(n, p1+p3, arma::fill::zeros);
  arma::mat curr_D23(n, p2+p3, arma::fill::zeros);
  arma::mat curr_D123(n, p1+p2+p3, arma::fill::zeros);  
  
  // Calculate the tuning parameters by the ratio formulation  
  double lambda12 = kappa*ratio_pair[0];
  double lambda13 = kappa*ratio_pair[1];
  double lambda23 = kappa*ratio_pair[2];
  
  double lambda1 = tau*ratio_indiv[0];
  double lambda2 = tau*ratio_indiv[1];
  double lambda3 = tau*ratio_indiv[2];
  
  // Some auxiliary objects used in obj
  int conv_indc = 0;
  int num_iter = -1;
  Rcpp::NumericVector diff_val(max_iter, NA_REAL);
  
  // Initial M is simply the concatenated data
  arma::mat curr_M = Z;
  
  // Run the dual-block forward-backward algorithm until convergence  
  
  while(conv_indc<1){
    
    num_iter = num_iter + 1;
    
    gam_M = gamma*curr_M;
    
    // Step 1 : Update the dual variables
    
    // Individuals
    upd_D1 = update_dual(curr_D1 + gam_M.cols(ind1), lambda1);
    upd_D2 = update_dual(curr_D2 + gam_M.cols(ind2), lambda2);
    upd_D3 = update_dual(curr_D3 + gam_M.cols(ind3), lambda3);
    
    // Pairwise
    upd_D12 = update_dual(curr_D12 + gam_M.cols(ind12), lambda12);
    upd_D13 = update_dual(curr_D13 + gam_M.cols(ind13), lambda13);
    upd_D23 = update_dual(curr_D23 + gam_M.cols(ind23), lambda23);
    
    // Concatenated
    upd_D123 = update_dual(curr_D123 + gam_M, lambda);
    
    // Step 2 : Update the primal variables
    
    upd_M.cols(ind1) = curr_M.cols(ind1) - (upd_D1 - curr_D1) - (upd_D12.cols(ind1) - curr_D12.cols(ind1)) - (upd_D13.cols(ind1) - curr_D13.cols(ind1)) - (upd_D123.cols(ind1) - curr_D123.cols(ind1));
    upd_M.cols(ind2) = curr_M.cols(ind2) - (upd_D2 - curr_D2) - (upd_D12.cols(ind2) - curr_D12.cols(ind2)) - (upd_D23.cols(ind23_first) - curr_D23.cols(ind23_first)) - (upd_D123.cols(ind2) - curr_D123.cols(ind2));
    upd_M.cols(ind3) = curr_M.cols(ind3) - (upd_D3 - curr_D3) - (upd_D13.cols(ind13_second) - curr_D13.cols(ind13_second)) - (upd_D23.cols(ind23_second) - curr_D23.cols(ind23_second)) - (upd_D123.cols(ind3) - curr_D123.cols(ind3));    
    
    // decide convergence of algorithm based on scaled Frobenius norm of upd_M-curr_M
    // (or calculate the convergence criterion which Irina used?)
    diff_val[num_iter] = norm(upd_M-curr_M, "fro");
    
    if(diff_val[num_iter]<=eps){
      conv_indc = 1;
      break;
    } else{
      if(num_iter==max_iter){
        // Stop the algorithm once we arrive at the maximum number of iteration        
        break;
      } else{
        // Get ready for the next iteration
        curr_M = upd_M; 
        curr_D1 = upd_D1; 
        curr_D2 = upd_D2;
        curr_D3 = upd_D3;
        curr_D12 = upd_D12; 
        curr_D13 = upd_D13;
        curr_D23 = upd_D23;
        curr_D123 = upd_D123;        
      }
    }
    
  }
  
  return Rcpp::List::create(Rcpp::Named("Mtilde") = upd_M, Rcpp::Named("diff_val") = diff_val, 
                            Rcpp::Named("num_iter") = num_iter, Rcpp::Named("conv_indc") = conv_indc);
  
}

