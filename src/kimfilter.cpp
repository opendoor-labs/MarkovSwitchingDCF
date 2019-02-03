#include <RcppArmadillo.h>
#include <Rcpp.h>

// Correctly setup the build environment
// [[Rcpp::depends(RcppArmadillo)]]

// Add a flag to enable OpenMP at compile time
// [[Rcpp::plugins(openmp)]]

// Protect against compilers without OpenMP
#ifdef _OPENMP
#include <omp.h>
#endif

// [[Rcpp::export]]
arma::mat gen_inv(arma::mat m){
  arma::mat out(m.n_rows, m.n_cols);
  try{
    out = inv(m);
  }catch(std::exception &ex){
    out = pinv(m);
  }
  return out;
}

// Finds the steady state probabilities from a transition matrix
// mat = |p_11 p_21 ... p_m1|
//       |p_12 p_22 ... p_m2|
//       |...            ...|
//       |p_1m p_2m ... p_mm|
// where the columns sum to 1
// [[Rcpp::export]]
arma::mat ss_prob(arma::mat mat){
  arma::mat Zero(mat.n_rows, 1);
  arma::mat One(1, mat.n_cols);
  One = One.ones();
  arma::mat A1(mat.n_rows, mat.n_rows);
  
  arma::mat A = join_cols(A1.eye() - mat, One);
  arma::mat B = gen_inv(A.t() * A) * A.t();
  return B * join_cols(Zero, One.submat(0, 0, 0, 0));;
}

// Likelihood function for the Kim filter
// [[Rcpp::export]]
arma::mat v_prob(arma::mat EV, arma::mat HE){
  return (1/sqrt(det(HE)))*exp(-0.5*EV.t() * gen_inv(HE) * EV); //Pr[Yt|S_t,Yt-1]
}

// [[Rcpp::export]]
Rcpp::List kim_filter(arma::cube B0, arma::cube P0, arma::cube Dt, arma::mat Ft, arma::mat Ht,
               arma::mat Qt, arma::mat Rt, arma::mat Tr_mat, arma::mat yt){
  
  //Define the storage lists
  int n_states = Tr_mat.n_cols;
  double tol = 1E-06;
  arma::cube B_tts(Ft.n_rows, yt.n_cols, n_states, arma::fill::zeros);
  arma::cube B_tlss(Ft.n_rows, yt.n_cols, n_states*n_states, arma::fill::zeros);
  arma::cube B_ttss(Ft.n_rows, yt.n_cols, n_states*n_states, arma::fill::zeros);
  arma::cube N_TLss(yt.n_rows, 1, n_states*n_states, arma::fill::zeros);
  arma::cube F_TLss(yt.n_rows, yt.n_rows, n_states*n_states, arma::fill::zeros);
  arma::cube Ks(Ft.n_rows, yt.n_rows, n_states, arma::fill::zeros);
  arma::cube Kss(Ft.n_rows, yt.n_rows, n_states*n_states, arma::fill::zeros);
  arma::cube K(Ft.n_rows, yt.n_rows, yt.n_cols, arma::fill::zeros);
  arma::cube PR_VLss(1, 1, n_states*n_states, arma::fill::zeros);
  arma::cube PROss(1, 1, n_states*n_states, arma::fill::zeros);
  
  arma::field<arma::cube> P_tts(n_states);
  for(int i = 0; i < n_states; i++){
    P_tts(i) = arma::cube(Ft.n_rows, Ft.n_cols, yt.n_cols, arma::fill::zeros);  
  }
  
  arma::field<arma::cube> P_tlss(n_states*n_states);
  arma::field<arma::cube> P_ttss(n_states*n_states);
  for(int i = 0; i < n_states*n_states; i++){
    P_tlss(i) = arma::cube(Ft.n_rows, Ft.n_cols, yt.n_cols, arma::fill::zeros);  
    P_ttss(i) = arma::cube(Ft.n_rows, Ft.n_cols, yt.n_cols, arma::fill::zeros);  
  }
  
  arma::mat Pr_tts(yt.n_cols, n_states, arma::fill::zeros);
  arma::mat Pr_tls(yt.n_cols, n_states, arma::fill::zeros);
  arma::mat Pr_tl(n_states, n_states, arma::fill::zeros);
  arma::mat D_tt(yt.n_cols, Ft.n_rows, arma::fill::zeros);
  arma::mat B_tt(yt.n_cols, Ft.n_rows, arma::fill::zeros);
  arma::uvec nonna_idx;
  arma::uvec na_idx;
  double lnl = 0.0;
  double Pr_val = 0.0;
  int s = 0;
  
  //Initialize the filter
  arma::cube B_lls = B0;
  arma::cube P_lls = P0;
  arma::mat Pr = ss_prob(Tr_mat);
    
  //Hamiliton + Kalman filter routine
  for(int i = 0; i < yt.n_cols; i++){
    Pr_tls.row(i) = (Tr_mat * Pr).t();

    // //Joint probabilities conditional on t-1
    if(Tr_mat.n_cols == 3){
      Pr_tl = Tr_mat % join_cols(join_cols(Pr.t(), Pr.t()), Pr.t()); //Pr[S_t=i,S_{t-1}=j|Y_{t-1}] = Pr[S_t=i|S_{t-1}=j,Y_{t-1}]*Pr[S_{t-1}=j|Y_{t-1}] 
    }else if(Tr_mat.n_cols == 2){
      Pr_tl = Tr_mat % join_cols(Pr.t(), Pr.t()); //Pr[S_t=i,S_{t-1}=j|Y_{t-1}] = Pr[S_t=i|S_{t-1}=j,Y_{t-1}]*Pr[S_{t-1}=j|Y_{t-1}]
    }else if(Tr_mat.n_cols == 1){
      Pr_tl = Tr_mat % Pr.t(); //Pr[S_t=i,S_{t-1}=j|Y_{t-1}] = Pr[S_t=i|S_{t-1}=j,Y_{t-1}]*Pr[S_{t-1}=j|Y_{t-1}]
    }
    
    s = 0;
    Pr_val = 0.0;
    for(int stl = 0; stl < n_states; stl++){
      for(int st = 0; st < n_states; st++){
        //When S_{t-1}=j, S_{t}=i
        //B^{i,j}_{t|t-1} = D_j + Fm %*% B^{j}_{t|t-1}
        B_tlss.slice(s).col(i) = Dt.slice(st) + Ft * B_lls.slice(stl);

        //Initial predictions of the unobserved componennt
        //P^{i,j}_{t|t-1} = P^{i,j}_{t|t-1} = Fm %*% P^{j}_{t-1|t-1} %*% t(Fm) + Qm
        P_tlss(s).slice(i) = Ft * P_lls.slice(stl) * Ft.t() + Qt;

        //Forecast errors
        //N^{i,j}_{t|t-1} = yti[, j] - Hm %*% B^{i,j}_{t|t-1}
        N_TLss.slice(s) = yt.col(i) - Ht * B_tlss.slice(s).col(i);
        //Ignore series with missing data
        nonna_idx = arma::find_finite(N_TLss.slice(s));
        na_idx = arma::find_nonfinite(N_TLss.slice(s));
        if(!na_idx.is_empty()){
          N_TLss.slice(s).rows(na_idx) = arma::mat(na_idx.n_elem, 1, arma::fill::zeros);
        }
        
        //Variance of forecast errors
        //F^{i,j}_{t|t-1} = Hm %*% P^{i,j}_{t|t-1} %*% t(Hm) + Rm
        F_TLss.slice(s) = Ht * P_tlss(s).slice(i) * Ht.t() + Rt;

        //Kalman gain
        //K^{i,j} = P^{i,j}_{t|t-1} %*% t(Hm) %*% solve(F^{i,j}_{t|t-1})
        Kss.slice(s) = P_tlss(s).slice(i) * Ht.t() * gen_inv(F_TLss.slice(s));

        //Updated predictions of unobserved component
        //B^{i,j}_{t|t} = B^{i,j}_{t|t-1} + K^{i,j} %*% N^{i,j}_{t|t-1}
        //Kalman gain from missing values is 0, same as if the model's prediction
        B_ttss.slice(s).col(i) = B_tlss.slice(s).col(i) + Kss.slice(s) * N_TLss.slice(s); 

        //Updated predictions of the covariance matrix of unobserved component
        //P^{i,j}_{t|t} = P^{i,j}_{t|t-1} - K^{i,j} %*% Hm %*% P^{i,j}_{t|t-1}
        P_ttss(s).slice(i) = P_tlss(s).slice(i) - Kss.slice(s) * Ht * P_tlss(s).slice(i);

        //f[y_t|S_t=j,S_{t-1}=i,Y_{t-1}]
        //PR_VL^{i,j} = v_prob(N^{i,j}_{t|t-1}, F^{i,j}_{t|t-1})*Pr_tl["d", "d"]
        //Ignore series with missing data
        // if(nonna_idx.is_empty()){ //If all series missing data
          PR_VLss.slice(s) = v_prob(N_TLss.slice(s), F_TLss.slice(s)) * Pr_tl(st, stl);
        // }else{
        //   PR_VLss.slice(s) = v_prob(N_TLss.slice(s).rows(nonna_idx), 
        //                 F_TLss.slice(s)(nonna_idx, nonna_idx)) * Pr_tl(st, stl);
        // }

        //Joint density conditional on t-1
        Pr_val = Pr_val + arma::as_scalar(PR_VLss.slice(s));
        s++;
      }
    }
    
    for(int s = 0; s < n_states*n_states; s++){
      //Joint probabilties conditional on t
      //Pr[S_t=j,S_{t-1}=i|Y_t] = f[y_t|S_t=j,S_{t-1}=i,Y_{t-1}]*Pr[S_t=j,S{t-1}=i|Yt]/f[y_t|Y_{t-1}]
      //PRO^{i,j} = PR_VL^{i,j}/Pr_val
      if(arma::as_scalar(Pr_val) > 0){
        PROss.slice(s) = PR_VLss.slice(s)/arma::as_scalar(Pr_val);  
      }else{
        PROss.slice(s) = PR_VLss.slice(s);
      }
    }

    for(int s = 0; s < n_states; s++){
      //Probabilities conditional on t: Pr[S_t=s|Yt]
      //Pr[j] = sum_{i}(PRO^{i,j})
      Pr.row(s) = 0;
      for(int j = 0; j < n_states; j++){
        Pr.row(s) += PROss.slice(s + j*n_states);
      }
    }
    Pr_tts.row(i) = Pr.t();

    //Collapsing terms
    K.slice(i) = arma::zeros(K.slice(i).n_rows, K.slice(i).n_cols);
    D_tt.row(i) = arma::zeros(D_tt.row(i).n_rows, D_tt.row(i).n_cols);
    B_tt.row(i) = arma::zeros(B_tt.row(i).n_rows, B_tt.row(i).n_cols);
    for(int s = 0; s < n_states; s++){
      B_tts.slice(s).col(i) = arma::zeros(B_tts.slice(s).col(i).n_rows, B_tts.slice(s).col(i).n_cols);
      for(int j = 0; j < n_states; j++){
        //B^{j}_{t|t} = (sum_{i}(PRO^{i,j}*B^{i,j}_{t|t}))/Pr[j]
        B_tts.slice(s).col(i) += arma::as_scalar(PROss.slice(s + j*n_states))*B_ttss.slice(s + j*n_states).col(i);
      }
      if(arma::as_scalar(Pr(s, 0)) < tol){
        B_tts.slice(s).col(i) /= arma::as_scalar(arma::datum::inf);
      }else{
        B_tts.slice(s).col(i) /= arma::as_scalar(Pr.row(s));
      }

      Ks.slice(s) = arma::zeros(Ks.slice(s).n_rows, Ks.slice(s).n_cols);
      for(int j = 0; j < n_states; j++){
        //K^{j} = (sum_{i}(PRO^{i,j}*K^{i,j}))/Pr[j]
        Ks.slice(s) += arma::as_scalar(PROss.slice(s + j*n_states))*Kss.slice(s + j*n_states);
      }
      if(arma::as_scalar(Pr(s, 0)) < tol){
        Ks.slice(s) /= arma::as_scalar(arma::datum::inf);
      }else{
        Ks.slice(s) /= arma::as_scalar(Pr.row(s));
      }

      P_tts(s).slice(i) = arma::zeros(P_tts(s).slice(i).n_rows, P_tts(s).slice(i).n_cols);
      for(int j = 0; j < n_states; j++){
        //P^{j}_{t|t} = (sum{i}(PRO^{i,j}*(P^{i,j}_{t|t} + (B^{j}_{t|t} - B^{i,j}_{t|t}) %*% t(B^{j}_{t|t} - B^{i,j}_{t|t}))))/Pr[j]
        arma::mat temp = B_tts.slice(s).col(i) - B_ttss.slice(s + j*n_states).col(i);
        P_tts(s).slice(i) += arma::as_scalar(PROss.slice(s + j*n_states))*(P_ttss(s + j*n_states).slice(i) + temp*temp.t());
      }
      if(arma::as_scalar(Pr(s, 0)) < tol){
        P_tts(s).slice(i) /= arma::as_scalar(arma::datum::inf);
      }else{
        P_tts(s).slice(i) /= arma::as_scalar(Pr.row(s));
      }

      
      K.slice(i) += arma::as_scalar(Pr.row(s)) * Ks.slice(s);
      D_tt.row(i) += (arma::as_scalar(Pr.row(s)) * Dt.slice(s)).t();
      B_tt.row(i) += (arma::as_scalar(Pr.row(s)) * B_tts.slice(s).col(i)).t();

      B_lls.slice(s) = B_tts.slice(s).col(i);
      P_lls.slice(s) = P_tts(s).slice(i);
    }

    lnl = lnl - log(Pr_val);
  }
  
  return Rcpp::List::create(Rcpp::Named("Pr_tls") = Pr_tls, 
                            Rcpp::Named("Pr_tts") = Pr_tts,
                            Rcpp::Named("B_tlss") = B_tlss, 
                            Rcpp::Named("B_ttss") = B_ttss, 
                            Rcpp::Named("B_tts") = B_tts,
                            Rcpp::Named("P_tlss") = P_tlss,
                            Rcpp::Named("P_ttss") = P_ttss,
                            Rcpp::Named("P_tts") = P_tts,
                            Rcpp::Named("K") = K,
                            Rcpp::Named("B_tt") = B_tt,
                            Rcpp::Named("D_tt") = D_tt,
                            Rcpp::Named("loglik") = -lnl); 
}
//library(Rcpp)
//library(RcppArmadillo)
//sourceCpp("/Users/alexhubbard/Dropbox (Opendoor)/R Codes/Packages/MarkovSwitchingDCF/src/kimfilter.cpp")

// [[Rcpp::export]]
Rcpp::List kim_smoother(arma::cube B_tlss, arma::cube B_tts, arma::mat B_tt, arma::field<arma::cube> P_tlss, arma::field<arma::cube> P_tts,
                       arma::mat Pr_tls, arma::mat Pr_tts, arma::mat Ft, arma::cube Dt, arma::mat Tr_mat){
  
  int n_states = Tr_mat.n_cols;
  int t = B_tts.slice(0).n_cols - 1;
  double tol = 1E-06;
  arma::cube B_tTss(Ft.n_rows, 1, n_states*n_states, arma::fill::zeros);
  arma::cube P_tTss(Ft.n_rows, Ft.n_rows, n_states*n_states, arma::fill::zeros);
  arma::cube Pr_tTss(1, 1, n_states*n_states, arma::fill::zeros);
  arma::mat D_tt(B_tts.n_cols, Ft.n_rows, arma::fill::zeros);
  
  for(int i = t - 1; i >= 0; i--){
    int s = 0;
    for(int st = 0; st < n_states; st++){
      for(int stf = 0; stf < n_states; stf++){
        //Full information inference on unobserved component and its covariance matrix
        //B^{j,k}_{t|T} = B^{j}_{t|t} + P^{j}_{t|t} %*% t(Fm) %*% solve(P^{j,k}_{t+1|t}) %*% (B^{k}_{t+1|T} - B^{j,k}_{t+1|t})
        B_tTss.slice(s) = B_tts.slice(st).col(i) + P_tts(st).slice(i) * Ft.t() * pinv(P_tlss(s).slice(i + 1)) * (B_tts.slice(st).col(i + 1) - B_tlss.slice(s).col(i + 1));
        
        // //P^{j,k} = P^{j}_{t|t} + P^{j}_{t|t} %*% t(Fm) %*% solve(P^{j,k}_{t+1|t}) %*% (P^{k}_{t+t|T} - P^{j,k}_{t+1|t}) %*% t(P^{j}_{t|t} %*% t(Fm) %*% solve(P^{j,k}_{t+1|t}))
        arma::mat temp = P_tts(st).slice(i) * Ft.t() * pinv(P_tlss(s).slice(i + 1));
        P_tTss.slice(s) = P_tts(st).slice(i) + temp * (P_tts(st).slice(i + 1) - P_tlss(s).slice(i + 1)) * temp.t();
        
        //Full information inference on probabilities
        //Pr[S_t|Y_T]: #Pr[S_{t+1}=k|Y_T]*Pr[S_{t+1}=k|S_t=j]*Pr[S_t=j|Y_t]/Pr[S_{t+1}=k|Y_t]
        Pr_tTss.slice(s) = Pr_tts(i + 1, stf) * Tr_mat(stf, st) * Pr_tts(i, st);
        if(arma::as_scalar(Pr_tls(i + 1, stf)) < tol){
          Pr_tTss.slice(s) /= arma::as_scalar(arma::datum::inf);
        }else{
          Pr_tTss.slice(s) /= Pr_tls(i + 1, stf);
        }
        s++;
      }
    }
    
    //Collapsing terms
    for(int s = 0; s < n_states; s++){
      //Pr[S_t=d|Y_T]
      //P^{j}_{t|t} = sum_{i}(Pr[S_t|Y_T]^{i,j})
      Pr_tts(i, s) = 0;
      for(int j = n_states - 1; j >= 0; j--){
        Pr_tts(i, s) += arma::as_scalar(Pr_tTss.slice((s + 1)*n_states - j - 1));
      }

      //B^{j}_{t|T} = (sum_{i}(Pr[S_t|Y_T]^{i,j}*B^{j,k}_{t|T}))/Pr[j]
      B_tts.slice(s).col(i) = arma::zeros(B_tts.slice(s).col(i).n_rows, B_tts.slice(s).col(i).n_cols);
      for(int j = n_states - 1; j >= 0; j--){
        B_tts.slice(s).col(i) += arma::as_scalar(Pr_tTss.slice((s + 1)*n_states - j - 1)) * B_tTss.slice((s + 1)*n_states - j - 1);
      }
      if(arma::as_scalar(Pr_tts(i, s)) < tol){
        B_tts.slice(s).col(i) /= arma::as_scalar(arma::datum::inf);
      }else{
        B_tts.slice(s).col(i) /= Pr_tts(i, s);
      }

      //P^{j}_{t|T} = (sum_{i}(Pr[S_t|Y_T]^{i,j}*(B^{j}_{t|T} - B^{j,k}_{t|T}) %*% t(B^{j}_{t|T} - B^{j,k}_{t|T})))/Pr[j]
      P_tts(s).slice(i) = arma::zeros(P_tts(s).slice(i).n_rows, P_tts(s).slice(i).n_cols);
      for(int j = n_states - 1; j >= 0; j--){
        arma::mat temp = (B_tts.slice(s).col(i) - B_tTss.slice((s + 1)*n_states - j - 1));
        P_tts(s).slice(i) += arma::as_scalar(Pr_tTss.slice((s + 1)*n_states - j - 1))*(P_tTss.slice((s + 1)*n_states - j - 1) + temp * temp.t());
      }
      if(arma::as_scalar(Pr_tts(i, s)) < tol){
        P_tts(s).slice(i) /= arma::as_scalar(arma::datum::inf);
      }else{
        P_tts(s).slice(i) /= Pr_tts(i, s);
      }

      B_tt.row(i) += (arma::as_scalar(Pr_tts(i, s)) * B_tts.slice(s).col(i)).t();
      D_tt.row(i) += (arma::as_scalar(Pr_tts(i, s)) * Dt.slice(s)).t();
    }
  }
  
  return Rcpp::List::create(Rcpp::Named("Pr_tts") = Pr_tts,
                            Rcpp::Named("P_tts") = P_tts,
                            Rcpp::Named("B_tts") = B_tts,
                            Rcpp::Named("B_tt") = B_tt,
                            Rcpp::Named("D_tt") = D_tt);
}

// B_TT = B_TL = matrix(0, ncol = ncol(yt), nrow = nrow(sp$Tt))
// P_TT = P_TL = as.list(NA, ncol(yti))
//
// #Initialize the filter
// B_LL = sp$a0
// P_LL = sp$P0
//
  // for(j in 1:ncol(yt)){
    //   ################## Kalman filter routine ##################
    //   B_TL[, j] = sp$dt + sp$Tt %*% B_LL  #Initial estimate of unobserved values conditional on t-1
    //   P_TL[[j]] = sp$Tt %*% P_LL %*% t(sp$Tt) + sp$HHt #Initial estimate of the covariance matrix conditional on t-1
    //   N_TL = yt[, j] - sp$ct - sp$Zt %*% B_TL[, j] #Prediction error conditoinal on t-1
    //   F_TL = sp$Zt %*% P_TL[[j]] %*% t(sp$Zt) + sp$GGt #Variance of the predictoin error conditional on t-1
    //   K_T = P_TL[[j]] %*% t(sp$Zt) %*% ginv(F_TL) #Kalman gain conditional on t-1
    //   B_TT[, j] = B_TL[, j] + K_T %*% N_TL #Final estimate of the unobserved values
    //   P_TT[[j]] = P_TL[[j]] - K_T %*% sp$Zt %*% P_TL[[j]] #Final estiamte of the covariance matrix
    //
      //   #Reinitialize for the next iteration
      //   B_LL = B_TT[, j]
      //   P_LL = P_TT[[j]]
      // }
//
  // #Kalman Smoother
  // for(j in (ncol(yt) - 1):1){
    //   B_TT[, j] = B_TT[, j] + P_TT[[j]] %*% t(sp$Tt) %*% ginv(P_TL[[j + 1]]) %*% (B_TT[, j + 1] - B_TL[, j + 1])
    //   P_TT[[j]] = P_TT[[j]] + P_TT[[j]] %*% t(sp$Tt) %*% ginv(P_TL[[j + 1]]) %*% (P_TT[[j + 1]] - P_TL[[j + 1]]) %*% t(P_TT[[j]] %*% t(sp$Tt) %*% ginv(P_TL[[j + 1]]))
    // }

//RcppArmadillo.package.skeleton(name = "kfdecomp", path = "/Users/alexhubbard/Dropbox (Opendoor)/R Codes/Packages")
//compileAttributes(verbose=TRUE)
//library(tools)
//package_native_routine_registration_skeleton("/Users/alexhubbard/Dropbox (Opendoor)/R Codes/Packages/hamiltonfilter")
//git config remote.origin.url git@github.com:opendoor-labs/hamiltonfilter.git

