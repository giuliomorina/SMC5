// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;
#include "smc.h"

// [[Rcpp::export]]
List KalmanFilterCpp(arma::vec& m_0, arma::mat& C_0, arma::mat& F_matrix, arma::mat& G, arma::mat& V, arma::mat& W, arma::mat& y){
  int n = y.n_rows;
  int d = m_0.n_elem;
  arma::mat m(d, n+1, fill::zeros);
  arma::cube C(d, d, n+1, fill::zeros);
  arma::vec a(d);
  arma::vec f(d);
  arma::vec m_current(d);
  arma::mat R(d,d);
  arma::mat Q(d,d);
  arma::mat C_current(d,d);
  arma::mat aux(d,d);
  m_current = m_0;
  C_current = C_0;


  m.col(0) = m_0;
  C.slice(0) = C_0;

  for (int i=0; i <n; i++){
    a = G*m_current;
    R = G*C_current*G + W;
    f = F_matrix*a;
    Q = F_matrix*R*F_matrix + V;
    aux = R*trans(F_matrix)*inv(Q);
    m_current = a + aux*(trans(y.row(i)) - F_matrix*a); // change here
    C_current = R - aux*F_matrix*R;
    m.col(i+1) = m_current;
    C.slice(i+1) = C_current;
  }
  return(List::create(Named("m") = m,
                      Named("C") = C));
}

// [[Rcpp::export]]
List KalmanSmoothingCpp(arma::colvec& m_0, arma::mat& C_0, arma::mat& F_matrix,
                        arma::mat& G, arma::mat& V, arma::mat& W, arma::mat& y) {
  //Setup variables
  int d = m_0.n_elem;
  int n = y.n_rows+1;
  arma::cube Kalman_smoothing_cov(d, d, n);
  arma::mat Kalman_smoothing_mean(d, n);

  //Perform Kalman Filter
  List Kalman_filter_res = KalmanFilterCpp(m_0,C_0,F_matrix,G,V,W,y);
  arma::mat Kalman_filter_mean = Kalman_filter_res["m"];
  arma::cube Kalman_filter_cov = Kalman_filter_res["C"];

  //Kalman smoothing at final time corresponds to filtering
  Kalman_smoothing_mean.col(n-1) = Kalman_filter_mean.col(n-1);
  Kalman_smoothing_cov.slice(n-1) = Kalman_filter_cov.slice(n-1);

  //Backward computation of the smoothing
  //This is based on this paper:"Movellan J. R. (2011) Discrete Time Kalman Filters and Smoothers. MPLab Tu- torials, University of California San Diego"
  arma::mat sigma_t_t1(C_0.n_rows, C_0.n_cols);
  arma::colvec mu_t_t1(m_0.n_elem);
  arma::mat g_t(C_0.n_rows, C_0.n_cols);
  for(int t = n-2; t>=0; t--) {
    sigma_t_t1 = G*Kalman_filter_cov.slice(t)*G.t() + W;
    mu_t_t1 = G*Kalman_filter_mean.col(t);
    g_t = arma::solve(sigma_t_t1, Kalman_filter_cov.slice(t)*G.t());
    Kalman_smoothing_mean.col(t) = Kalman_filter_mean.col(t) + g_t*(Kalman_smoothing_mean.col(t+1)-mu_t_t1);
    Kalman_smoothing_cov.slice(t) = Kalman_filter_cov.slice(t) - g_t*(sigma_t_t1-Kalman_smoothing_cov.slice(t+1))*g_t.t();
  }

  return(List::create(Named("m") = Kalman_smoothing_mean,
                      Named("C") = Kalman_smoothing_cov));
}

// [[Rcpp::export]]
arma::mat KalmanSamplerCpp(int n, List kalmanRes, int time) {
  arma::mat m = kalmanRes["m"];
  arma::cube C = kalmanRes["C"];
  arma::rowvec kalmanMean = m.col(time).t();
  arma::mat kalmanCovMatrix = C.slice(time);
  return(mvrnormArma(n,kalmanMean,kalmanCovMatrix,checkDiagonal(kalmanCovMatrix)));
}

