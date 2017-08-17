// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;

const double log2pi = std::log(2.0 * M_PI);


//' Multivariate Normal Density
//'
//' This function computes the density of a multivariate normal.
//'
//' @param x a nxp matrix whose columns are the components of the multivariate normal.
//' @param mean a p-dimensional mean vector
//' @param sigma a pxp covariance matrix
//' @param logd if \code{TRUE} log densities are returned
//' @details This function is taken from \url{http://gallery.rcpp.org/articles/dmvnorm_arma/}
//' @return A n-dimensional column vector with the values of the density.
//' @author Nino Hardt, Dicko Ahmadou
//' @export
// [[Rcpp::export]]
arma::vec dmvnrmArma(arma::mat &x,
                     arma::rowvec &mean,
                     arma::mat &sigma,
                     bool logd = false,
                     bool diag = false) {
  int n = x.n_rows;
  int xdim = x.n_cols;
  arma::vec out(n);
  arma::mat rooti(n,n);
  double rootisum;

  if(diag == false){
    rooti = arma::trans(arma::inv(trimatu(arma::chol(sigma))));
    rootisum = arma::sum(log(rooti.diag()));
  }else{
    arma::vec aux = pow(diagvec(sigma), -0.5);
    rooti = diagmat(aux);
    rootisum = arma::sum(log(aux));
  }

  double constants = -(static_cast<double>(xdim)/2.0) * log2pi;

  for (int i=0; i < n; i++) {
    arma::vec z = rooti * arma::trans( x.row(i) - mean) ;
    out(i)      = constants - 0.5 * arma::sum(z%z) + rootisum;
  }

  if (logd == false) {
    out = exp(out);
  }
  return(out);
}


arma::vec dmvnrmArma(arma::rowvec &x,
                     arma::rowvec &mean,
                     arma::mat &sigma,
                     bool logd = false,
                     bool diag = false) {
  arma::mat xMatrix(1,x.size());
  xMatrix.row(0) = x;
  return(dmvnrmArma(xMatrix,mean,sigma, logd, diag));
}

// [[Rcpp::export]]
double dnrmArma(double x,
                double mean,
                double sigma,
                bool logd = false) {
  arma::mat aux(1,1);
  arma::rowvec aux2(1,1);
  arma::mat aux3(1,1);
  double res;
  aux(0,0) = x;
  aux2(0,0) = mean;
  aux3(0,0) = pow(sigma,2);
  res = dmvnrmArma(aux,aux2,aux3,logd,true).at(0);
  return(res);
}


//' Simulate from a Multivariate Normal Density
//'
//' This function produces samples from a Multivariate Normal distribution.
//'
//' @param n size of sample required
//' @param mean a p-dimensional mean vector
//' @param sigma a pxp covariance matrix
//' @details This function is taken from \url{http://gallery.rcpp.org/articles/simulate-multivariate-normal/}
//' @return A nxp matrix with the generated sample
//' @author Dicko Ahmadou
//' @export
// [[Rcpp::export]]
arma::mat mvrnormArma(int n, arma::rowvec &mean, arma::mat &sigma, bool diag = false) {
  int ncols = sigma.n_cols;
  arma::mat Y = arma::randn(n, ncols);
  arma::vec diagonal = diagvec(sigma);
  if(diag == true){
    if(n ==1){
      arma::mat out(1,ncols);
      out.row(0) = mean + Y.row(0)%sqrt(diagonal.t());
      return(out);
    }else{
      return arma::repmat(mean.t(), 1, n).t() + Y * arma::diagmat(sqrt(diagonal));
    }

  }
  return arma::repmat(mean.t(), 1, n).t() + Y * arma::chol(sigma);
}

// [[Rcpp::export]]
arma::colvec rnormArma(int n, double mean, double sigma) {
  arma::rowvec aux(1);
  arma::mat aux2(1,1);
  aux(0) = mean;
  aux2(0,0) = pow(sigma,2);
  arma::mat res = mvrnormArma(n,aux,aux2,true);
  return(res.col(0));
}
